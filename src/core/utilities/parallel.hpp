#ifndef PHARE_CORE_UTILITIES_PARALLEL_HPP
#define PHARE_CORE_UTILITIES_PARALLEL_HPP

#include <mutex>
#include <atomic>
#include <chrono>
#include <string>
#include <thread>
#include <vector>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <algorithm>
#include <exception>
#include <condition_variable>

#if defined(__linux__)
#include <sched.h>
#include <fstream>
#include <unistd.h>
#endif

// Threads of one MPI process: parallel_for at the end of this file is the interface.

namespace PHARE::core
{
/*
  The threads parallel_for runs on: the thread that calls it from outside any loop (the main
  thread, index 0) and size() - 1 workers started by resize.

  A loop is a count of indices and a function. The thread that starts one publishes it, takes
  indices from it until none is left, then waits for the threads that took some to be done. A
  worker with nothing to do joins a published loop and takes indices from it until none is left.

  Loops nest: the function of a loop may start loops of its own, typically the pieces of the work
  of one patch. A worker joins the outermost loop that has indices left, so that the pieces of a
  patch are only shared once every patch has a thread. A thread waiting for its loop helps the
  loops nested at least as deep in the meantime, never an outer one: it would start another patch
  while its own is unfinished.

  Idle workers spin for a short while, then sleep until a loop is published.

  When the launcher (mpirun, srun) restricted the process to a part of the machine with at least
  as many cores as threads, each thread is pinned to one of these cores (Linux): two hardware
  threads of one core are given to the process as two cpus, and threads left free to move end up
  sharing a core while another is idle. A process that may run anywhere is left alone: several
  of them would pin their threads to the same cores.
*/
class ThreadPool
{
public:
    static ThreadPool& instance()
    {
        static ThreadPool pool;
        return pool;
    }

    ThreadPool()                             = default;
    ThreadPool(ThreadPool const&)            = delete;
    ThreadPool& operator=(ThreadPool const&) = delete;
    ~ThreadPool() { resize(1); }

    // threads loops run on, the calling one included
    std::size_t size() const { return workers_.size() + 1; }

    // from outside any loop only
    void resize(std::size_t const nThreads)
    {
        if (std::max<std::size_t>(nThreads, 1) == size())
            return;

        {
            std::lock_guard const lock{mutex_};
            stop_ = true;
        }
        wakeUp_.notify_all();
        for (auto& worker : workers_)
            worker.join();
        workers_.clear();
        stop_ = false;

        auto const cores = coresFor_(nThreads);
        pinTo_(cores, 0);
        for (std::size_t index = 1; index < nThreads; ++index)
            workers_.emplace_back([this, index, cores] {
                pinTo_(cores, index);
                work_(index);
            });
    }

    // index of the calling thread among the threads of the pool
    static std::size_t threadIndex() { return index_; }

    // fn(i) has run for i in [0, n) on return, the first exception it threw being rethrown
    template<typename Fn>
    void run(std::size_t const n, Fn& fn)
    {
        Loop loop{n, depth_, [](void* const f, std::size_t const i) { (*static_cast<Fn*>(f))(i); },
                  const_cast<void*>(static_cast<void const*>(&fn))};

        if (n == 1) // nothing to share
            takeIndices_(loop);
        else
        {
            {
                std::lock_guard const lock{mutex_};
                published_.push_back(&loop);
                ++unfinished_;
            }
            if (sleepers_ > 0)
                wakeUp_.notify_all();

            takeIndices_(loop);

            {
                std::lock_guard const lock{mutex_};
                published_.erase(std::find(published_.begin(), published_.end(), &loop));
            }

            // no thread can join the loop any more, those that did are running its last indices
            while (loop.helpers > 0)
                if (unfinished_ == 0 or !help_(loop.depth))
                    relax_();
        }

        if (loop.error)
            std::rethrow_exception(loop.error);
    }

private:
    struct Loop
    {
        std::size_t n;
        std::size_t depth; // loops it is nested in
        void (*call)(void*, std::size_t);
        void* fn;
        std::atomic<std::size_t> next{0};    // index to take, none left from n on
        std::atomic<std::size_t> helpers{0}; // threads that joined and may still use the loop
        std::atomic<bool> finished{false};   // someone saw that no index is left
        std::exception_ptr error{};
    };

    // a published loop is counted in unfinished_ until the first thread sees it has no index left
    void sawNoIndexLeft_(Loop& loop)
    {
        if (loop.n > 1 and !loop.finished.exchange(true))
            --unfinished_;
    }

    // runs indices of the loop until none is left
    void takeIndices_(Loop& loop)
    {
        auto const outerDepth = depth_;
        depth_                = loop.depth + 1;

        for (auto i = loop.next.fetch_add(1); i < loop.n; i = loop.next.fetch_add(1))
        {
            try
            {
                loop.call(loop.fn, i);
            }
            catch (...)
            {
                std::lock_guard const lock{mutex_};
                if (!loop.error)
                    loop.error = std::current_exception();
            }
        }

        depth_ = outerDepth;
        sawNoIndexLeft_(loop);
    }

    // joins the outermost published loop nested at least minDepth deep that has indices left
    bool help_(std::size_t const minDepth)
    {
        Loop* loop = nullptr;
        {
            std::lock_guard const lock{mutex_};
            for (auto* const candidate : published_)
            {
                if (candidate->next >= candidate->n)
                    sawNoIndexLeft_(*candidate); // its threads are busy with its last indices
                else if (candidate->depth >= minDepth and (!loop or candidate->depth < loop->depth))
                    loop = candidate;
            }
            if (loop)
                ++loop->helpers; // keeps the loop alive: its owner waits for helpers to leave
        }
        if (!loop)
            return false;

        takeIndices_(*loop);
        --loop->helpers;
        return true;
    }

    void work_(std::size_t const index)
    {
        using Clock = std::chrono::steady_clock;

        index_         = index;
        auto idleSince = Clock::now();

        for (std::size_t spins = 1; !stop_; ++spins)
        {
            if (unfinished_ > 0 and help_(0))
            {
                idleSince = Clock::now();
                continue;
            }

            relax_();

            if (spins % 1024 == 0 and Clock::now() - idleSince > spinTime_)
            {
                std::unique_lock lock{mutex_};
                ++sleepers_;
                wakeUp_.wait(lock, [&] { return stop_ or unfinished_ > 0; });
                --sleepers_;
                idleSince = Clock::now();
            }
        }
    }

#if defined(__linux__)
    // One cpu per core for nThreads threads to be pinned on, none if they are to be left free.
    // A calling thread pinned by a previous call gets back the cpus it had first.
    std::vector<int> coresFor_(std::size_t const nThreads)
    {
        if (callerCpus_)
        {
            sched_setaffinity(0, sizeof(cpu_set_t), &*callerCpus_);
            callerCpus_.reset();
        }

        cpu_set_t allowed;
        if (nThreads < 2 or sched_getaffinity(0, sizeof(cpu_set_t), &allowed) != 0)
            return {};

        std::vector<int> cores, seen;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        {
            if (!CPU_ISSET(cpu, &allowed))
                continue;

            // the first cpu of the core this one belongs to, e.g. 3 from "3,51"
            int first = cpu;
            std::ifstream{"/sys/devices/system/cpu/cpu" + std::to_string(cpu)
                          + "/topology/thread_siblings_list"}
                >> first;

            if (std::find(seen.begin(), seen.end(), first) == seen.end())
            {
                seen.push_back(first);
                cores.push_back(cpu);
            }
        }

        bool const restricted = CPU_COUNT(&allowed) < sysconf(_SC_NPROCESSORS_ONLN);
        if (!restricted or cores.size() < nThreads)
            return {};

        callerCpus_ = allowed;
        return cores;
    }

    static void pinTo_(std::vector<int> const& cores, std::size_t const index)
    {
        if (index < cores.size())
        {
            cpu_set_t cpu;
            CPU_ZERO(&cpu);
            CPU_SET(cores[index], &cpu);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpu); // 0: the calling thread
        }
    }

    std::optional<cpu_set_t> callerCpus_; // of the calling thread before it was pinned
#else
    static std::vector<int> coresFor_(std::size_t) { return {}; }
    static void pinTo_(std::vector<int> const&, std::size_t) {}
#endif

    static void relax_()
    {
#if defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#else
        std::this_thread::yield();
#endif
    }

    // the gaps between two loops of a time step are shorter than this
    static constexpr std::chrono::milliseconds spinTime_{2};

    static inline thread_local std::size_t index_ = 0;
    static inline thread_local std::size_t depth_ = 0; // loops the thread is running an index of

    std::vector<std::thread> workers_;

    std::mutex mutex_; // published_, Loop::error, and what wakeUp_ waits for
    std::condition_variable wakeUp_;
    std::vector<Loop*> published_;
    std::atomic<std::size_t> unfinished_{0}; // published loops not yet seen without index left
    std::atomic<std::size_t> sleepers_{0};
    std::atomic<bool> stop_{false};
};



// threads parallel_for may use on this process, 1 unless set
inline std::size_t nbrThreads()
{
    return ThreadPool::instance().size();
}

// from outside any parallel_for only
inline void setNbrThreads(std::size_t const n)
{
    ThreadPool::instance().resize(n);
}

// PHARE_THREADS, 1 if unset: an MPI run that does not ask for threads never starts any
inline void setNbrThreadsFromEnv()
{
    auto const* const env = std::getenv("PHARE_THREADS");
    setNbrThreads(env ? std::stoul(env) : 1);
}

// index of the calling thread among the threads of the process, below nbrThreads()
inline std::size_t threadIndex()
{
    return ThreadPool::threadIndex();
}


// fn(i) has run for every i in [0, n) on return. Indices are handed one at a time to idle
// threads, so the order is not defined: fn(i) and fn(j) must not write the same data.
//
// Called from fn itself, the inner indices are also offered to the threads that have no outer
// index left: the pieces of the last items of a loop are shared instead of leaving threads
// idle. An inner fn may therefore run on any thread, and must not use what belongs to the thread
// of its caller.
//
// The first exception thrown by fn is rethrown here, once all indices are done. With one thread
// this is a plain loop.
template<typename Fn>
void parallel_for(std::size_t const n, Fn&& fn)
{
    auto& pool = ThreadPool::instance();

    if (pool.size() > 1 and n > 0)
        pool.run(n, fn);
    else
        for (std::size_t i = 0; i < n; ++i)
            fn(i);
}


// number of pieces of at least grain items [0, n) is cut in for parallel_for, 1 with one thread
inline std::size_t nbrChunks(std::size_t const n, std::size_t const grain)
{
    return nbrThreads() > 1 ? std::max<std::size_t>(n / grain, 1) : 1;
}

// [first, last) of piece c among the nChunks consecutive pieces of [0, n)
inline auto chunk(std::size_t const n, std::size_t const nChunks, std::size_t const c)
{
    struct Bounds
    {
        std::size_t first, last;
    };
    return Bounds{c * n / nChunks, (c + 1) * n / nChunks};
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_PARALLEL_HPP
