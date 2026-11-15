#include "core/utilities/parallel.hpp"

#include "gtest/gtest.h"

#include <atomic>
#include <chrono>
#include <thread>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <stdexcept>

#if defined(__linux__)
#include <sched.h>
#include <string>
#include <fstream>
#include <unistd.h>
#endif

using namespace PHARE::core;

namespace
{
// restores the serial default so tests do not depend on their order
struct Threads
{
    explicit Threads(std::size_t const n) { setNbrThreads(n); }
    ~Threads() { setNbrThreads(1); }
};
} // namespace


TEST(ParallelFor, isSerialByDefault)
{
    EXPECT_EQ(nbrThreads(), 1u);

    std::vector<std::size_t> order;
    parallel_for(5, [&](std::size_t const i) {
        EXPECT_EQ(threadIndex(), 0u);
        order.push_back(i);
    });
    EXPECT_EQ(order, (std::vector<std::size_t>{0, 1, 2, 3, 4}));
}


TEST(ParallelFor, visitsEveryIndexOnce)
{
    Threads const threads{4};
    std::size_t const n = 1000;

    std::vector<int> visits(n, 0);
    parallel_for(n, [&](std::size_t const i) { ++visits[i]; });

    for (auto const v : visits)
        EXPECT_EQ(v, 1);
}


TEST(ParallelFor, threadIndexIsBelowTheNumberOfThreads)
{
    Threads const threads{4};

    for (std::size_t const n : {std::size_t{2}, std::size_t{100}})
    {
        std::atomic<std::size_t> outOfRange{0};
        parallel_for(n, [&](std::size_t) {
            if (threadIndex() >= nbrThreads())
                ++outOfRange;
        });
        EXPECT_EQ(outOfRange, 0u);
    }
}


// loops over the same items give an item to the same thread, whose cache holds its data
TEST(ParallelFor, outermostLoopGivesEachThreadItsOwnIndicesFirst)
{
    std::size_t const nThreads = 4, perThread = 3;
    Threads const threads{nThreads};

    for (int repeat = 0; repeat < 3; ++repeat)
    {
        std::vector<std::size_t> threadOf(nThreads * perThread, nThreads);
        std::vector<char> counted(nThreads, 0);
        std::atomic<std::size_t> arrived{0};

        parallel_for(threadOf.size(), [&](std::size_t const i) {
            // nobody goes on to the indices of another before every thread has one of its own
            if (!counted[threadIndex()])
            {
                counted[threadIndex()] = 1;
                ++arrived;
            }
            while (arrived < nThreads)
                std::this_thread::yield();
            threadOf[i] = threadIndex();
        });

        // a thread that is done early takes the last indices of another: the first ones stay
        for (std::size_t t = 0; t < nThreads; ++t)
            EXPECT_EQ(threadOf[t * perThread], t) << "repeat " << repeat;
    }
}


TEST(ParallelFor, nestedCallVisitsEveryIndexOnce)
{
    Threads const threads{4};
    std::size_t const nOuter = 3, nInner = 50;

    std::vector<int> visits(nOuter * nInner, 0);
    std::atomic<std::size_t> outOfRange{0};
    parallel_for(nOuter, [&](std::size_t const i) {
        parallel_for(nInner, [&](std::size_t const j) {
            ++visits[i * nInner + j];
            if (threadIndex() >= nbrThreads())
                ++outOfRange;
        });
    });

    for (auto const v : visits)
        EXPECT_EQ(v, 1);
    EXPECT_EQ(outOfRange, 0u);
}


TEST(ParallelFor, nestedCallIsSharedWithTheThreadsLeftWithoutAnOuterIndex)
{
    Threads const threads{4};

    // one outer index: three threads have nothing else to do than the inner indices
    std::vector<std::size_t> ranOn(32);
    parallel_for(1, [&](std::size_t) {
        parallel_for(ranOn.size(), [&](std::size_t const j) {
            ranOn[j] = threadIndex();
            volatile double x = 0;
            for (int k = 0; k < 2000000; ++k)
                x = x + k;
        });
    });

    std::vector<int> used(nbrThreads(), 0);
    for (auto const t : ranOn)
        used[t] = 1;
    EXPECT_GT(std::count(used.begin(), used.end(), 1), 1);
}


TEST(ParallelFor, manyNestedLoopsLoseNoIndex)
{
    Threads const threads{8};
    std::size_t const nOuter = 37, nInner = 23, nRepeats = 200;

    for (std::size_t r = 0; r < nRepeats; ++r)
    {
        std::atomic<std::size_t> sum{0};
        parallel_for(nOuter, [&](std::size_t const i) {
            parallel_for(nInner, [&](std::size_t const j) { sum += i * nInner + j + 1; });
            parallel_for(2, [&](std::size_t) { ++sum; });
        });
        auto const n = nOuter * nInner;
        ASSERT_EQ(sum, n * (n + 1) / 2 + 2 * nOuter);
    }
}


TEST(ParallelFor, workersThatWentToSleepComeBack)
{
    Threads const threads{4};

    std::this_thread::sleep_for(std::chrono::milliseconds{20}); // longer than they spin
    std::vector<int> visits(1000, 0);
    parallel_for(visits.size(), [&](std::size_t const i) { ++visits[i]; });
    for (auto const v : visits)
        EXPECT_EQ(v, 1);
}


TEST(ParallelFor, nestedCallRethrows)
{
    Threads const threads{4};

    EXPECT_THROW(parallel_for(2,
                              [&](std::size_t) {
                                  parallel_for(16, [&](std::size_t const j) {
                                      if (j == 5)
                                          throw std::runtime_error("inner 5");
                                  });
                              }),
                 std::runtime_error);
}


#if defined(__linux__)
namespace
{
// first cpu of the core a cpu belongs to
int coreOf(int const cpu)
{
    int first = cpu;
    std::ifstream{"/sys/devices/system/cpu/cpu" + std::to_string(cpu)
                  + "/topology/thread_siblings_list"}
        >> first;
    return first;
}

std::vector<int> cpusOf(cpu_set_t const& set)
{
    std::vector<int> cpus;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        if (CPU_ISSET(cpu, &set))
            cpus.push_back(cpu);
    return cpus;
}
} // namespace


TEST(ThreadPool, pinsEachThreadToItsOwnCoreWhenTheProcessIsRestricted)
{
    cpu_set_t launch;
    ASSERT_EQ(sched_getaffinity(0, sizeof(cpu_set_t), &launch), 0);

    // restrict this process to every cpu of three cores, as a launcher would
    std::vector<int> cores;
    cpu_set_t restricted;
    CPU_ZERO(&restricted);
    for (auto const cpu : cpusOf(launch))
    {
        auto const core = coreOf(cpu);
        if (std::find(cores.begin(), cores.end(), core) == cores.end() and cores.size() < 3)
            cores.push_back(core);
        if (std::find(cores.begin(), cores.end(), core) != cores.end())
            CPU_SET(cpu, &restricted);
    }
    if (cores.size() < 3 or CPU_COUNT(&restricted) == CPU_COUNT(&launch))
        GTEST_SKIP() << "needs more than three cores";

    ASSERT_EQ(sched_setaffinity(0, sizeof(cpu_set_t), &restricted), 0);
    {
        Threads const threads{3};

        std::vector<int> coreOfThread(3, -1), cpusOfThread(3, 0);
        for (int repeat = 0; repeat < 20; ++repeat)
            parallel_for(64, [&](std::size_t) {
                cpu_set_t mine;
                sched_getaffinity(0, sizeof(cpu_set_t), &mine);
                cpusOfThread[threadIndex()] = CPU_COUNT(&mine);
                coreOfThread[threadIndex()] = coreOf(sched_getcpu());
                volatile double x = 0;
                for (int k = 0; k < 20000; ++k)
                    x = x + k;
            });

        std::sort(coreOfThread.begin(), coreOfThread.end());
        EXPECT_EQ(coreOfThread, cores) << "one thread per core";
        EXPECT_EQ(cpusOfThread, (std::vector<int>{1, 1, 1}));
    }

    // back to one thread: the caller may run where the process was launched again
    cpu_set_t after;
    sched_getaffinity(0, sizeof(cpu_set_t), &after);
    EXPECT_TRUE(CPU_EQUAL(&after, &restricted));
    sched_setaffinity(0, sizeof(cpu_set_t), &launch);
}


TEST(ThreadPool, leavesThreadsFreeWhenTheProcessMayRunAnywhere)
{
    cpu_set_t launch;
    ASSERT_EQ(sched_getaffinity(0, sizeof(cpu_set_t), &launch), 0);
    if (CPU_COUNT(&launch) < sysconf(_SC_NPROCESSORS_ONLN))
        GTEST_SKIP() << "this process is restricted already";

    Threads const threads{2};
    std::atomic<int> pinned{0};
    parallel_for(64, [&](std::size_t) {
        cpu_set_t mine;
        sched_getaffinity(0, sizeof(cpu_set_t), &mine);
        if (CPU_COUNT(&mine) != CPU_COUNT(&launch))
            ++pinned;
    });
    EXPECT_EQ(pinned, 0);
}
#endif


TEST(Chunks, coverTheRangeOnceInOrder)
{
    Threads const threads{4};

    EXPECT_EQ(nbrChunks(10, 100), 1u);
    EXPECT_EQ(nbrChunks(1000, 100), 10u);

    std::size_t const n = 1003, nChunks = nbrChunks(n, 100);
    std::size_t expectedFirst = 0;
    for (std::size_t c = 0; c < nChunks; ++c)
    {
        auto const [first, last] = chunk(n, nChunks, c);
        EXPECT_EQ(first, expectedFirst);
        EXPECT_GE(last - first, 100u);
        expectedFirst = last;
    }
    EXPECT_EQ(expectedFirst, n);
}


TEST(ParallelForPieces, coverTheRangeOnceWithConsecutivePieces)
{
    Threads const threads{4};

    std::size_t const n = 10000, grain = 1000;
    std::vector<int> visits(n, 0);
    std::atomic<std::size_t> nPieces{0};

    parallel_for_pieces(n, grain, [&](std::size_t const first, std::size_t const last) {
        EXPECT_LT(first, last);
        EXPECT_GE(last - first, grain);
        ++nPieces;
        for (auto i = first; i < last; ++i)
            ++visits[i];
    });

    EXPECT_EQ(nPieces, 10u);
    EXPECT_TRUE(std::all_of(visits.begin(), visits.end(), [](int const v) { return v == 1; }));
}

TEST(ParallelForPieces, isOneCallWithOneThreadOrFewItems)
{
    std::size_t calls = 0;
    parallel_for_pieces(10000, 1000, [&](std::size_t const first, std::size_t const last) {
        ++calls;
        EXPECT_EQ(first, 0u);
        EXPECT_EQ(last, 10000u);
    });
    EXPECT_EQ(calls, 1u);

    Threads const threads{4};
    calls = 0;
    parallel_for_pieces(999, 1000, [&](std::size_t, std::size_t) { ++calls; });
    EXPECT_EQ(calls, 1u);
}

TEST(Chunks, areOneWithOneThread)
{
    EXPECT_EQ(nbrChunks(1000000, 100), 1u);
}


TEST(ParallelFor, rethrowsAnExceptionAfterAllIndicesRan)
{
    Threads const threads{4};
    std::size_t const n = 64;

    std::atomic<std::size_t> ran{0};
    EXPECT_THROW(parallel_for(n,
                              [&](std::size_t const i) {
                                  ++ran;
                                  if (i == 7)
                                      throw std::runtime_error("index 7");
                              }),
                 std::runtime_error);

    if (nbrThreads() > 1) // serial stops at the throw, as a plain loop does
    {
        EXPECT_EQ(ran, n);
    }
}


TEST(ParallelFor, emptyRangeCallsNothing)
{
    Threads const threads{4};
    parallel_for(0, [](std::size_t) { FAIL(); });
}


TEST(ParallelFor, threadsComeFromTheEnvironment)
{
    setenv("PHARE_THREADS", "3", 1);
    setNbrThreadsFromEnv();
    EXPECT_EQ(nbrThreads(), 3u);

    unsetenv("PHARE_THREADS");
    setNbrThreadsFromEnv();
    EXPECT_EQ(nbrThreads(), 1u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
