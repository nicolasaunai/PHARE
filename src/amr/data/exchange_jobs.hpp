#ifndef PHARE_AMR_DATA_EXCHANGE_JOBS_HPP
#define PHARE_AMR_DATA_EXCHANGE_JOBS_HPP

#include "core/utilities/parallel.hpp"

#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/ScheduleOpsStrategy.h>

#include <new>
#include <memory>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

namespace PHARE::amr
{
/*
  Runs the copies, packs and unpacks of a SAMRAI schedule on the threads of the process.

  A schedule walks its transactions one by one, on one thread, in three phases: pack what goes
  to other processes, copy between local patches, unpack what was received. SAMRAI lets patch
  data postpone the work of a phase: the calls it makes through PatchData::copyFuseable,
  packStreamFuseable and unpackStreamFuseable only have to be done when the phase ends, which is
  when it calls the matching post hook below (tbox::Schedule: postPack, postCopy, postUnpack).
  Messages are sent after postPack (deferMessageSend) and receive buffers are released after
  postUnpack, so a postponed pack or unpack still finds its buffer.

  With more than one thread, patch data therefore only record a job per call, and the post hook
  runs them with core::parallel_for:
    - jobs writing the same destination run on one thread, in the order SAMRAI made the calls,
      so "the last writer wins" holds as in the serial walk;
    - jobs of different destinations run concurrently: a job must read nothing that a job of
      another destination writes. Ghost fills write ghosts and read interiors, border sums and
      maxima read a variable they do not write.
  With one thread nothing is recorded and every call runs at once, as without this class.

  Attached to the schedules where they are stored (Communicator::add).

  The same is done for the work SAMRAI does patch after patch between two levels: interpolating
  a fine level from the next coarser one, and coarsening a fine level. There the phase is opened
  and closed by begin() and end(), called by the patch strategies of
  concurrent_level_strategies.hpp, and the jobs are the calls to PHARE's refine and coarsen
  operators, grouped by the patch they write.
*/
class ExchangeJobs : public SAMRAI::tbox::ScheduleOpsStrategy
{
public:
    // What a job runs, kept in the job itself: thousands are recorded and dropped per time step
    // by the one thread that walks the schedules, which must not allocate for each.
    // fn holds pointers and sizes only, to data that outlives the phase.
    class Job
    {
    public:
        template<typename Fn>
        Job(Fn const& fn)
            : run_{[](void const* const stored) { (*static_cast<Fn const*>(stored))(); }}
        {
            static_assert(sizeof(Fn) <= sizeof(stored_) and alignof(Fn) <= alignof(void*)
                          and std::is_trivially_copyable_v<Fn>
                          and std::is_trivially_destructible_v<Fn>);
            ::new (static_cast<void*>(stored_)) Fn(fn);
        }

        void operator()() const { run_(stored_); }

    private:
        void (*run_)(void const*);
        alignas(void*) unsigned char stored_[48];
    };

    static ExchangeJobs& instance()
    {
        static ExchangeJobs jobs;
        return jobs;
    }

    // true while SAMRAI walks the transactions of a phase and jobs are recorded
    bool recording() const { return recording_ and !running_; }

    // job runs when the phase ends, after the jobs added before it for the same destination.
    // nullptr: the job shares its destination with no other job.
    // Jobs given the place of their patch in its level run in that order, so that the threads
    // share them out as they share out a loop over the patches, see core::ThreadPool
    void add(void const* const destination, Job const& job, std::size_t const place = noPlace)
    {
        auto group = nGroups_;
        if (destination)
        {
            auto& known = groupOf_[destination];
            if (known.phase != phase_)
                known = {nGroups_, phase_};
            group = known.group;
        }
        if (group == nGroups_)
        {
            if (++nGroups_ > groups_.size())
            {
                groups_.emplace_back();
                places_.emplace_back();
            }
            places_[group] = place;
            placed_        = placed_ or place != noPlace;
        }
        groups_[group].push_back(job);
    }

    static constexpr std::size_t noPlace = static_cast<std::size_t>(-1);

    // a copy of value that lives until the jobs of the phase have run: for what a job needs
    // and its caller only lends, such as a local variable of a SAMRAI loop
    template<typename T>
    T const& keep(T const& value)
    {
        auto kept = std::make_shared<T>(value);
        kept_.push_back(kept);
        return *kept;
    }

    // opens a phase, in which jobs are recorded if there are several threads
    void begin() { recording_ = core::nbrThreads() > 1; }

    // closes the phase: its jobs have run on return
    void end()
    {
        if (nGroups_ > 0)
            run_();
        recording_ = false;
    }

    // runs the recorded jobs now: for a call that must not overtake them
    void settle()
    {
        if (recording() and nGroups_ > 0)
            run_();
    }

    void preCommunicate() override {}
    void postCommunicate() override {}

    void preCopy() override { begin(); }
    void postCopy() override { end(); }

    void prePack() override { begin(); }
    void postPack() override { end(); }

    void preUnpack() override { begin(); }
    void postUnpack() override { end(); }

    // a recorded pack fills its part of the message after SAMRAI's walk
    bool deferMessageSend() override { return core::nbrThreads() > 1; }

private:
    void run_()
    {
        struct Running
        {
            ExchangeJobs& self;
            ~Running()
            {
                for (std::size_t g = 0; g < self.nGroups_; ++g)
                    self.groups_[g].clear();
                self.nGroups_ = 0;
                self.placed_  = false;
                self.kept_.clear();
                ++self.phase_;
                // destinations stay known from one phase to the next, patch data of past
                // hierarchies included: forget them all from time to time
                if (self.groupOf_.size() > maxKnownDestinations_)
                    self.groupOf_.clear();
                self.running_ = false;
            }
        };

        running_ = true;
        Running const running{*this};

        order_.resize(nGroups_);
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        if (placed_)
            std::stable_sort(order_.begin(), order_.end(),
                             [&](auto const a, auto const b) { return places_[a] < places_[b]; });

        core::parallel_for(nGroups_, [&](std::size_t const i) {
            for (auto const& job : groups_[order_[i]])
                job();
        });
    }

    struct Group
    {
        std::size_t group = 0;
        std::size_t phase = 0; // the group is of that phase only
    };

    static constexpr std::size_t maxKnownDestinations_ = 1 << 20;

    bool recording_      = false;
    bool running_        = false;
    std::size_t nGroups_ = 0;
    std::size_t phase_   = 1;
    std::vector<std::vector<Job>> groups_; // kept between phases, only the first nGroups_ are used
    std::vector<std::size_t> places_;      // of the groups, see add
    std::vector<std::size_t> order_;       // the groups in the order they are shared out
    bool placed_ = false;                  // some group has a place
    std::unordered_map<void const*, Group> groupOf_;
    std::vector<std::shared_ptr<void>> kept_;
};



// value itself, or a copy of it kept for the jobs being recorded
template<typename T>
T const& forJobs(T const& value)
{
    auto& jobs = ExchangeJobs::instance();
    return jobs.recording() ? jobs.keep(value) : value;
}


// copy() writes destination: now, or when the phase SAMRAI is walking ends
void copyOrDefer(void const* const destination, auto const& copy)
{
    auto& jobs = ExchangeJobs::instance();
    if (jobs.recording())
        jobs.add(destination, copy);
    else
        copy();
}


// pack(stream) appends exactly nbytes to its stream: now, or when the phase ends, into a part of
// the message reserved now so that the message keeps the order of SAMRAI's calls
void packOrDefer(SAMRAI::tbox::MessageStream& stream, std::size_t const nbytes, auto const& pack)
{
    auto& jobs = ExchangeJobs::instance();

    if (!jobs.recording() or nbytes == 0)
        return pack(stream);

    if (!stream.canCopyIn(nbytes))
    {
        // the stream is smaller than announced: growing it moves it, with the recorded packs
        // still pointing into it
        jobs.settle();
        return pack(stream);
    }

    auto* const reserved = stream.getWriteBuffer<char>(nbytes);
    jobs.add(nullptr, [reserved, nbytes, pack]() {
        SAMRAI::tbox::MessageStream packed{nbytes, SAMRAI::tbox::MessageStream::Write};
        pack(packed);
        if (packed.getCurrentSize() != nbytes)
            throw std::runtime_error("ExchangeJobs: a pack did not write the announced size");
        std::memcpy(reserved, packed.getBufferStart(), nbytes);
    });
}


// unpack(stream) reads exactly nbytes from its stream and writes destination: now, or when the
// phase ends, from the part of the message it would have read now
void unpackOrDefer(void const* const destination, SAMRAI::tbox::MessageStream& stream,
                   std::size_t const nbytes, auto const& unpack)
{
    auto& jobs = ExchangeJobs::instance();

    if (!jobs.recording() or nbytes == 0)
        return unpack(stream);

    auto const* const received = stream.getReadBuffer<char>(nbytes);
    jobs.add(destination, [received, nbytes, unpack]() {
        SAMRAI::tbox::MessageStream part{nbytes, SAMRAI::tbox::MessageStream::Read, received,
                                         false /* read in place */};
        unpack(part);
    });
}

} // namespace PHARE::amr

#endif // PHARE_AMR_DATA_EXCHANGE_JOBS_HPP
