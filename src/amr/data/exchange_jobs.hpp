#ifndef PHARE_AMR_DATA_EXCHANGE_JOBS_HPP
#define PHARE_AMR_DATA_EXCHANGE_JOBS_HPP

#include "core/utilities/parallel.hpp"

#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/ScheduleOpsStrategy.h>

#include <new>
#include <vector>
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
      another destination writes. Ghost fills write ghosts and read interiors, border sums read
      a variable they do not write.
  With one thread nothing is recorded and every call runs at once, as without this class.

  Attached to the schedules where they are stored (Communicator::add).
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
    // nullptr: the job shares its destination with no other job
    void add(void const* const destination, Job const& job)
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
                groups_.emplace_back();
        }
        groups_[group].push_back(job);
    }

    // runs the recorded jobs now: for a call that must not overtake them
    void settle()
    {
        if (recording() and nGroups_ > 0)
            run_();
    }

    void preCommunicate() override {}
    void postCommunicate() override {}

    void preCopy() override { start_(); }
    void postCopy() override { stop_(); }

    void prePack() override { start_(); }
    void postPack() override { stop_(); }

    void preUnpack() override { start_(); }
    void postUnpack() override { stop_(); }

    // a recorded pack fills its part of the message after SAMRAI's walk
    bool deferMessageSend() override { return core::nbrThreads() > 1; }

private:
    void start_() { recording_ = core::nbrThreads() > 1; }

    void stop_()
    {
        if (nGroups_ > 0)
            run_();
        recording_ = false;
    }

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
        core::parallel_for(nGroups_, [&](std::size_t const g) {
            for (auto const& job : groups_[g])
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
    std::unordered_map<void const*, Group> groupOf_;
};



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
