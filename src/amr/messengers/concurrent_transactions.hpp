#ifndef PHARE_AMR_MESSENGERS_CONCURRENT_TRANSACTIONS_HPP
#define PHARE_AMR_MESSENGERS_CONCURRENT_TRANSACTIONS_HPP

#include "core/utilities/parallel.hpp"

#include "amr/data/exchange_jobs.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/BoxOverlap.h>
#include <SAMRAI/tbox/Transaction.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/hier/ComponentSelector.h>
#include <SAMRAI/tbox/TransactionFuseable.h>
#include <SAMRAI/xfer/RefineTransactionFactory.h>
#include <SAMRAI/xfer/StandardRefineTransactionFactory.h>

#include <map>
#include <tuple>
#include <memory>
#include <vector>
#include <cstddef>
#include <ostream>
#include <utility>

/*
  A refine schedule holds one transaction per pair of overlapping patches and per quantity. Each
  time it fills its level it walks them all on one thread: thousands of calls that look patch
  data up, only to record a job for the threads of the process (ExchangeJobs).

  SAMRAI asks a factory for each transaction when it builds the schedule. The factory below gives
  it, for the copies between two patches of this process, one transaction per destination patch:
  it holds SAMRAI's own transactions for that patch and runs them all as one job. The walk then
  records one job per patch, and the patch data are looked up by the threads. Time interpolated
  copies, which SAMRAI makes through a temporary that cannot be recorded, run on the threads too.

  Transactions with another process are handed over unchanged: their packs and unpacks are
  recorded one by one.

  With one thread schedules are given the factory they had, see transactionsFor.
*/

namespace PHARE::amr
{
// The copies a schedule makes into one patch from the patches of this process: SAMRAI's
// transactions, run in the order the schedule would have run them.
class PatchCopies : public SAMRAI::tbox::Transaction
{
    using Transaction = SAMRAI::tbox::Transaction;

public:
    // place: of the patch among those of its level on this process
    PatchCopies(int const rank, std::size_t const place)
        : rank_{rank}
        , place_{place}
    {
    }

    // The order of tbox::Schedule (addTransaction, appendTransaction, performLocalCopies):
    // transactions SAMRAI may fuse first, then the others; in each kind those within one patch
    // first, the last created first, then the others as they were created.
    void add(std::shared_ptr<Transaction> copy, bool const withinPatch)
    {
        auto& kind = std::dynamic_pointer_cast<SAMRAI::tbox::TransactionFuseable>(copy) ? fuseable_
                                                                                        : others_;
        if (withinPatch)
            kind.insert(kind.begin(), std::move(copy));
        else
            kind.push_back(std::move(copy));
    }

    // this transaction is the only one that writes its patch in the schedule
    void copyLocalData() override
    {
        auto const copyAll = [this]() {
            for (auto const& copy : fuseable_)
                copy->copyLocalData();
            for (auto const& copy : others_)
                copy->copyLocalData();
        };

        auto& jobs = ExchangeJobs::instance();
        if (jobs.recording())
            jobs.add(nullptr, copyAll, place_);
        else
            copyAll();
    }

    int getSourceProcessor() override { return rank_; }
    int getDestinationProcessor() override { return rank_; }

    // nothing goes through a message
    bool canEstimateIncomingMessageSize() override { return true; }
    std::size_t computeIncomingMessageSize() override { return 0; }
    std::size_t computeOutgoingMessageSize() override { return 0; }
    void packStream(SAMRAI::tbox::MessageStream&) override {}
    void unpackStream(SAMRAI::tbox::MessageStream&) override {}

    void printClassData(std::ostream& stream) const override
    {
        stream << "Copies into one patch\n";
        for (auto const* const kind : {&fuseable_, &others_})
            for (auto const& copy : *kind)
                copy->printClassData(stream);
    }

private:
    int rank_;
    std::size_t place_;
    std::vector<std::shared_ptr<Transaction>> fuseable_, others_;
};


// What a schedule is given in place of a transaction that a PatchCopies holds
class NoCopy : public SAMRAI::tbox::Transaction
{
public:
    explicit NoCopy(int const rank)
        : rank_{rank}
    {
    }

    void copyLocalData() override {}

    int getSourceProcessor() override { return rank_; }
    int getDestinationProcessor() override { return rank_; }

    bool canEstimateIncomingMessageSize() override { return true; }
    std::size_t computeIncomingMessageSize() override { return 0; }
    std::size_t computeOutgoingMessageSize() override { return 0; }
    void packStream(SAMRAI::tbox::MessageStream&) override {}
    void unpackStream(SAMRAI::tbox::MessageStream&) override {}
    void printClassData(std::ostream&) const override {}

private:
    int rank_;
};


// Makes the transactions of the factory it wraps (SAMRAI's own if none), and gathers those
// between two patches of this process in one PatchCopies per destination patch.
// One factory per schedule: it finds the PatchCopies of a patch again from one call to the next.
class ConcurrentRefineTransactionFactory : public SAMRAI::xfer::RefineTransactionFactory
{
    using Factory = SAMRAI::xfer::RefineTransactionFactory;
    using Level   = SAMRAI::hier::PatchLevel;
    using Box     = SAMRAI::hier::Box;

public:
    explicit ConcurrentRefineTransactionFactory(std::shared_ptr<Factory> wrapped = nullptr)
        : wrapped_{wrapped ? std::move(wrapped)
                           : std::make_shared<SAMRAI::xfer::StandardRefineTransactionFactory>()}
    {
    }

    std::shared_ptr<SAMRAI::tbox::Transaction>
    allocate(std::shared_ptr<Level> const& dst_level, std::shared_ptr<Level> const& src_level,
             std::shared_ptr<SAMRAI::hier::BoxOverlap> const& overlap, Box const& dst_box,
             Box const& src_box, SAMRAI::xfer::RefineClasses::Data const** refine_data,
             int const item_id, Box const& box, bool const use_time_interpolation) const override
    {
        auto transaction = wrapped_->allocate(dst_level, src_level, overlap, dst_box, src_box,
                                              refine_data, item_id, box, use_time_interpolation);

        auto const rank = dst_level->getBoxLevel()->getMPI().getRank();
        if (dst_box.getOwnerRank() != rank or src_box.getOwnerRank() != rank)
            return transaction;

        // A refine schedule (refine_data is its own) files its transactions in two tbox
        // schedules that it runs one after the other, depending on this flag of the quantity:
        // see RefineSchedule::constructScheduleTransactions, from which the next lines are
        bool const finePriority = refine_data[item_id]->d_fine_bdry_reps_var;
        bool const withinPatch
            = dst_level == src_level and dst_box.getGlobalId() == src_box.getGlobalId();

        auto& copies     = copies_[Key{refine_data, finePriority, dst_level.get(),
                                   dst_box.getLocalId().getValue()}];
        bool const first = !copies;
        // the patches of a level are walked in the order of their local ids
        if (first)
            copies = std::make_shared<PatchCopies>(
                rank, static_cast<std::size_t>(dst_box.getLocalId().getValue()));
        copies->add(std::move(transaction), withinPatch);

        if (first)
            return copies;
        if (!noCopy_)
            noCopy_ = std::make_shared<NoCopy>(rank);
        return noCopy_;
    }

    void setTransactionTime(double const fill_time) override
    {
        wrapped_->setTransactionTime(fill_time);
    }

    void
    preprocessScratchSpace(std::shared_ptr<Level> const& level, double const fill_time,
                           SAMRAI::hier::ComponentSelector const& preprocess_vector) const override
    {
        wrapped_->preprocessScratchSpace(level, fill_time, preprocess_vector);
    }

private:
    // quantities of a schedule, their tbox schedule, level and patch written
    using Key = std::tuple<void const*, bool, Level const*, int>;

    std::shared_ptr<Factory> wrapped_;
    mutable std::map<Key, std::shared_ptr<PatchCopies>> copies_;
    mutable std::shared_ptr<NoCopy> noCopy_;
};


// the factory to create a refine schedule with: the one it had with one thread (none: SAMRAI's
// own transactions), as without threads; otherwise one that gathers its copies patch by patch
inline std::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>
transactionsFor(std::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> factory = nullptr)
{
    if (core::nbrThreads() == 1)
        return factory;
    return std::make_shared<ConcurrentRefineTransactionFactory>(std::move(factory));
}

} // namespace PHARE::amr

#endif // PHARE_AMR_MESSENGERS_CONCURRENT_TRANSACTIONS_HPP
