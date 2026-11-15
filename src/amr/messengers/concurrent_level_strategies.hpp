#ifndef PHARE_AMR_MESSENGERS_CONCURRENT_LEVEL_STRATEGIES_HPP
#define PHARE_AMR_MESSENGERS_CONCURRENT_LEVEL_STRATEGIES_HPP

#include "core/utilities/parallel.hpp"

#include "amr/data/exchange_jobs.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/xfer/RefinePatchStrategy.h>
#include <SAMRAI/xfer/CoarsenPatchStrategy.h>

/*
  SAMRAI interpolates a fine level from the next coarser one patch after patch, on one thread
  (RefineSchedule::refineScratchData): for each patch it calls the refine operator of each
  quantity, and the patch strategy of the schedule before and after them. It coarsens a level the
  same way (CoarsenSchedule::coarsenSourceData). The patch strategy is also told when the loop
  over the patches starts and when it is over.

  The strategies below use these two calls to open and close a phase of ExchangeJobs. In between,
  PHARE's operators record their work instead of doing it, and the phase runs it on the threads
  of the process when the loop is over: what is done for a patch on one thread, in the order of
  SAMRAI's calls. Operators only write the patch they are called for, and read the other level.

  With one thread nothing is recorded: schedules are given the strategy they had, see the
  functions at the end.
*/

namespace PHARE::amr
{
// Hands every call to the strategy it wraps, if any. What that strategy does before and after
// the operators of a patch is recorded with them.
class ConcurrentRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
    using Strategy  = SAMRAI::xfer::RefinePatchStrategy;
    using Patch     = SAMRAI::hier::Patch;
    using Level     = SAMRAI::hier::PatchLevel;
    using IntVector = SAMRAI::hier::IntVector;
    using Boxes     = SAMRAI::hier::BoxContainer;
    using Connector = SAMRAI::hier::Connector;

public:
    // wrapped is not owned, and may be nullptr
    explicit ConcurrentRefinePatchStrategy(Strategy* const wrapped = nullptr)
        : wrapped_{wrapped}
    {
    }

    void setPhysicalBoundaryConditions(Patch& patch, double const fill_time,
                                       IntVector const& ghost_width_to_fill) override
    {
        if (wrapped_)
            wrapped_->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }

    IntVector getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return wrapped_ ? wrapped_->getRefineOpStencilWidth(dim) : IntVector::getZero(dim);
    }

    void preprocessRefine(Patch& fine, Patch const& coarse, SAMRAI::hier::Box const& fine_box,
                          IntVector const& ratio) override
    {
        if (wrapped_)
            wrapped_->preprocessRefine(fine, coarse, fine_box, ratio);
    }

    void postprocessRefine(Patch& fine, Patch const& coarse, SAMRAI::hier::Box const& fine_box,
                           IntVector const& ratio) override
    {
        if (wrapped_)
            wrapped_->postprocessRefine(fine, coarse, fine_box, ratio);
    }

    // RefineSchedule calls the box versions, once per patch: the boxes and the ratio are locals
    // of its loop
    void preprocessRefineBoxes(Patch& fine, Patch const& coarse, Boxes const& fine_boxes,
                               IntVector const& ratio) override
    {
        if (wrapped_)
            copyOrDefer(&fine, [strategy = wrapped_, &fine, &coarse,
                                boxes = &forJobs(fine_boxes), kept = &forJobs(ratio)]() {
                strategy->preprocessRefineBoxes(fine, coarse, *boxes, *kept);
            });
    }

    void postprocessRefineBoxes(Patch& fine, Patch const& coarse, Boxes const& fine_boxes,
                                IntVector const& ratio) override
    {
        if (wrapped_)
            copyOrDefer(&fine, [strategy = wrapped_, &fine, &coarse,
                                boxes = &forJobs(fine_boxes), kept = &forJobs(ratio)]() {
                strategy->postprocessRefineBoxes(fine, coarse, *boxes, *kept);
            });
    }

    void preprocessRefineLevel(
        Level& fine_level, Level const& coarse_level, Connector const& coarse_to_fine,
        Connector const& coarse_to_unfilled,
        std::vector<std::vector<std::shared_ptr<SAMRAI::hier::BoxOverlap>>> const& overlaps,
        SAMRAI::xfer::RefineClasses::Data const** refine_items) override
    {
        if (wrapped_)
            wrapped_->preprocessRefineLevel(fine_level, coarse_level, coarse_to_fine,
                                            coarse_to_unfilled, overlaps, refine_items);
        ExchangeJobs::instance().begin();
    }

    void postprocessRefineLevel(Level& fine_level, Level const& coarse_level,
                                Connector const& coarse_to_fine,
                                Connector const& coarse_to_unfilled) override
    {
        ExchangeJobs::instance().end();
        if (wrapped_)
            wrapped_->postprocessRefineLevel(fine_level, coarse_level, coarse_to_fine,
                                             coarse_to_unfilled);
        setNeedRefineSynchronize(false);
    }

private:
    Strategy* wrapped_;
};


// PHARE gives no strategy to its coarsen schedules: this one only opens and closes the phase
class ConcurrentCoarsenPatchStrategy : public SAMRAI::xfer::CoarsenPatchStrategy
{
    using Patch     = SAMRAI::hier::Patch;
    using Level     = SAMRAI::hier::PatchLevel;
    using IntVector = SAMRAI::hier::IntVector;

public:
    IntVector getCoarsenOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return IntVector::getZero(dim);
    }

    void preprocessCoarsen(Patch&, Patch const&, SAMRAI::hier::Box const&,
                           IntVector const&) override
    {
    }

    void postprocessCoarsen(Patch&, Patch const&, SAMRAI::hier::Box const&,
                            IntVector const&) override
    {
    }

    void preprocessCoarsenLevel(Level&, Level const&) override
    {
        ExchangeJobs::instance().begin();
        setNeedCoarsenSynchronize(false);
    }

    void postprocessCoarsenLevel(Level&, Level const&) override
    {
        ExchangeJobs::instance().end();
        setNeedCoarsenSynchronize(false);
    }
};


// the strategy to create a refine schedule with: strategy itself with one thread, as without
// threads, otherwise onThreads, which wraps it
inline SAMRAI::xfer::RefinePatchStrategy*
refineStrategyFor(SAMRAI::xfer::RefinePatchStrategy* const strategy,
                  ConcurrentRefinePatchStrategy& onThreads)
{
    return core::nbrThreads() > 1 ? &onThreads : strategy;
}

// the strategy to create a coarsen schedule with: none with one thread, as without threads
inline SAMRAI::xfer::CoarsenPatchStrategy*
coarsenStrategyFor(ConcurrentCoarsenPatchStrategy& onThreads)
{
    return core::nbrThreads() > 1 ? &onThreads : nullptr;
}

// schedule with its copies, packs and unpacks run on the process threads, see ExchangeJobs
template<typename Schedule>
std::shared_ptr<Schedule> onThreads(std::shared_ptr<Schedule> schedule)
{
    schedule->setScheduleOpsStrategy(&ExchangeJobs::instance());
    return schedule;
}

} // namespace PHARE::amr

#endif // PHARE_AMR_MESSENGERS_CONCURRENT_LEVEL_STRATEGIES_HPP
