#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/utilities/parallel.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "initializer/data_provider.hpp"


#include <span>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <algorithm>


namespace PHARE::core
{
enum class UpdaterMode { domain_only = 1, all = 2 };

template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
    using This = IonUpdater<Ions, Electromag, GridLayout>;

public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField          = Ions::vecfield_type;
    using Field             = Ions::field_type;
    using ParticleArray     = Ions::particle_array_type;
    using Particle_t        = ParticleArray::Particle_t;
    using PartIterator      = ParticleArray::iterator;
    using ParticleRange     = IndexRange<ParticleArray>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::Pusher<dimension, ParticleRange, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleRange, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }

    template<typename Boxing_t>
    void updatePopulations(Ions& ions, Electromag const& em, Boxing_t const& boxing, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions);


    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(ParticleArray{Box{}});
    }


private:
    template<typename Boxing_t>
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Boxing_t>
    void updateAndDepositAll_(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Population, typename Boxing_t>
    void pushCopyAndDepositDomain_(Population& pop, Electromag const& em, Boxing_t const& boxing);


    // adds the particles of range to the particle density and flux of their population
    template<typename Population>
    void deposit_(ParticleRange range, Population& pop, GridLayout const& layout);

    // pushes the domain particles of the population and adds to its density and flux those that
    // end in a cell deposited(cell) accepts
    template<typename Population>
    void moveAndDeposit_(Population& pop, Electromag const& em, GridLayout const& layout,
                         bool keepPushed, auto const& deposited);

    // density and flux a deposit piece writes on, see depositTargets_
    struct DepositTarget
    {
        Field density;
        VecField flux;
        double* copy; // what both are views on, null for the fields of the population
    };

    template<typename Population>
    std::vector<DepositTarget> depositTargets_(Population& pop, std::size_t nPieces);

    template<typename Population>
    void addDepositCopies_(Population& pop, std::size_t nPieces);

    // pieces a deposit of nParticles is cut in, 1 with one thread. It must not depend on the
    // number of threads, see depositTargets_
    static std::size_t depositPieces_(std::size_t const nParticles, std::size_t const nNodes)
    {
        if (nbrThreads() == 1)
            return 1;

        // a piece zeroes and adds one copy of the fields: it has to hold many more particles
        // than the fields have nodes
        auto const particlesPerPiece = std::max<std::size_t>(16384, 4 * nNodes);
        return std::clamp<std::size_t>(nParticles / particlesPerPiece, 1, 64);
    }

    // dealloced on regridding/load balancing coarsest
    ParticleArray tmp_particles_{Box{}}; //{std::make_unique<ParticleArray>(Box{})};

    // density and flux copies of the deposit pieces
    std::vector<double> depositBuffer_;
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 Boxing_t const& boxing, double dt,
                                                                 UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);
    pusher_->setMeshAndTimeStep(boxing.layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
    {
        updateAndDepositDomain_(ions, em, boxing);
    }
    else
    {
        updateAndDepositAll_(ions, em, boxing);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
    ions.computeChargeFlux();
}

/**
 * What the pieces of a deposit write on: the density and flux of the population for the first
 * one, zeroed copies of them for the others, since two pieces may write the same nodes.
 * addDepositCopies_ then adds the copies to the fields, in the order of the pieces.
 *
 * A node sums its particles piece by piece instead of one by one: the result differs from the
 * one-piece deposit by rounding. The pieces depend on the number of particles and of nodes
 * only, see depositPieces_, so the result is the same whatever the number of threads above one.
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Population>
auto IonUpdater<Ions, Electromag, GridLayout>::depositTargets_(Population& pop,
                                                               std::size_t const nPieces)
    -> std::vector<DepositTarget>
{
    auto& density = pop.particleDensity();
    auto& flux    = pop.flux();

    auto constexpr nFields = std::size_t{4}; // density and flux
    auto const nNodes      = density.size();
    auto const copySize    = nFields * nNodes;
    depositBuffer_.resize((nPieces - 1) * copySize);

    std::vector<DepositTarget> targets;
    targets.reserve(nPieces);
    targets.push_back({density, flux, nullptr});

    for (std::size_t p = 1; p < nPieces; ++p)
    {
        auto* const copy = depositBuffer_.data() + (p - 1) * copySize;
        targets.push_back(
            {Field{density.name(), density.physicalQuantity(), copy, density.shape()}, flux, copy});
        for (std::size_t i = 0; i < 3; ++i)
            targets.back().flux[i] = Field{flux[i].name(), flux[i].physicalQuantity(),
                                           copy + (i + 1) * nNodes, flux[i].shape()};
    }
    return targets;
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename Population>
void IonUpdater<Ions, Electromag, GridLayout>::addDepositCopies_(Population& pop,
                                                                 std::size_t const nPieces)
{
    auto& density       = pop.particleDensity();
    auto& flux          = pop.flux();
    auto const nNodes   = density.size();
    auto const copySize = 4 * nNodes;

    // nodes are independent of each other: pieces of them can be summed concurrently
    auto const nNodePieces = nPieces > 1 ? nbrChunks(nNodes, 4096) : 0;
    parallel_for(nNodePieces, [&](std::size_t const np) {
        auto const nodes = chunk(nNodes, nNodePieces, np);
        for (std::size_t p = 1; p < nPieces; ++p)
        {
            auto const* const copy = depositBuffer_.data() + (p - 1) * copySize;
            auto const add         = [&](auto& field, double const* const from) {
                auto* const to = field.data();
                for (auto i = nodes.first; i < nodes.last; ++i)
                    to[i] += from[i];
            };
            add(density, copy);
            for (std::size_t i = 0; i < 3; ++i)
                add(flux[i], copy + (i + 1) * nNodes);
        }
    });
}


/**
 * A long range is cut in consecutive pieces deposited by the threads of the process, see
 * core::parallel_for and depositTargets_.
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Population>
void IonUpdater<Ions, Electromag, GridLayout>::deposit_(ParticleRange range, Population& pop,
                                                        GridLayout const& layout)
{
    auto const n       = range.size();
    auto const nPieces = depositPieces_(n, pop.particleDensity().size());

    if (nPieces == 1)
    {
        interpolator_(range, pop.particleDensity(), pop.flux(), layout);
        return;
    }

    auto targets = depositTargets_(pop, nPieces);

    parallel_for(nPieces, [&](std::size_t const p) {
        auto const piece = chunk(n, nPieces, p);
        auto pieceRange
            = makeRange(range.array(), range.ibegin() + piece.first, range.ibegin() + piece.last);
        Interpolator pieceInterpolator; // its work arrays are per piece

        auto& [density, flux, copy] = targets[p];
        if (copy)
            std::fill_n(copy, 4 * density.size(), 0.);
        pieceInterpolator(pieceRange, density, flux, layout);
    });

    addDepositCopies_(pop, nPieces);
}


/**
 * The push and the deposit of the domain particles in one pass over them: the pusher hands
 * each block of particles it has pushed to the deposit while the block is still in the cache
 * of the core, instead of the deposit reading all the particles again from memory.
 *
 * Without keepPushed the pushed particles are not stored at all: the predictor only needs
 * them for the deposit, and the particles of the population stay as they are.
 *
 * The pieces are those of deposit_, each writing on its own target: see depositTargets_.
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Population>
void IonUpdater<Ions, Electromag, GridLayout>::moveAndDeposit_(Population& pop,
                                                               Electromag const& em,
                                                               GridLayout const& layout,
                                                               bool const keepPushed,
                                                               auto const& deposited)
{
    auto range         = makeIndexRange(pop.domainParticles());
    auto const nPieces = depositPieces_(range.size(), pop.particleDensity().size());
    auto targets       = depositTargets_(pop, nPieces);

    // a piece zeroes its copy when it starts: its blocks come one after the other
    std::vector<char> started(nPieces, 0);

    auto const deposit = [&](std::size_t const p, std::span<Particle_t const> const block) {
        auto& [density, flux, copy] = targets[p];
        if (copy and !started[p])
        {
            std::fill_n(copy, 4 * density.size(), 0.);
            started[p] = 1;
        }

        Interpolator blockInterpolator; // its work arrays are per block
        blockInterpolator(
            std::views::filter(block,
                               [&](auto const& particle) { return deposited(particle.iCell); }),
            density, flux, layout);
    };

    if (keepPushed)
        pusher_->moveInPieces(range, em, pop.charge(), pop.mass(), interpolator_, layout, nPieces,
                              deposit);
    else
        pusher_->pushCopiesInPieces(range, em, pop.charge(), pop.mass(), interpolator_, layout,
                                    nPieces, deposit);

    // a piece without particles never started
    for (std::size_t p = 1; p < nPieces; ++p)
        if (!started[p])
            std::fill_n(targets[p].copy, 4 * targets[p].density.size(), 0.);

    addDepositCopies_(pop, nPieces);
}


// this is to detach how we partition particles from the updater directly
template<typename IonUpdater_t, typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = IonUpdater_t::Box;
    using Selector_t                     = IonUpdater_t::Pusher::ParticleSelector;

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    // domain cells owned by an overlapping same-level patch, usually empty
    std::vector<Box_t> const nonOwnedBoxes{};
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    static bool isOwned(auto const& cell, Box_t const& domainBox,
                        std::vector<Box_t> const& nonOwnedBoxes)
    {
        return isIn(cell, domainBox)
               and (nonOwnedBoxes.empty() or !isIn(Point{cell}, nonOwnedBoxes));
    }

    // a particle belongs to this patch if and only if its cell is in the domain and owned
    bool isOwned(auto const& cell) const { return isOwned(cell, domainBox, nonOwnedBoxes); }

    // what inNonLevelGhostBox selects, for one cell: the domain is part of these boxes
    bool isInNonLevelGhostBox(auto const& cell) const
    {
        return isIn(cell, domainBox) or isIn(Point{cell}, nonLevelGhostBox);
    }

    // no selection: tells the pusher that it keeps every particle between its two position
    // updates, see Pusher::move
    Selector_t const allParticles{};

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inOwnedDomain
        = [domainBox = domainBox, nonOwnedBoxes = nonOwnedBoxes](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isOwned(cell, domainBox, nonOwnedBoxes);
              });
          };

    Selector_t const inGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return isIn(cell, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [nonLevelGhostBox = nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(Point{cell}, nonLevelGhostBox);
              });
          };

    Selector_t const inGhostLayer
        = [ghostBox = ghostBox, domainBox = domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, domainBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return !isIn(cell, ghostBox); });
    };
};

/**
 * The predictor's push and deposit of the domain particles with one thread: a copy of the
 * particles is pushed, selected and deposited, as three passes over them.
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Population, typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::pushCopyAndDepositDomain_(Population& pop,
                                                                         Electromag const& em,
                                                                         Boxing_t const& boxing)
{
    auto const& layout = boxing.layout;

    // work on a local copy, made by the pusher as it reads the domain particles.
    // Its cell map starts empty instead of as a copy: the pusher registers in it the
    // particles that change cell, and those are all the selection after the push looks
    // for, since no domain particle starts in a cell that inNonLevelGhostBox rejects.
    tmp_particles_.resize_like(pop.domainParticles());

    auto const inRange = makeIndexRange(pop.domainParticles());
    auto outRange      = makeIndexRange(tmp_particles_);
    auto allowed       = outRange
        = pusher_->move(inRange, outRange, em, pop.charge(), pop.mass(), interpolator_, layout,
                        boxing.allParticles, boxing.inNonLevelGhostBox);

    deposit_(allowed, pop, layout);
}


/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                       Electromag const& em,
                                                                       Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositDomain_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        // both ways below rely on it
        PHARE_DEBUG_DO({
            for (auto const& particle : pop.domainParticles())
                if (!isIn(Point{particle.iCell}, boxing.nonLevelGhostBox))
                    throw std::runtime_error(
                        "IonUpdater: a domain particle is outside the non level ghost box");
        })

        // push all domain particles, accumulate those inNonLevelGhostBox
        if (nbrThreads() > 1)
            // the pushed domain particles are only needed by the deposit
            moveAndDeposit_(pop, em, layout, /*keepPushed=*/false,
                            [&](auto const& cell) { return boxing.isInNonLevelGhostBox(cell); });
        else
            pushCopyAndDepositDomain_(pop, em, boxing);

        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inOwnedDomain

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(tmp_particles_);

            auto enteredInDomain
                = pusher_->move(outRange, outRange, em, pop.charge(), pop.mass(), interpolator_,
                                layout, boxing.inGhostBox, boxing.inOwnedDomain);

            deposit_(enteredInDomain, pop, layout);
        };

        // After this function is done domain particles overlaping ghost layers of neighbor patches
        // are sent to these neighbor's patchghost particle array.
        // After being pushed, some patch ghost particles may enter the domain. These need to be
        // copied into the domain array so they are transfered to the neighbor patch
        // ghost array and contribute to moments there too.
        // On the contrary level ghost particles entering the domain here do not need to be copied
        // since they contribute to nodes that are not shared with neighbor patches an since
        // level border nodes will receive contributions from levelghost old and new particles

        if (pop.levelGhostParticles().size())
            pushAndAccumulateGhosts(pop.levelGhostParticles());
    }
}


/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
                                                                    Electromag const& em,
                                                                    Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositAll_");

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        // with threads the domain particles are deposited as they are pushed, those that end
        // where the deposits below would find them: in the domain or among the new patch ghosts
        bool const depositedInPush = nbrThreads() > 1;
        auto& patchGhost           = pop.patchGhostParticles();
        auto const oldPatchGhosts  = patchGhost.size();

        if (depositedInPush)
            moveAndDeposit_(pop, em, layout, /*keepPushed=*/true,
                            [&](auto const& cell) { return boxing.isInNonLevelGhostBox(cell); });

        auto inDomain
            = depositedInPush
                  ? boxing.inOwnedDomain(domainPartRange)
                  : pusher_->move(domainPartRange, domainPartRange, em, pop.charge(), pop.mass(),
                                  interpolator_, layout, boxing.allParticles, boxing.inOwnedDomain);

        auto now_ghosts = makeRange(domainParticles, inDomain.iend(), domainParticles.size());
        auto const not_level_ghosts = boxing.inNonLevelGhostBox(now_ghosts);

        // copy out new patch ghosts
        patchGhost.reserve(patchGhost.size() + not_level_ghosts.size());
        std::copy(not_level_ghosts.begin(), not_level_ghosts.end(), std::back_inserter(patchGhost));

        PHARE_DEBUG_DO({
            auto const outsideGhostBox = boxing.outsideGhostBox(now_ghosts);
            for (auto const& particle : outsideGhostBox)
            {
                PHARE_LOG_LINE_SS(particle);
                auto const nearbyBox = grow(Box(particle.iCell, particle.iCell), 3);
                for (auto const& xyz : em.E)
                    if (auto const overlap = nearbyBox * layout.AMRGhostBoxFor(xyz))
                        layout.amr_lcl_idx(*overlap, [&](auto const& bix, auto const& lix) {
                            PHARE_LOG_LINE_SS(xyz.name() << " at:" << bix << ":" << xyz(lix));
                        });
            }
            if (outsideGhostBox.size())
                throw core::DictionaryException{}("ID", "Updater::outsideGhostBox");
        })

        domainParticles.erase(now_ghosts); // drop all ghosts
        auto const pushedInDomain = domainParticles.size();

        if (pop.levelGhostParticles().size())
        {
            auto particleRange = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange
                = pusher_->move(particleRange, particleRange, em, pop.charge(), pop.mass(),
                                interpolator_, layout, boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(domainParticles,
                                           [&](auto const& cell) { return boxing.isOwned(cell); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        if (depositedInPush)
        {
            // what the push has not seen: level ghosts that entered, older patch ghosts
            deposit_(makeRange(domainParticles, pushedInDomain, domainParticles.size()), pop,
                     layout);
            deposit_(makeRange(patchGhost, 0, oldPatchGhosts), pop, layout);
        }
        else
        {
            deposit_(makeIndexRange(domainParticles), pop, layout);
            deposit_(makeIndexRange(patchGhost), pop, layout);
        }
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
