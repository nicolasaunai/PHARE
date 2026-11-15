#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/def.hpp"
#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/utilities/parallel.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <exception>

namespace PHARE::core
{


template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
    : public Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
    struct MoveTwoCellException : std::exception
    {
        MoveTwoCellException(double const d, double const v)
            : delta{d}
            , vel{v}
        {
        }

        double delta, vel;
    };


public:
    using Super
        = Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>;

private:
    using ParticleSelector = Super::ParticleSelector;
    using Particle_t       = Super::Particle_t;
    using ParticleArray_t  = ParticleRange::array_t;
    using PushedBlock      = Super::PushedBlock;

public:
    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    /** see Pusher::move() documentation*/
#if 0
    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double charge, double mass,
                       Interpolator& interpolator,
                       ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                       GridLayout const& layout) override
    {
            // push the particles of half a step
            // rangeIn : t=n, rangeOut : t=n+1/Z
            // get a pointer on the first particle of rangeOut that leaves the patch
            auto firstLeaving
                = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

            // apply boundary condition on the particles in [firstLeaving, rangeOut.end[
            // that actually leave through a physical boundary condition
            // get a pointer on the new end of rangeOut. Particles passed newEnd
            // are those that have left the patch through a non-physical boundary
            // they should be discarded now
            auto newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

            rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

            // get electromagnetic fields interpolated on the particles of rangeOut
            // stop at newEnd.
            interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

            // get the particle velocity from t=n to t=n+1
            accelerate_(rangeOut, rangeOut, mass);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

            // apply BC on the leaving particles that leave through physical BC
            // and get pointer on new End, discarding particles leaving elsewhere
            newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

            rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

            return rangeOut.end();
    }
#endif


    /** see Pusher::move() documentation
     *
     * With an empty firstSelector no particle is set aside between the two position updates:
     * the range is then pushed in pieces, see moveInPieces_, with the same result.
     */
    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double charge, double mass,
                       Interpolator& interpolator, GridLayout const& layout,
                       ParticleSelector firstSelector, ParticleSelector secondSelector) override
    {
        PHARE_LOG_SCOPE(3, "Boris::move_no_bc");

        double const dto2m = 0.5 * dt_ / mass;
        double const coef1 = charge * dto2m;

        if (!firstSelector)
        {
            moveInPieces_(rangeIn, rangeOut, emFields, coef1, interpolator, layout);
            return secondSelector(rangeOut);
        }

        auto& outParticles    = rangeOut.array();
        auto const changeCell = [&](auto const& newCell, std::size_t const idx) {
            outParticles.change_icell(newCell, idx);
        };

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
             ++inIdx, ++outIdx)
            prePushStep_(rangeIn.array()[inIdx], outParticles[outIdx], outIdx, changeCell);

        rangeOut = firstSelector(rangeOut);

        for (auto idx = rangeOut.ibegin(); idx < rangeOut.iend(); ++idx)
            accelerateAndPostPushStep_(outParticles[idx], idx, emFields, coef1, interpolator,
                                       layout, changeCell);

        return secondSelector(rangeOut);
    }



    /** see Pusher::moveInPieces() documentation */
    void moveInPieces(ParticleRange& range, Electromag const& emFields, double charge, double mass,
                      Interpolator& interpolator, GridLayout const& layout, std::size_t nPieces,
                      PushedBlock const& pushed) override
    {
        blocksOfPieces_(range, &range.array(), emFields, charge * 0.5 * dt_ / mass, interpolator,
                        layout, nPieces, pushed);
    }

    /** see Pusher::pushCopiesInPieces() documentation */
    void pushCopiesInPieces(ParticleRange const& range, Electromag const& emFields, double charge,
                            double mass, Interpolator& interpolator, GridLayout const& layout,
                            std::size_t nPieces, PushedBlock const& pushed) override
    {
        blocksOfPieces_(range, nullptr, emFields, charge * 0.5 * dt_ / mass, interpolator, layout,
                        nPieces, pushed);
    }



    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename Particle>
    PHARE_ALWAYS_INLINE auto advancePosition_(Particle const& partIn, Particle& partOut)
    {
        std::array<int, dim> newCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double const delta
                = partIn.delta[iDim] + static_cast<double>(halfDtOverDl_[iDim] * partIn.v[iDim]);

            if (std::abs(delta) > 2)
                throw MoveTwoCellException{delta, partIn.v[iDim]};

            auto const iCell    = static_cast<int>(std::floor(delta));
            partOut.delta[iDim] = delta - iCell;
            newCell[iDim]       = iCell + partIn.iCell[iDim];
        }
        return newCell;
    }


    // particles of a piece: they fit in the cache of a core (256 kB in 2D), so that the piece
    // goes through memory once for the whole push
    static constexpr std::size_t particlesPerPiece_ = 4096;

    // particles of one piece that changed cell in one of the two position updates
    struct alignas(64) CellChanges
    {
        struct Change
        {
            std::size_t idx;
            std::array<int, dim> oldCell, newCell;
        };
        std::vector<Change> changes;
    };

    /** Pushes rangeIn into rangeOut as the two loops of move do, piece by piece: the first
     * position update of a piece of particles, then their velocity and second position update,
     * while they are still in the cache of the core.
     *
     * The pieces run on the threads of the process, see core::parallel_for. The cell map of
     * the output array is shared by the pieces, so a piece notes the particles that change
     * cell, and the map learns the changes afterwards in the order one thread running the two
     * loops over the whole range would have made them: all those of the first position update
     * in particle order, then all those of the second. The map and the particles end the same.
     *
     * When the two ranges are in different arrays, the output particles are copies of the input
     * ones made here: the output array only needs the size and the cell map of the input one.
     */
    void moveInPieces_(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double const coef1, Interpolator& interpolator,
                       GridLayout const& layout)
    {
        auto& inParticles   = rangeIn.array();
        auto& outParticles  = rangeOut.array();
        auto const firstIn  = rangeIn.ibegin();
        auto const firstOut = rangeOut.ibegin();
        auto const n        = rangeIn.size();
        auto const nPieces  = std::max<std::size_t>(n / particlesPerPiece_, 1);

        std::vector<CellChanges> prePush(nPieces), postPush(nPieces);

        parallel_for(nPieces, [&](std::size_t const p) {
            auto const piece       = chunk(n, nPieces, p);
            auto pieceInterpolator = interpolator; // its work arrays are per piece

            if (piece.first < piece.last)
                pushBlock_(&inParticles[firstIn + piece.first],
                           &outParticles[firstOut + piece.first], piece.last - piece.first,
                           firstOut + piece.first, emFields, coef1, pieceInterpolator, layout,
                           noteIn_(prePush[p], outParticles), noteIn_(postPush[p], outParticles));
        });

        replay_(prePush, postPush, outParticles);
    }


    /** moveInPieces and pushCopiesInPieces: the pushed particles are those of the range when
     * kept is their array, copies that are not kept when it is null.
     *
     * Particles pushed in their array change the cell map as in moveInPieces_.
     */
    void blocksOfPieces_(ParticleRange const& range, ParticleArray_t* const kept,
                         Electromag const& emFields, double const coef1, Interpolator& interpolator,
                         GridLayout const& layout, std::size_t const nPieces,
                         PushedBlock const& pushed)
    {
        auto const& particles = range.array();
        auto const first      = range.ibegin();
        auto const n          = range.size();

        std::vector<CellChanges> prePush(kept ? nPieces : 0), postPush(kept ? nPieces : 0);

        parallel_for(nPieces, [&](std::size_t const p) {
            auto const piece       = chunk(n, nPieces, p);
            auto pieceInterpolator = interpolator; // its work arrays are per piece

            for (auto b = piece.first; b < piece.last; b += particlesPerPiece_)
            {
                auto const count     = std::min(particlesPerPiece_, piece.last - b);
                auto const* const in = &particles[first + b];

                if (kept)
                {
                    auto* const block = &(*kept)[first + b];
                    pushBlock_(in, block, count, first + b, emFields, coef1, pieceInterpolator,
                               layout, noteIn_(prePush[p], *kept), noteIn_(postPush[p], *kept));
                    pushed(p, {block, count});
                }
                else
                {
                    // one block of copies per thread: a block is done with before the thread
                    // starts another
                    static thread_local std::vector<Particle_t> copies(particlesPerPiece_);

                    auto const setCell = [&](auto const& newCell, std::size_t const idx) {
                        copies[idx].iCell = newCell;
                    };
                    pushBlock_(in, copies.data(), count, 0, emFields, coef1, pieceInterpolator,
                               layout, setCell, setCell);
                    pushed(p, {copies.data(), count});
                }
            }
        });

        if (kept)
            replay_(prePush, postPush, *kept);
    }


    /** Pushes count particles from in to out (the same, or copies): their first position update,
     * then their velocity and second position update. out[i] is at index firstIdx + i for the
     * two changeCell(newCell, index), called for the particles that change cell in the first and
     * in the second position update.
     */
    PHARE_ALWAYS_INLINE void pushBlock_(Particle_t const* const in, Particle_t* const out,
                                        std::size_t const count, std::size_t const firstIdx,
                                        Electromag const& emFields, double const coef1,
                                        Interpolator& interpolator, GridLayout const& layout,
                                        auto const& changeCellPrePush,
                                        auto const& changeCellPostPush)
    {
        for (std::size_t i = 0; i < count; ++i)
        {
            if (in != out)
                out[i] = in[i];
            prePushStep_(in[i], out[i], firstIdx + i, changeCellPrePush);
        }

        for (std::size_t i = 0; i < count; ++i)
            accelerateAndPostPushStep_(out[i], firstIdx + i, emFields, coef1, interpolator, layout,
                                       changeCellPostPush);
    }

    // notes in changes the particles of the array that change cell, for replay_
    static auto noteIn_(CellChanges& changes, ParticleArray_t& particles)
    {
        return [&](auto const& newCell, std::size_t const idx) {
            changes.changes.push_back({idx, particles[idx].iCell, newCell});
            particles[idx].iCell = newCell;
        };
    }

    // the cell map of the array learns the changes of cell noted by the pieces, in the order one
    // thread pushing the whole range in two loops would have made them
    static void replay_(std::vector<CellChanges> const& prePush,
                        std::vector<CellChanges> const& postPush, ParticleArray_t& particles)
    {
        for (auto const* const pass : {&prePush, &postPush})
            for (auto const& piece : *pass)
                for (auto const& [idx, oldCell, newCell] : piece.changes)
                    particles.cell_changed(oldCell, newCell, idx);
    }


    /** advance partIn of half a time step and store it in partOut, at index outIdx of its array
     *  changeCell(newCell, outIdx) is called if it changes cell
     */
    template<typename Particle_t>
    PHARE_ALWAYS_INLINE void prePushStep_(Particle_t const& partIn, Particle_t& partOut,
                                          std::size_t const outIdx, auto const& changeCell)
    {
        // in the first push, this is the first time
        // we push to rangeOut, which contains crap
        // the push will only touch the particle position
        // but the next step being the acceleration of
        // rangeOut, we need to copy rangeIn weights
        // and velocity. This is done here although
        // not strictly speaking this function's business
        // to take advantage that we're already looping
        // over rangeIn particles.

        partOut.weight = partIn.weight;
        partOut.v      = partIn.v;

        try
        {
            auto newCell = advancePosition_(partIn, partOut);
            if (newCell != partIn.iCell)
                changeCell(newCell, outIdx);
        }
        catch (MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PrePush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            DictionaryException ex{"cause", ss.str()};
            throw ex;
        }
    }

    /** get the particle velocity from t=n to t=n+1 with the electromagnetic fields interpolated
     *  at its position, then advance it from t=n+1/2 to t=n+1 using v_{n+1} just calculated.
     *  changeCell(newCell, idx) is called if it changes cell
     */
    template<typename Particle_t>
    PHARE_ALWAYS_INLINE void
    accelerateAndPostPushStep_(Particle_t& particle, std::size_t const idx,
                               Electromag const& emFields, double const coef1,
                               Interpolator& interpolator, GridLayout const& layout,
                               auto const& changeCell)
    {
        auto const& local_em = interpolator(particle, emFields, layout);
        accelerate_(particle, local_em, coef1);

        try
        {
            auto newCell = advancePosition_(particle, particle);
            if (newCell != particle.iCell)
                changeCell(newCell, idx);
        }
        catch (MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PostPush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            DictionaryException ex{"cause", ss.str()};
            auto const& [E, B] = local_em;
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("E_" + std::to_string(i), std::to_string(E[i]));
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("B_" + std::to_string(i), std::to_string(B[i]));
            ex("level", std::to_string(layout.levelNumber()));
            throw ex;
        }
    }


    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename Particle_t, typename ParticleEB>
    PHARE_ALWAYS_INLINE void accelerate_(Particle_t& part, ParticleEB const& particleEB,
                                         double const coef1)
    {
        auto& [pE, pB]        = particleEB;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        // We now apply the 3 steps of the BORIS PUSHER

        // 1st half push of the electric field
        double velx1 = part.v[0] + coef1 * pEx;
        double vely1 = part.v[1] + coef1 * pEy;
        double velz1 = part.v[2] + coef1 * pEz;


        // preparing variables for magnetic rotation
        double const rx = coef1 * pBx;
        double const ry = coef1 * pBy;
        double const rz = coef1 * pBz;

        double const rx2  = rx * rx;
        double const ry2  = ry * ry;
        double const rz2  = rz * rz;
        double const rxry = rx * ry;
        double const rxrz = rx * rz;
        double const ryrz = ry * rz;

        double const invDet = 1. / (1. + rx2 + ry2 + rz2);

        // preparing rotation matrix due to the magnetic field
        // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
        double const mxx = 1. + rx2 - ry2 - rz2;
        double const mxy = 2. * (rxry + rz);
        double const mxz = 2. * (rxrz - ry);

        double const myx = 2. * (rxry - rz);
        double const myy = 1. + ry2 - rx2 - rz2;
        double const myz = 2. * (ryrz + rx);

        double const mzx = 2. * (rxrz + ry);
        double const mzy = 2. * (ryrz - rx);
        double const mzz = 1. + rz2 - rx2 - ry2;

        // magnetic rotation
        double const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
        double const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
        double const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;


        // 2nd half push of the electric field
        velx1 = velx2 + coef1 * pEx;
        vely1 = vely2 + coef1 * pEy;
        velz1 = velz2 + coef1 * pEz;

        // Update particle velocity
        part.v[0] = velx1;
        part.v[1] = vely1;
        part.v[2] = velz1;
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
