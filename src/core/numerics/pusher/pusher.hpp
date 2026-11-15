#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP

#include <span>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <functional>

#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle.hpp"

namespace PHARE
{
namespace core
{
    template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
             typename BoundaryCondition, typename GridLayout, bool UseParticleIndex = true>
    class Pusher
    {
    protected:
        static auto constexpr dimension = GridLayout::dimension;

    public:
        using ParticleSelector = std::function<ParticleRange(ParticleRange&)>;

        // TODO : to really be independant on boris which has 2 push steps
        // we should have an arbitrary number of selectors, 1 per push step
        // charge and mass are those of the population the particles belong to
        // an empty firstSelector selects every particle
        virtual ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                                   Electromag const& emFields, double charge, double mass,
                                   Interpolator& interpolator, GridLayout const& layout,
                                   ParticleSelector firstSelector, ParticleSelector secondSelector)
            = 0;


        using Particle_t = typename ParticleRange::array_t::Particle_t;

        // what moveInPieces and pushCopiesInPieces do with the particles of a block once pushed
        using PushedBlock = std::function<void(std::size_t piece, std::span<Particle_t const>)>;

        /** Pushes the particles of range in their array, as move does when no particle is
         * selected, and hands them over while they are still in the cache of the core.
         *
         * The range is cut in nPieces consecutive pieces run by the threads of the process, see
         * core::parallel_for. A piece is pushed block by block of particles that fit in the
         * cache of a core, and pushed(piece, block) is called on each block right after its
         * push, by the thread that pushed it: the blocks of a piece one after the other in the
         * order of the particles, different pieces concurrently.
         */
        virtual void moveInPieces(ParticleRange& range, Electromag const& emFields, double charge,
                                  double mass, Interpolator& interpolator, GridLayout const& layout,
                                  std::size_t nPieces, PushedBlock const& pushed)
            = 0;

        /** Same as moveInPieces on copies of the particles that only exist for pushed: the
         * particles of range and their array are left as they are.
         */
        virtual void pushCopiesInPieces(ParticleRange const& range, Electromag const& emFields,
                                        double charge, double mass, Interpolator& interpolator,
                                        GridLayout const& layout, std::size_t nPieces,
                                        PushedBlock const& pushed)
            = 0;


        virtual void setMeshAndTimeStep(std::array<double, dim> ms, double ts) = 0;

        virtual ~Pusher() {}
    };

} // namespace core
} // namespace PHARE

#endif
