#ifndef PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_H
#define PHARE_CORE_NUMERIC_PUSHER_PUSHER_FACTORY_H

#include <cstddef>
#include <memory>
#include <string>

#include "boris.h"
#include "pusher.h"

namespace PHARE
{
namespace core
{
    class PusherFactory
    {
    public:
        template<std::size_t dim, typename ParticleIterator, typename Electromag,
                 typename Interpolator, typename BoundaryCondition, typename GridLayout>
        static auto makePusher(std::string pusherName)
        {
            if (pusherName == "modified_boris")
            {
                return std::make_unique<BorisPusher<dim, ParticleIterator, Electromag, Interpolator,
                                                    BoundaryCondition, GridLayout>>();
            }

            throw std::runtime_error("Error : Invalid Pusher name");
        }
    };
} // namespace core
} // namespace PHARE

#endif // PUSHER_FACTORY_H
