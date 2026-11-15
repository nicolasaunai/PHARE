#ifndef PHARE_IONS_HPP
#define PHARE_IONS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <stdexcept>


#include "core/def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/parallel.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"

namespace PHARE
{
namespace core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions
    {
    public:
        using value_type                = IonPopulation;
        using field_type                = IonPopulation::field_type;
        using vecfield_type             = IonPopulation::vecfield_type;
        using Float                     = field_type::type;
        using tensorfield_type          = IonPopulation::tensorfield_type;
        using particle_array_type       = IonPopulation::particle_array_type;
        using gridlayout_type           = GridLayout;
        static constexpr auto dimension = GridLayout::dimension;


        Ions(Ions const&) = default;
        Ions(Ions&&)      = default;


        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : massDensity_{massDensityName(), HybridQuantity::Scalar::rho}
            , chargeDensity_{chargeDensityName(), HybridQuantity::Scalar::rho}
            , bulkVelocity_{"bulkVel", HybridQuantity::Vector::V}
            , populations_{generate_from(
                  [&dict](auto ipop) { //
                      return IonPopulation{dict["pop" + std::to_string(ipop)]};
                  },
                  dict["nbrPopulations"].template to<std::size_t>())}
            , momentumTensor_{"momentumTensor", HybridQuantity::Tensor::M}
            , chargeFlux_{"chargeFlux", HybridQuantity::Vector::V}
        {
        }


        NO_DISCARD auto nbrPopulations() const { return populations_.size(); }
        NO_DISCARD auto size() const { return nbrPopulations(); }

        NO_DISCARD field_type const& massDensity() const { return massDensity_; }
        NO_DISCARD field_type& massDensity() { return massDensity_; }

        NO_DISCARD field_type const& chargeDensity() const { return chargeDensity_; }
        NO_DISCARD field_type& chargeDensity() { return chargeDensity_; }

        NO_DISCARD vecfield_type const& velocity() const { return bulkVelocity_; }
        NO_DISCARD vecfield_type& velocity() { return bulkVelocity_; }

        NO_DISCARD vecfield_type const& chargeFlux() const { return chargeFlux_; }
        NO_DISCARD vecfield_type& chargeFlux() { return chargeFlux_; }

        NO_DISCARD std::string static chargeDensityName() { return "chargeDensity"; }
        NO_DISCARD std::string static massDensityName() { return "massDensity"; }

        tensorfield_type const& momentumTensor() const { return momentumTensor_; }
        tensorfield_type& momentumTensor() { return momentumTensor_; }

        // the charge density of a population is its charge times its particle density,
        // which is complete (border sums, level ghost particles) when this is called
        // out = op(a) and out = op(a, b) node by node: in pieces run by the threads of the
        // process when the fields are large, see core::parallel_for_pieces
        static void transform_(auto const& a, auto& out, auto const& op)
        {
            parallel_for_pieces(a.size(), nodesPerPiece,
                                [&](std::size_t const first, std::size_t const last) {
                                    std::transform(std::begin(a) + first, std::begin(a) + last,
                                                   std::begin(out) + first, op);
                                });
        }

        static void combine_(auto const& a, auto const& b, auto& out, auto const& op)
        {
            parallel_for_pieces(
                a.size(), nodesPerPiece, [&](std::size_t const first, std::size_t const last) {
                    std::transform(std::begin(a) + first, std::begin(a) + last,
                                   std::begin(b) + first, std::begin(out) + first, op);
                });
        }

        void computeChargeDensity()
        {
            chargeDensity_.zero();

            for (auto& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto const& particleDensity = pop.particleDensity();
                auto& popDensity            = pop.chargeDensity();
                transform_(particleDensity, popDensity,
                           [q = pop.charge()](Float const& n) { return q * n; });
                combine_(chargeDensity_, popDensity, chargeDensity_, std::plus<Float>{});
            }
        }

        void computeMassDensity()
        {
            massDensity_.zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.particleDensity();
                combine_(
                    massDensity_, popDensity, massDensity_,
                    [&pop](auto const& n, auto const& pop_n) { return n + pop_n * pop.mass(); });
            }
        }


        void computeBulkVelocity()
        {
            computeMassDensity();

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                // account for mass only if populations have different masses
                std::function<Float(Float, Float)> plusMass
                    = [&pop](Float const& v, Float const& f) { return v + f * pop.mass(); };

                auto const& flux    = pop.flux();
                auto&& [fx, fy, fz] = flux();

                combine_(vx, fx, vx, plusMass);
                combine_(vy, fy, vy, plusMass);
                combine_(vz, fz, vz, plusMass);
            }


            combine_(vx, massDensity_, vx, std::divides<Float>{});
            combine_(vy, massDensity_, vy, std::divides<Float>{});
            combine_(vz, massDensity_, vz, std::divides<Float>{});
        }


        // sum_s q_s F_s, the charge flux entering the electron velocity
        void computeChargeFlux()
        {
            chargeFlux_.zero();

            for (auto const& pop : populations_)
            {
                auto const plusCharge
                    = [q = pop.charge()](Float const& c, Float const& f) { return c + q * f; };

                for (std::size_t i = 0; i < chargeFlux_.size(); ++i)
                {
                    auto& c       = chargeFlux_[i];
                    auto const& f = pop.flux()[i];
                    combine_(c, f, c, plusCharge);
                }
            }
        }


        void computeFullMomentumTensor()
        {
            momentumTensor_.zero();
            auto& mom = momentumTensor_;

            for (auto& pop : populations_)
            {
                auto& p_mom = pop.momentumTensor();
                for (auto p_mij = p_mom.begin(), mij = mom.begin(); p_mij != p_mom.end();
                     ++p_mij, ++mij)
                {
                    std::transform(std::begin(*mij), std::end(*mij), std::begin(*p_mij),
                                   std::begin(*mij), std::plus<typename field_type::type>{});
                }
            }
        }


        // drops the domain particles of all populations located in `boxes`
        template<typename Boxes>
        void eraseDomainParticlesIn(Boxes const& boxes)
        {
            if (boxes.empty())
                return;

            for (auto& pop : populations_)
            {
                auto& particles = pop.domainParticles();
                auto range      = makeIndexRange(particles);
                auto const kept = particles.partition(
                    range, [&](auto const& cell) { return !isIn(Point{cell}, boxes); });
                particles.erase(makeRange(particles, kept.iend(), particles.size()));
            }
        }


        NO_DISCARD auto begin() { return std::begin(populations_); }
        NO_DISCARD auto end() { return std::end(populations_); }
        NO_DISCARD auto begin() const { return std::begin(populations_); }
        NO_DISCARD auto end() const { return std::end(populations_); }

        NO_DISCARD auto& population(std::size_t const i)
        {
            if (i >= populations_.size())
                throw std::out_of_range("Ions population index out of range");
            return populations_[i];
        }

        NO_DISCARD auto const& population(std::size_t const i) const
        {
            if (i >= populations_.size())
                throw std::out_of_range("Ions population index out of range");
            return populations_[i];
        }

        NO_DISCARD bool isUsable() const
        {
            return core::isUsable(chargeDensity_, bulkVelocity_, momentumTensor_, massDensity_,
                                  chargeFlux_)
                   and core::all(populations_, [](auto const& pop) { return pop.isUsable(); });
        }



        NO_DISCARD bool isSettable() const
        {
            return core::isSettable(chargeDensity_, bulkVelocity_, momentumTensor_, massDensity_,
                                    chargeFlux_)
                   and core::all(populations_, [](auto const& pop) { return pop.isSettable(); });
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::vector<IonPopulation>& getRunTimeResourcesViewList()
        {
            return populations_;
        }

        NO_DISCARD auto getCompileTimeResourcesViewList()
        {
            return std::forward_as_tuple(bulkVelocity_, momentumTensor_, chargeDensity_,
                                         massDensity_, chargeFlux_);
        }



        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------


        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << "Ions\n";
            ss << "------------------------------------\n";
            ss << "number of populations  : " << nbrPopulations() << "\n";
            for (auto& pop : populations_)
                ss << core::to_str(pop);
            return ss.str();
        }


        auto& operator[](std::size_t const i) const { return population(i); }
        auto& operator[](std::size_t const i) { return population(i); }

    private:
        field_type massDensity_;
        field_type chargeDensity_;
        vecfield_type bulkVelocity_;
        std::vector<IonPopulation> populations_;
        tensorfield_type momentumTensor_;
        vecfield_type chargeFlux_;
    };

} // namespace core
} // namespace PHARE



#endif
