#ifndef PHARE_LOHNER_PARAMETERS_HPP
#define PHARE_LOHNER_PARAMETERS_HPP

#include "initializer/data_provider.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{
/**
 * @brief The per-quantity parameters of both forms of Löhner's criterion.
 */
struct LohnerParameters
{
    double reltol;       ///< weight of the noise filter, relative to the local values
    double refMagnitude; ///< added to the denominator, in the units of the quantity

    /**
     * @brief Read the parameters of one tagged quantity.
     *
     * Both keys are required: their defaults live in pharein (pyphare/pyphare/pharein/tagging.py),
     * which always writes them.
     *
     * @param quantity config dict of the quantity, at path "simulation/AMR/refinement/tagging/Q<i>"
     * @param name The name of the quantity, for error messages
     * @return the parameters
     * @throws std::runtime_error if reltol or ref_magnitude is absent, negative or not finite
     */
    static LohnerParameters read(initializer::PHAREDict const& quantity, std::string const& name)
    {
        auto const value = [&](std::string const& key) {
            if (!quantity.contains(key))
                throw std::runtime_error("tagging quantity '" + name + "': " + key
                                         + " is required");
            auto const v = quantity[key].template to<double>();
            if (!(std::isfinite(v) and v >= 0.))
                throw std::runtime_error("tagging quantity '" + name + "': " + key
                                         + " must be finite and >= 0, got " + std::to_string(v));
            return v;
        };
        return {value("reltol"), value("ref_magnitude")};
    }
};

} // namespace PHARE::amr

#endif
