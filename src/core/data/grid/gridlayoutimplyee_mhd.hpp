#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_HPP_MHD
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_HPP_MHD

#include <array>
#include <cstddef>
#include <vector>

#include "core/def.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/types.hpp"
#include "gridlayoutdefs.hpp"

namespace PHARE {
namespace core {
/**
 * @brief GridLayoutNdArrayImplYee class is a concrete GridLayoutNdArrayImpl used a Yee
 * type grid layout is needed.
 *
 * It provides methods related to grid layout operations:
 * - physical domain start/end indexes
 * - indexes of the first and last ghost nodes
 * - allocation sizes for Field attributes of other classes
 * - partial derivative operator (Faraday)
 * - physical coordinate given a field and a primal point (ix, iy, iz)
 * - cell centered coordinate given a primal point (ix, iy, iz)
 */
template <std::size_t dim, std::size_t interpOrder>
class GridLayoutImplYeeMHD {
    // ------------------------------------------------------------------------
    //                              PRIVATE
    // ------------------------------------------------------------------------
   public:
    static constexpr std::size_t      dimension    = dim;
    static constexpr std::size_t      interp_order = interpOrder;
    static constexpr std::string_view type         = "yee";
    using quantity_type                            = MHDQuantity;
    /*
    void constexpr initLinearCombinations_();

    LinearCombination momentsToEx_;
    LinearCombination momentsToEy_;
    LinearCombination momentsToEz_;
    LinearCombination BxToEy_;
    LinearCombination BxToEz_;
    LinearCombination ByToEx_;
    LinearCombination ByToEz_;
    LinearCombination BzToEx_;
    LinearCombination BzToEy_;
    LinearCombination ExToMoment_;
    LinearCombination EyToMoment_;
    LinearCombination EzToMoment_;
    */

    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize
     * the table MHDQuantityCentering_. This is THE important array in the GridLayout module.
     * This table knows which quantity is primal/dual along each direction. It is **this** array
     * that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee
     * interface are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_() {
        const gridDataT data{};

        const std::array<QtyCentering, NBR_COMPO> Rho = {{data.dual, data.dual, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Vx = {{data.dual, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> Vy = {{data.dual, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> Vz = {{data.dual, data.dual, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Bx = {{data.primal, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> By = {{data.dual, data.primal, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> Bz = {{data.dual, data.dual, data.primal}};

        const std::array<QtyCentering, NBR_COMPO> P = {{data.dual, data.dual, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Mx = {{data.dual, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> My = {{data.dual, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> Mz = {{data.dual, data.dual, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Etot = {{data.dual, data.dual, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Ex = {{data.dual, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Ey = {{data.primal, data.dual, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Ez = {{data.primal, data.primal, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Jx = {{data.dual, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Jy = {{data.primal, data.dual, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Jz = {{data.primal, data.primal, data.dual}};

        const std::array<std::array<QtyCentering, NBR_COMPO>,
                         static_cast<std::size_t>(MHDQuantity::Scalar::count)>
            _QtyCentering{Rho, Bx, By, Bz, Vx, Vy, Vz, P, Mx, My, Mz, Etot, Ex, Ey, Ez, Jx, Jy, Jz};

        return _QtyCentering;
    }

    //! says for each MHDQuantity::Quantity whether it is primal or dual, in each direction
    constexpr const static std::array<std::array<QtyCentering, NBR_COMPO>,
                                      static_cast<std::size_t>(MHDQuantity::Scalar::count)>
        _QtyCentering_{initLayoutCentering_()};

    static const std::size_t dim_{dim};

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
   public:
    NO_DISCARD constexpr static std::array<QtyCentering, dim> centering(
        MHDQuantity::Scalar MHDQuantity) {
        constexpr gridDataT_mhd gridData_{};
        if constexpr (dim == 1) {
            switch (MHDQuantity) {
                case MHDQuantity::Scalar::Bx:
                    return {{_QtyCentering_[gridData_.iBx][gridData_.idirX]}};
                case MHDQuantity::Scalar::By:
                    return {{_QtyCentering_[gridData_.iBy][gridData_.idirX]}};
                case MHDQuantity::Scalar::Bz:
                    return {{_QtyCentering_[gridData_.iBz][gridData_.idirX]}};
                case MHDQuantity::Scalar::Ex:
                    return {{_QtyCentering_[gridData_.iEx][gridData_.idirX]}};
                case MHDQuantity::Scalar::Ey:
                    return {{_QtyCentering_[gridData_.iEy][gridData_.idirX]}};
                case MHDQuantity::Scalar::Ez:
                    return {{_QtyCentering_[gridData_.iEz][gridData_.idirX]}};
                case MHDQuantity::Scalar::Jx:
                    return {{_QtyCentering_[gridData_.iJx][gridData_.idirX]}};
                case MHDQuantity::Scalar::Jy:
                    return {{_QtyCentering_[gridData_.iJy][gridData_.idirX]}};
                case MHDQuantity::Scalar::Jz:
                    return {{_QtyCentering_[gridData_.iJz][gridData_.idirX]}};
                case MHDQuantity::Scalar::rho:
                    return {{_QtyCentering_[gridData_.irho][gridData_.idirX]}};
                case MHDQuantity::Scalar::Vx:
                    return {{_QtyCentering_[gridData_.iVx][gridData_.idirX]}};
                case MHDQuantity::Scalar::Vy:
                    return {{_QtyCentering_[gridData_.iVy][gridData_.idirX]}};
                case MHDQuantity::Scalar::Vz:
                    return {{_QtyCentering_[gridData_.iVz][gridData_.idirX]}};
                case MHDQuantity::Scalar::P:
                    return {{_QtyCentering_[gridData_.iP][gridData_.idirX]}};
                case MHDQuantity::Scalar::Mx:
                    return {{_QtyCentering_[gridData_.iMx][gridData_.idirX]}};
                case MHDQuantity::Scalar::My:
                    return {{_QtyCentering_[gridData_.iMy][gridData_.idirX]}};
                case MHDQuantity::Scalar::Mz:
                    return {{_QtyCentering_[gridData_.iMz][gridData_.idirX]}};
                default:
                    throw std::runtime_error("Wrong MHDQuantity");
            }
        }

        else if constexpr (dim == 2) {
            switch (MHDQuantity) {
                case MHDQuantity::Scalar::Bx:
                    return {{_QtyCentering_[gridData_.iBx][gridData_.idirX],
                             _QtyCentering_[gridData_.iBx][gridData_.idirY]}};
                case MHDQuantity::Scalar::By:
                    return {{_QtyCentering_[gridData_.iBy][gridData_.idirX],
                             _QtyCentering_[gridData_.iBy][gridData_.idirY]}};
                case MHDQuantity::Scalar::Bz:
                    return {{_QtyCentering_[gridData_.iBz][gridData_.idirX],
                             _QtyCentering_[gridData_.iBz][gridData_.idirY]}};
                case MHDQuantity::Scalar::Ex:
                    return {{_QtyCentering_[gridData_.iEx][gridData_.idirX],
                             _QtyCentering_[gridData_.iEx][gridData_.idirY]}};
                case MHDQuantity::Scalar::Ey:
                    return {{_QtyCentering_[gridData_.iEy][gridData_.idirX],
                             _QtyCentering_[gridData_.iEy][gridData_.idirY]}};
                case MHDQuantity::Scalar::Ez:
                    return {{_QtyCentering_[gridData_.iEz][gridData_.idirX],
                             _QtyCentering_[gridData_.iEz][gridData_.idirY]}};
                case MHDQuantity::Scalar::Jx:
                    return {{_QtyCentering_[gridData_.iJx][gridData_.idirX],
                             _QtyCentering_[gridData_.iJx][gridData_.idirY]}};
                case MHDQuantity::Scalar::Jy:
                    return {{_QtyCentering_[gridData_.iJy][gridData_.idirX],
                             _QtyCentering_[gridData_.iJy][gridData_.idirY]}};
                case MHDQuantity::Scalar::Jz:
                    return {{_QtyCentering_[gridData_.iJz][gridData_.idirX],
                             _QtyCentering_[gridData_.iJz][gridData_.idirY]}};
                case MHDQuantity::Scalar::rho:
                    return {{_QtyCentering_[gridData_.irho][gridData_.idirX],
                             _QtyCentering_[gridData_.irho][gridData_.idirY]}};
                case MHDQuantity::Scalar::Vx:
                    return {{_QtyCentering_[gridData_.iVx][gridData_.idirX],
                             _QtyCentering_[gridData_.iVx][gridData_.idirY]}};
                case MHDQuantity::Scalar::Vy:
                    return {{_QtyCentering_[gridData_.iVy][gridData_.idirX],
                             _QtyCentering_[gridData_.iVy][gridData_.idirY]}};
                case MHDQuantity::Scalar::Vz:
                    return {{_QtyCentering_[gridData_.iVz][gridData_.idirX],
                             _QtyCentering_[gridData_.iVz][gridData_.idirY]}};
                case MHDQuantity::Scalar::P:
                    return {{_QtyCentering_[gridData_.iP][gridData_.idirX],
                             _QtyCentering_[gridData_.iP][gridData_.idirY]}};
                case MHDQuantity::Scalar::Mx:
                    return {{_QtyCentering_[gridData_.iMx][gridData_.idirX],
                             _QtyCentering_[gridData_.iMx][gridData_.idirY]}};
                case MHDQuantity::Scalar::My:
                    return {{_QtyCentering_[gridData_.iMy][gridData_.idirX],
                             _QtyCentering_[gridData_.iMy][gridData_.idirY]}};
                case MHDQuantity::Scalar::Mz:
                    return {{_QtyCentering_[gridData_.iMz][gridData_.idirX],
                             _QtyCentering_[gridData_.iMz][gridData_.idirY]}};
                default:
                    throw std::runtime_error("Wrong MHDQuantity");
            }
        }

        else if constexpr (dim == 3) {
            switch (MHDQuantity) {
                case MHDQuantity::Scalar::Bx:
                    return {{_QtyCentering_[gridData_.iBx][gridData_.idirX],
                             _QtyCentering_[gridData_.iBx][gridData_.idirY],
                             _QtyCentering_[gridData_.iBx][gridData_.idirZ]}};
                case MHDQuantity::Scalar::By:
                    return {{_QtyCentering_[gridData_.iBy][gridData_.idirX],
                             _QtyCentering_[gridData_.iBy][gridData_.idirY],
                             _QtyCentering_[gridData_.iBy][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Bz:
                    return {{_QtyCentering_[gridData_.iBz][gridData_.idirX],
                             _QtyCentering_[gridData_.iBz][gridData_.idirY],
                             _QtyCentering_[gridData_.iBz][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Ex:
                    return {{_QtyCentering_[gridData_.iEx][gridData_.idirX],
                             _QtyCentering_[gridData_.iEx][gridData_.idirY],
                             _QtyCentering_[gridData_.iEx][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Ey:
                    return {{_QtyCentering_[gridData_.iEy][gridData_.idirX],
                             _QtyCentering_[gridData_.iEy][gridData_.idirY],
                             _QtyCentering_[gridData_.iEy][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Ez:
                    return {{_QtyCentering_[gridData_.iEz][gridData_.idirX],
                             _QtyCentering_[gridData_.iEz][gridData_.idirY],
                             _QtyCentering_[gridData_.iEz][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Jx:
                    return {{_QtyCentering_[gridData_.iJx][gridData_.idirX],
                             _QtyCentering_[gridData_.iJx][gridData_.idirY],
                             _QtyCentering_[gridData_.iJx][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Jy:
                    return {{_QtyCentering_[gridData_.iJy][gridData_.idirX],
                             _QtyCentering_[gridData_.iJy][gridData_.idirY],
                             _QtyCentering_[gridData_.iJy][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Jz:
                    return {{_QtyCentering_[gridData_.iJz][gridData_.idirX],
                             _QtyCentering_[gridData_.iJz][gridData_.idirY],
                             _QtyCentering_[gridData_.iJz][gridData_.idirZ]}};
                case MHDQuantity::Scalar::rho:
                    return {{_QtyCentering_[gridData_.irho][gridData_.idirX],
                             _QtyCentering_[gridData_.irho][gridData_.idirY],
                             _QtyCentering_[gridData_.irho][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Vx:
                    return {{_QtyCentering_[gridData_.iVx][gridData_.idirX],
                             _QtyCentering_[gridData_.iVx][gridData_.idirY],
                             _QtyCentering_[gridData_.iVx][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Vy:
                    return {{_QtyCentering_[gridData_.iVy][gridData_.idirX],
                             _QtyCentering_[gridData_.iVy][gridData_.idirY],
                             _QtyCentering_[gridData_.iVy][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Vz:
                    return {{_QtyCentering_[gridData_.iVz][gridData_.idirX],
                             _QtyCentering_[gridData_.iVz][gridData_.idirY],
                             _QtyCentering_[gridData_.iVz][gridData_.idirZ]}};
                case MHDQuantity::Scalar::P:
                    return {{_QtyCentering_[gridData_.iP][gridData_.idirX],
                             _QtyCentering_[gridData_.iP][gridData_.idirY],
                             _QtyCentering_[gridData_.iP][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Mx:
                    return {{_QtyCentering_[gridData_.iMx][gridData_.idirX],
                             _QtyCentering_[gridData_.iMx][gridData_.idirY],
                             _QtyCentering_[gridData_.iMx][gridData_.idirZ]}};
                case MHDQuantity::Scalar::My:
                    return {{_QtyCentering_[gridData_.iMy][gridData_.idirX],
                             _QtyCentering_[gridData_.iMy][gridData_.idirY],
                             _QtyCentering_[gridData_.iMy][gridData_.idirZ]}};
                case MHDQuantity::Scalar::Mz:
                    return {{_QtyCentering_[gridData_.iMz][gridData_.idirX],
                             _QtyCentering_[gridData_.iMz][gridData_.idirY],
                             _QtyCentering_[gridData_.iMz][gridData_.idirZ]}};
                default:
                    throw std::runtime_error("Wrong MHDQuantity");
            }
        }
    }

    NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3> centering(
        MHDQuantity::Vector MHDQuantity) {
        switch (MHDQuantity) {
            case MHDQuantity::Vector::B:
                return {{centering(MHDQuantity::Scalar::Bx), centering(MHDQuantity::Scalar::By),
                         centering(MHDQuantity::Scalar::Bz)}};

            case MHDQuantity::Vector::V:
                return {{centering(MHDQuantity::Scalar::Vx), centering(MHDQuantity::Scalar::Vy),
                         centering(MHDQuantity::Scalar::Vz)}};

            case MHDQuantity::Vector::J:
                return {{centering(MHDQuantity::Scalar::Jx), centering(MHDQuantity::Scalar::Jy),
                         centering(MHDQuantity::Scalar::Jz)}};

            case MHDQuantity::Vector::E:
                return {{centering(MHDQuantity::Scalar::Ex), centering(MHDQuantity::Scalar::Ey),
                         centering(MHDQuantity::Scalar::Ez)}};

            default:
                throw std::runtime_error("Wrong MHDQuantity");
        }
    }

    NO_DISCARD auto static constexpr dualToPrimal() { return -1; }

    NO_DISCARD auto static constexpr primalToDual() { return 1; }
};  // namespace core

}  // namespace core
}  // namespace PHARE

#endif  // PHARE_CORE_GRID_GRIDLAYOUTYEE_HPP
