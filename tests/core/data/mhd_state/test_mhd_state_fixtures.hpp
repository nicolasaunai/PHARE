#ifndef PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_MHDSTATE_MHDSTATE_FIXTURES_HPP

#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"
#include "phare_core.hpp"
#include "tests/core/data/field/test_field_fixtures_mhd.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace PHARE::core {

template <typename State, typename GridLayout>
auto default_mhd_init(State &state, GridLayout const &layout) {
    auto constexpr static dim = GridLayout::dimension;

    auto setter = [&](auto &v) {
        auto box = layout.ghostBoxFor(v);
        for (std::size_t i = 0; i < box.upper[0]; i += 2) {
            if constexpr (dim == 1) {
                v(std::array{i}) += .05;
            } else {
                for (std::size_t j = 0; j < box.upper[1]; j += 2) {
                    if constexpr (dim == 2) {
                        v(i, j) += .05;
                    } else {
                        if constexpr (dim == 3) {
                            for (std::size_t k = 0; k < box.upper[2]; k += 2) {
                                v(i, j, k) += .05;
                            }
                        }
                    }
                }
            }
        }
    };
    for (auto &xyz : state.V) setter(xyz);
    for (auto &xyz : state.B) setter(xyz);
    setter(state.rho);
    setter(state.P);
}

template <std::size_t dim>
class UsableMHDState : public MHDState<VecFieldMHD<dim>> {
    using Array_t = NdArrayVector<dim, double, /*c_ordering*/ true>;
    using Grid_t  = Grid<Array_t, HybridQuantity::Scalar>;
    using Super   = MHDState<VecFieldMHD<dim>>;

    void _set() {
        auto &&[_rho, _V, _B, _P, _M, _Etot, _J] = Super::getCompileTimeResourcesViewList();
        _rho.setBuffer(&rho);
        V.set_on(_V);
        B.set_on(_B);
        _P.setBuffer(&P);
    }

   public:
    template <typename GridLayout>
    UsableMHDState(GridLayout const &layout)
        : Super{"MHD"},
          rho{"MHD_rho", layout, MHDQuantity::Scalar::rho},
          V{"MHD_V", layout, MHDQuantity::Vector::V},
          B{"MHD_B", layout, MHDQuantity::Vector::B},
          P{"MHD_P", layout, MHDQuantity::Scalar::P} {
        _set();
        default_mhd_init(super(), layout);
    }

    UsableMHDState(UsableMHDState &&that)
        : Super{std::forward<Super>(that)},
          rho{std::move(that.rho)},
          V{std::move(that.V)},
          B{std::move(that.B)},
          P{std::move(that.P)} {
        _set();
    }

    Super       &super() { return *this; }
    Super const &super() const { return *this; }
    auto        &operator*() { return super(); }
    auto        &operator*() const { return super(); }

    Grid_t                 rho, P;
    UsableVecFieldMHD<dim> V, B;
};

}  // namespace PHARE::core

#endif
