#include <gtest/gtest.h>

#include <core/data/grid/grid.hpp>
#include <core/data/grid/gridlayoutimplyee.hpp>
#include <core/mhd/mhd_quantities.hpp>

#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "phare_core.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"

constexpr std::uint32_t cells = 65;
constexpr std::size_t   dim = 1, interp = 1;

using YeeLayout_t       = PHARE::core::GridLayoutImplYee<dim, interp>;
using GridLayout_t      = TestGridLayout<YeeLayout_t>;
using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;

/*template <std::size_t dimension_>*/
/*struct DummyLayout {*/
/*    static constexpr std::size_t        dimension = dimension_;*/
/*    std::array<unsigned int, dimension> nbrCells_;*/
/*    auto                                nbrCells() const { return nbrCells_; }*/
/*    auto AMRBox() const { return PHARE::core::emptyBox<int, dimension>(); }*/
/*};*/

TEST(UsableMHDStateTest, ConstructionTest) {
    GridLayout_t                     layout{cells};
    PHARE::core::UsableMHDState<dim> state(layout);

    ASSERT_NO_THROW();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
