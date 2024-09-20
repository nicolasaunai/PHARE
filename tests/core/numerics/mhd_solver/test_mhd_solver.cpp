#include <gtest/gtest.h>

#include <core/data/grid/grid.hpp>
#include <core/data/grid/gridlayoutimplyee.hpp>
#include <core/mhd/mhd_quantities.hpp>

#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "phare_core.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"

constexpr std::uint32_t cells = 65;
constexpr std::size_t   dim = 1, interp = 1;

using YeeLayout_t       = PHARE::core::GridLayoutImplYeeMHD<dim, interp>;
using GridLayoutMHD     = PHARE::core::GridLayout<YeeLayout_t>;
using GridLayout_t      = TestGridLayout<GridLayoutMHD>;
using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;

TEST(UsableMHDStateTest, ConstructionTest) {
    GridLayout_t                     layout{cells};
    PHARE::core::UsableMHDState<dim> state(layout);

    ASSERT_NO_THROW();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
