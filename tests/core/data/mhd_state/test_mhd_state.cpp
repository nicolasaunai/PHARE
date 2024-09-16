#include <type_traits>

#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/vecfield/vecfield.hpp"

#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"

using namespace PHARE::core;
using namespace PHARE::initializer;
using namespace PHARE::initializer::test_fn::func_1d;

using Field_t          = Field<1, MHDQuantity::Scalar>;
using VecField1D       = VecField<Field_t, MHDQuantity>;


PHAREDict getDict()
{
    using initfunc = InitFunction<1>;
    PHAREDict dict;
    
    dict["density"]["initializer"] = static_cast<initfunc>(density);
    
    dict["velocity"]["initializer"]["x_component"] = static_cast<initfunc>(vx);
    dict["velocity"]["initialiser"]["y_component"] = static_cast<initfunc>(vy);
    dict["velocity"]["initialiser"]["z_component"] = static_cast<initfunc>(vz);
    
    dict["magnetic"]["initialiser"]["x_component"] = static_cast<initfunc>(bx);
    dict["magnetic"]["initialiser"]["y_component"] = static_cast<initfunc>(by);
    dict["magnetic"]["initialiser"]["z_component"] = static_cast<initfunc>(bz);
    
    dict["pressure"]["initialiser"] = static_cast<initfunc>(pressure);
    
    return dict;
}

struct AnMHDState : public ::testing::Test
{
    MHDState<VecField1D> state{getDict()};
    virtual ~AnMHDState();
};

AnMHDState::~AnMHDState() {}

TEST_F(AnMHDState, noUsableFieldsUponConstruction)
{
    EXPECT_FALSE(state.isUsable());
}

TEST_F(AnMHDState, fieldsAreSettable)
{
    EXPECT_TRUE(state.isSettable());
}

TEST_F(AnMHDState, hasTupleResourceList)
{
    auto resources = state.getCompileTimeResourcesViewList();
    [[maybe_unused]] auto& rho = std::get<0>(resources);
    [[maybe_unused]] auto& v = std::get<1>(resources);
    [[maybe_unused]] auto& b = std::get<2>(resources);
    [[maybe_unused]] auto& p = std::get<3>(resources);
    [[maybe_unused]] auto& m = std::get<4>(resources);
    [[maybe_unused]] auto& etot = std::get<5>(resources);
    [[maybe_unused]] auto& j = std::get<6>(resources);

    EXPECT_FALSE(rho.isUsable());
    EXPECT_FALSE(v.isUsable());
    EXPECT_FALSE(b.isUsable());
    EXPECT_FALSE(p.isUsable());
    EXPECT_FALSE(m.isUsable());
    EXPECT_FALSE(etot.isUsable());
    EXPECT_FALSE(j.isUsable());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
