

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/utilities/box/box.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <numeric>
#include <algorithm>

using namespace PHARE::core;

std::size_t constexpr static dim = 1;
using GridLayoutT                = PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t;

// field loops cut in pieces for the threads of a process visit the same points, once
TEST(GridLayout, evalOnBoxInPiecesVisitsThePointsOfEvalOnBoxOnce)
{
    using GridLayout2D = PHARE_Types<PHARE::SimOpts{2, 1}>::Hybrid::GridLayout_t;
    GridLayout2D layout{{0.1, 0.1}, {300u, 200u}, {{0., 0.}}, Box{Point{0, 0}, Point{299, 199}}};

    struct
    {
        auto physicalQuantity() const { return HybridQuantity::Scalar::rho; }
    } rho;

    std::size_t const ny = 400;
    std::vector<int> onOneThread(ny * ny, 0), inPieces(ny * ny, 0), bigger(ny * ny, 0);

    layout.evalOnBox(rho, [&](auto const i, auto const j) { ++onOneThread[i * ny + j]; });

    setNbrThreads(4);
    layout.evalOnBoxInPieces(rho, [&](auto const i, auto const j) { ++inPieces[i * ny + j]; });
    layout.evalOnBiggerBoxInPieces(rho, Point<std::uint32_t, 2>{1u, 1u},
                                   [&](auto const i, auto const j) { ++bigger[i * ny + j]; });
    setNbrThreads(1);

    EXPECT_EQ(onOneThread, inPieces);
    EXPECT_EQ(std::accumulate(inPieces.begin(), inPieces.end(), 0), 301 * 201);
    EXPECT_EQ(std::accumulate(bigger.begin(), bigger.end(), 0), 303 * 203);
    EXPECT_EQ(*std::max_element(bigger.begin(), bigger.end()), 1);
}


TEST(GridLayout, isGivenAnAMRIndexBoxAtConstruction)
{
    GridLayoutT layout({0.1}, {50u}, {{0.}}, Box{Point{0}, Point{49}});
}


auto badLayout()
{
    auto nbrCells = 50;
    return GridLayoutT{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells}}};
}


auto goodLayout()
{
    auto nbrCells = 50;
    return GridLayoutT{
        {0.1}, {static_cast<std::uint32_t>(nbrCells)}, {{0.}}, Box{Point{0}, Point{nbrCells - 1}}};
}



TEST(GridLayout, AMRBoxHasNbrCellsCells)
{
    EXPECT_ANY_THROW({ badLayout(); });
    EXPECT_NO_THROW({ goodLayout(); });
}


TEST(GridLayout, canTransformALocalIndexIntoAnAMRIndex)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts = layout.options.field_ghost_width;

    EXPECT_EQ(Point{50}, layout.localToAMR(Point{nGhosts}));
    EXPECT_EQ(Point{85}, layout.localToAMR(Point{nGhosts + 35}));
}


TEST(GridLayout, canTransformALocalBoxIntoAnAMRBox)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts         = layout.options.field_ghost_width;
    auto localBox       = Box{Point{nGhosts + 5}, Point{nGhosts + 15}};
    auto expectedAMRBox = Box{Point{55}, Point{65}};

    EXPECT_EQ(expectedAMRBox, layout.localToAMR(localBox));
}



TEST(GridLayout, canTransformAnAMRIndexIntoALocalIndex)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    int nGhosts = layout.options.field_ghost_width;
    EXPECT_EQ(Point{nGhosts}, layout.AMRToLocal(Point{50}));
    EXPECT_EQ(Point{nGhosts + 35}, layout.AMRToLocal(Point{85}));
}




TEST(GridLayout, canTransformAnAMRBoxIntoALocalBox)
{
    GridLayoutT layout{{0.1}, {50u}, {{0.}}, Box{Point{50}, Point{99}}};

    std::uint32_t nGhosts = layout.options.field_ghost_width;
    auto AMRBox           = Box<int, dim>{Point<int, dim>{55}, Point<int, dim>{65}};
    auto expectedLocalBox = Box<std::uint32_t, dim>{Point<std::uint32_t, dim>{nGhosts + 5},
                                                    Point<std::uint32_t, dim>{nGhosts + 15}};

    EXPECT_EQ(expectedLocalBox, layout.AMRToLocal(AMRBox));
}
