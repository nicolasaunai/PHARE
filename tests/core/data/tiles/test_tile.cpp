#include <ctype.h>
#include <string>

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "core/data/tiles/tile_set.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>
#include <unordered_set>


using namespace PHARE::core;


template<typename TileSet>
class TileTestBase
{
public:
    static auto constexpr dimension = TileSet::dimension;

    TileTestBase(Box<int, dimension> box_, std::array<std::size_t, dimension> const& tile_size)
        : box{box_}
        , tileSet{box, tile_size}
    {
    }

    Box<int, dimension> box;
    TileSet tileSet;
};

template<typename TileSet>
class TileTestBoxShapeNotMultipleTileSize : public TileTestBase<TileSet>, public ::testing::Test
{
public:
    TileTestBoxShapeNotMultipleTileSize()
        : TileTestBase<TileSet>{
              Box<int, TileSet::dimension>{ConstArray<int, TileSet::dimension>(0),
                                           ConstArray<int, TileSet::dimension>(54)},
              ConstArray<std::size_t, TileSet::dimension>(4)}
    {
    }
};

template<typename TileSet>
class TileTestBoxShapeMultipleTileSize : public TileTestBase<TileSet>, public ::testing::Test
{
public:
    TileTestBoxShapeMultipleTileSize()
        : TileTestBase<TileSet>{
              Box<int, TileSet::dimension>{ConstArray<int, TileSet::dimension>(0),
                                           ConstArray<int, TileSet::dimension>(47)},
              ConstArray<std::size_t, TileSet::dimension>(4)}
    {
    }
};


template<typename TileSet>
class TileTest : public ::testing::Test
{
};


template<std::size_t dim>
class TileMock : public Box<int, dim>
{
};

using DimTiles = testing::Types<TileSet<TileMock<1>>, TileSet<TileMock<2>>, TileSet<TileMock<3>>>;

TYPED_TEST_SUITE(TileTestBoxShapeNotMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTestBoxShapeMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTest, DimTiles);


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, expectedNbrOfTilesPerDimToCoverTheBox)
{
    auto const& shape = this->tileSet.shape();
    for (auto i = 0u; i < this->dimension; ++i)
        EXPECT_EQ(shape[i], 14);
}

TYPED_TEST(TileTestBoxShapeMultipleTileSize, cluserSetSizeIsCorrect)
{
    EXPECT_EQ(this->tileSet.size(), std::pow(12, this->dimension));
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, totalTileSetSurfaceIsEqualToBoxSurface)
{
    auto surface = 0.;
    for (auto i = 0u; i < this->tileSet.size(); ++i)
    {
        auto current_surface = 1.;
        for (auto d = 0u; d < this->dimension; ++d)
        {
            auto l = (this->tileSet[i].upper[d] - this->tileSet[i].lower[d] + 1);
            current_surface *= l;
        }
        surface += current_surface;
    }
    EXPECT_EQ(surface, this->box.size());
}



TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, tileHasNoOverlapWithOthers)
{
    auto constexpr dim = TypeParam::dimension;
    for (auto const& tile : this->tileSet)
    {
        for (auto const& other : this->tileSet)
        {
            if (&tile != &other)
            {
                auto const box1 = Box<int, dim>{tile.lower, tile.upper};
                auto const box2 = Box<int, dim>{other.lower, other.upper};
                auto overlap    = box1 * box2;
                EXPECT_FALSE(overlap.has_value());
            }
        }
    }
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, retrieveTilesFromBoxOverlap)
{
    auto constexpr dim = TypeParam::dimension;
    Box<int, dim> selection_box{ConstArray<int, dim>(11), ConstArray<int, dim>(34)};

    auto expected_nbr = std::pow(7, this->dimension);
    auto overlapeds   = this->tileSet.overlaped_with(selection_box);
    EXPECT_EQ(overlapeds.size(), expected_nbr);

    auto completes   = 0.;
    auto incompletes = 0.;
    for (auto const& overlaped : overlapeds)
    {
        auto const& [is_complete, tile] = overlaped;
        if (is_complete)
            ++completes;
        else
            ++incompletes;
    }
    EXPECT_EQ(completes, std::pow(5, dim));
    EXPECT_EQ(incompletes, std::pow(7, dim) - std::pow(5, dim));
}


TYPED_TEST(TileTest, cannotCreateTileWithTileSizeBiggerThanBox)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(5)};
    auto const tile_size = PHARE::core::ConstArray<std::size_t, dim>(7); // larger than box shape
    EXPECT_THROW(std::make_unique<TypeParam>(box, tile_size), std::runtime_error);
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, canRetrieveTileFromCell)
{
    auto constexpr dim = TypeParam::dimension;
    auto tile          = [&]() {
        if constexpr (dim == 1)
            return this->tileSet.at(13);
        else if constexpr (dim == 2)
            return this->tileSet.at(13, 13);
        else if constexpr (dim == 3)
            return this->tileSet.at(13, 13, 13);
    }();
    auto const expected_box = Box<int, dim>{ConstArray<int, dim>(12), ConstArray<int, dim>(15)};
    EXPECT_TRUE(*tile == expected_box);
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, useTileSetView)
{
    auto view  = this->tileSet.make_view();
    auto shape = view.shape();
    for (auto const& tile : view)
    {
        EXPECT_LE(tile.size(), std::pow(4, this->dimension));
    }
}

struct PatchData
{
    using Box2D     = Box<int, 2>;
    using TileSet2D = TileSet<Box2D>;
    Box<int, 2> box;
    TileSet2D tileSet;

    PatchData(Box2D box_,
              std::array<std::size_t, 2> const& tile_size = ConstArray<std::size_t, 2>(4))
        : box{box_}
        , tileSet{box, tile_size}
    {
    }
};


TEST(TileSetViewSpan, fromManyPatches)
{
    using Box2D = Box<int, 2>;
    std::array<Box2D, 3> boxes{Box2D{{0, 10}, {24, 58}}, Box2D{{33, 57}, {65, 83}},
                               Box2D{{102, 99}, {134, 128}}};

    std::vector<PatchData> patchdatas;
    std::vector<TileSetView<Box2D>> views;

    for (auto const& box : boxes)
    {
        patchdatas.emplace_back(box);
    }


    // now making views
    for (auto /*const? make_view fails if so...*/& patch : patchdatas)
    {
        views.push_back(patch.tileSet.make_view());
    }
}




TYPED_TEST(TileTest, InnerTileHaveCorrectNbrOfNeighbors)
{
    // build the tileSet for a given box
    auto constexpr tile_size = 4u;
    auto constexpr dim       = TypeParam::dimension;
    using BoxND              = Box<int, dim>;

    auto const lower = ConstArray<int, dim>(0);
    auto const upper = ConstArray<int, dim>(51);
    BoxND box{lower, upper};
    TileSet<BoxND> tileSet{box, ConstArray<std::size_t, dim>(tile_size)};

    // get the tiles that overlap the inner box
    auto inner_box{box};
    inner_box.grow(-tile_size);
    auto inner_tiles = tileSet.overlaped_with(inner_box);

    auto constexpr expected_neighbor_nbr = [&]() {
        if constexpr (dim == 1)
            return 2;
        else if constexpr (dim == 2)
            return 8;
        else if constexpr (dim == 3)
            return 26;
    }();

    for (auto const& [complete_overlap, tile] : inner_tiles)
    {
        if (complete_overlap)
        {
            std::unordered_set<BoxND*> neighbors;
            BoxND ghost_box{*tile};
            ghost_box.grow(1);
            for (auto const& cell : ghost_box)
            {
                auto tile_ptr = [&]() {
                    if constexpr (dim == 1)
                        return tileSet.at(cell[0]);
                    else if constexpr (dim == 2)
                        return tileSet.at(cell[0], cell[1]);
                    else
                        return tileSet.at(cell[0], cell[1], cell[2]);
                }();
                if (!isIn(cell, *tile))
                    neighbors.insert(tile_ptr);
            }
            EXPECT_EQ(neighbors.size(), expected_neighbor_nbr);
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
