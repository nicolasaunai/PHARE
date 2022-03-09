#include <string>
#include <vector>

#include "core/utilities/bucketlist.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;



TEST(BucketList, canBeGivenAParticle)
{
    BucketList<100> bl;
    bl.add(0);
}



TEST(BucketList, registerMoreThanBucketSize)
{
    BucketList<2> bl;
    bl.add(0);
    bl.add(1);
    bl.add(2);
}



TEST(BucketList, sizeEqualsNumberOfRegisteredItems)
{
    BucketList<3> bl;
    bl.add(0);
    bl.add(1);
    bl.add(2);
    bl.add(3);

    EXPECT_EQ(4, bl.size());

    bl.add(4);
    bl.add(5);
    bl.add(6);
    bl.add(7);

    EXPECT_EQ(8, bl.size());
}


TEST(BucketList, capcityEqualsTotalNbrOfMemorySlots)
{
    BucketList<3> bl;
    EXPECT_EQ(3, bl.capacity());
    BucketList<100> bl100;
    EXPECT_EQ(100, bl100.capacity());

    bl.add(0);
    bl.add(1);
    bl.add(2);
    EXPECT_EQ(3, bl.capacity());
    EXPECT_EQ(bl.size(), bl.capacity());

    bl.add(4);
    EXPECT_EQ(6, bl.capacity());
}



TEST(BucketList, emptySetSizeZeroLeavingCapacityUnchanged)
{
    BucketList<3> bl;
    EXPECT_EQ(3, bl.capacity());
    EXPECT_EQ(0, bl.size());

    bl.add(0);
    bl.add(1);
    bl.add(2);
    bl.add(3);

    EXPECT_EQ(4, bl.size());
    EXPECT_FALSE(bl.is_empty());

    bl.empty();
    EXPECT_TRUE(bl.is_empty());

    EXPECT_EQ(0, bl.size());
    EXPECT_EQ(6, bl.capacity());
}




TEST(BucketList, beginReturnsIteratorOnFirstElement)
{
    BucketList<3> bl;
    std::array<int, 4> values = {3, 4, 5, 6};
    bl.add(0);
    bl.add(1);
    bl.add(2);
    bl.add(3);
    int first = values[*std::begin(bl)];
    EXPECT_EQ(3, first);
}


TEST(BucketList, endReturnsIteratorOnLastPlusOneElement)
{
    std::vector<int> values({1, 2, 3, 4, 5});
    BucketList<3> bl;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        bl.add(itemIndex);

    int actual;
    for (auto const& itemIndex : bl)
    {
        actual = values[itemIndex];
    }
    EXPECT_EQ(actual, 5);
}

TEST(BucketList, singleElementBucketListHasCorrectEnd)
{
    std::vector<int> values({18});
    BucketList<3> bl;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        bl.add(itemIndex);

    int actual;
    for (auto const& itemIndex : bl)
    {
        actual = values[itemIndex];
    }
    EXPECT_EQ(actual, 18);
}

TEST(BucketList, emptyBucketListHasEndOnBegin)
{
    BucketList<3> bl;
    EXPECT_EQ(*(std::end(bl)), *(std::begin(bl)));
}

TEST(BucketList, stdFindUsageOnBucketList)
{
    std::vector<int> values({18, 22, 43, 24});
    BucketList<3> bl;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        bl.add(itemIndex);

    for (std::size_t idx = 0; idx < values.size(); ++idx)
    {
        auto it = std::find(std::begin(bl), std::end(bl), idx);
        auto i  = *it;
        EXPECT_EQ(values[i], values[idx]);
    }

    bl.empty();
    bl.add(2);
    auto it = std::find(std::begin(bl), std::end(bl), 2);
    EXPECT_EQ(*it, 2);
    auto it2 = std::find(std::begin(bl), std::end(bl), 3);
    EXPECT_FALSE(it2 != std::end(bl));
}


TEST(BucketList, stdFindOnSingleItemBucket)
{
    BucketList<3> bl;
    bl.add(2);
    auto it = std::find(std::begin(bl), std::end(bl), 2);
    EXPECT_EQ(*it, 2);
}

TEST(BucketList, loopOverEmptyBucketMakesNoIteration)
{
    std::vector<int> values({1, 2, 3, 4, 5});
    BucketList<3> bl;
    auto cpt = 0;
    for (auto const& itemIndex : bl)
    {
        ++cpt;
    }
    EXPECT_EQ(0, cpt);
}




TEST(BucketList, trimRemovesTheNLastEmptyBuckets)
{
    std::vector<int> values({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
    BucketList<3> bl;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        bl.add(itemIndex);

    EXPECT_EQ(4 * 3, bl.capacity());
    EXPECT_EQ(bl.size(), bl.capacity());

    bl.empty();
    std::vector<int> others({1, 2, 3, 4});
    for (std::size_t itemIndex = 0; itemIndex < others.size(); ++itemIndex)
        bl.add(itemIndex);

    EXPECT_EQ(4, bl.size());
    EXPECT_EQ(12, bl.capacity());

    bl.trim(1);
    EXPECT_EQ(9, bl.capacity());


    bl.empty();
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        bl.add(itemIndex);
    EXPECT_EQ(12, bl.capacity());
    bl.trim(0);
    EXPECT_EQ(12, bl.capacity());
}



TEST(BucketList, removeAnElement)
{
    BucketList<3> bl;
    std::array<int, 4> values = {3, 4, 5, 6};

    bl.add(0);
    bl.add(1);
    bl.add(2);
    bl.add(3);
    bl.remove(1);
    std::array<int, 3> expected{3, 6, 5};
    std::size_t i = 0;
    for (auto index : bl)
    {
        EXPECT_EQ(values[index], expected[i++]);
    }
    EXPECT_EQ(bl.size(), expected.size());
    EXPECT_EQ(bl.capacity(), 3 * 2);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
