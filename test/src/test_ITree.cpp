/**
 * @file test_ITree.cpp 
 * @author Robin Dowell
 * @brief Testing the centered Interval Tree: \ref CITree 
 * Notice that the testing of \ref Inode is only implicit
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"
#include "ITree.h"
#include "Intervals.h"

class ITreeTest: public :: testing::Test {
    protected:
    void SetUp() override {
      setofIntervals.push_back(&i1);
      setofIntervals.push_back(&i2);
      setofIntervals.push_back(&i5);
      setofIntervals.push_back(&i3);
      setofIntervals.push_back(&i4);
      sut.constructTree(setofIntervals);
    }
    void TearDown() override {
      sut.destroyTree();
    }
    // Arrange: bring SUT to desired state
    std::vector<bed12 *>setofIntervals;
    bed12 i1 = bed12("chr1", 1, 10, "temp1");
    bed12 i2 = bed12("chr1", 15, 24, "temp2");
    bed12 i3 = bed12("chr1", 17, 29, "temp3");
    bed12 i4 = bed12("chr1", 31, 38, "temp4");
    bed12 i5 = bed12("chr1", 30, 40, "temp5");
    CITree sut;
};

TEST_F(ITreeTest, treeConstruction)
{
    std::string output = "\nL:\nL:\nIntervals overlapping 10.000000\n";
    output += "#chr1:1-10,temp1,0.0000,.,None\nIntervals overlapping ";
    output += "24.000000\n#chr1:15-24,temp2,0.0000,.,None\nIntervals ";
    output += "overlapping 29.000000\n#chr1:17-29,temp3,0.0000,.,None\n";
    output += "R:\nL:\nIntervals overlapping 38.000000\n#chr1:31-38,temp4";
    output += ",0.0000,.,None\nIntervals overlapping 40.000000\n#chr1:30-";
    output += "40,temp5,0.0000,.,None";

    
    // Assert: Verify the outcome
    EXPECT_EQ(output, sut.write_Full_Tree());
}

TEST_F(ITreeTest, SearchPointContains)
{
    // Arrange
    std::vector<bed12 *> results;
    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)33);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchPointEdge)
{
    // Arrange: bring SUT to desired state
    std::vector<bed12 *> results;
  
    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)17);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchPointMissing)
{
    // Arrange: bring SUT to desired state
    std::vector<bed12 *> results;

    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)12);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

TEST_F(ITreeTest, SearchIntervalExact) 
{
    // Arrange: bring SUT to desired state
    std::vector<bed12 *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&i4);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchIntervalMissing) 
{
    // Arrange: bring SUT to desired state
    bed12 queryI("chr1", 11, 14, "missing");
    std::vector<bed12 *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&queryI);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

TEST_F(ITreeTest, SearchIntervalContained) 
{
    // Arrange: bring SUT to desired state
    bed12 queryI("chr1", 33, 35, "contained");
    std::vector<bed12 *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&queryI);

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

