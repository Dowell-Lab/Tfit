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
    std::vector<gInterval *>setofIntervals;
    gInterval i1 = gInterval("chr1", 1, 10, "temp1");
    gInterval i2 = gInterval("chr1", 15, 24, "temp2");
    gInterval i3 = gInterval("chr1", 17, 29, "temp3");
    gInterval i4 = gInterval("chr1", 31, 38, "temp4");
    gInterval i5 = gInterval("chr1", 30, 40, "temp5");
    CITree sut;
};

TEST_F(ITreeTest, treeConstruction)
{
    std::string output = "\nL:\nL:\nIntervals overlapping 10.000000\n#chr1:1-10,temp1";
    output += "\nIntervals overlapping 24.000000\n#chr1:15-24,temp2\nIntervals overlapping ";
    output += "29.000000\n#chr1:17-29,temp3\nR:\nL:\nIntervals overlapping 38.000000\n";
    output += "#chr1:31-38,temp4\nIntervals overlapping 40.000000\n#chr1:30-40,temp5";

    // Assert: Verify the outcome
    EXPECT_EQ(output, sut.write_Full_Tree());
}

// Need to test a known result, an edge case (edge of interval), and not in tree.
TEST_F(ITreeTest, SearchPointContains)
{
    // Arrange
    std::vector<gInterval *> results;
    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)33);

    //std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
     //   std::cout << (*it)->write_out() << std::endl;
    //}

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchPointEdge)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *> results;
  
    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)17);

    //  std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
    //     std::cout << (*it)->write_out() << std::endl;
    // }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchPointMissing)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *> results;

    // Act: call methods on SUT, capture output
    results = sut.searchPoint((double)12);

    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
      //  std::cout << (*it)->write_out() << std::endl;
    //}

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

TEST_F(ITreeTest, SearchIntervalExact) 
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&i4);

    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
    //   std::cout << (*it)->write_out() << std::endl;
    // }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST_F(ITreeTest, SearchIntervalMissing) 
{
    // Arrange: bring SUT to desired state
    gInterval queryI("chr1", 11, 14, "missing");
    std::vector<gInterval *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&queryI);

    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
    //   std::cout << (*it)->write_out() << std::endl;
    // }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}

TEST_F(ITreeTest, SearchIntervalContained) 
{
    // Arrange: bring SUT to desired state
    gInterval queryI("chr1", 33, 35, "contained");
    std::vector<gInterval *> results;

    // Act: call methods on SUT, capture output
    results = sut.overlapSearch(&queryI);

    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
    //   std::cout << (*it)->write_out() << std::endl;
    // }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

