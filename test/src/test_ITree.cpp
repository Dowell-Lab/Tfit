/**
 * @file test_ITree.cpp 
 * @author Robin Dowell
 * @brief Testing the centered Interval Tree
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"
#include "ITree.h"
#include "Intervals.h"

// Creating a fixture with a simple example to test against.
class CITreeTest:public::testing::Test {
  protected:

   void SetUp() override {
    gInterval i1("chr1", 1, 10, "temp1");
    gInterval i2("chr1", 15, 24, "temp2");
    gInterval i5("chr1", 30, 40, "temp5");
    gInterval i3("chr1", 17, 29, "temp3");
    gInterval i4("chr1", 31, 38, "temp4");
    setofIntervals.push_back(&i1);
    setofIntervals.push_back(&i2);
    setofIntervals.push_back(&i5);
    setofIntervals.push_back(&i3);
    setofIntervals.push_back(&i4);
   }

   void TearDown() override { }
   std::vector<gInterval *> setofIntervals;
   // gInterval i1, i2, i3, i4, i5;
};

    /*
TEST_F(CITreeTest, treeConstruction)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *>setofIntervals;
    gInterval i1("chr1", 1, 10, "temp1");
    setofIntervals.push_back(&i1);
    gInterval i2("chr1", 15, 24, "temp2");
    setofIntervals.push_back(&i2);
    gInterval i5("chr1", 30, 40, "temp5");
    setofIntervals.push_back(&i5);
    gInterval i3("chr1", 17, 29, "temp3");
    setofIntervals.push_back(&i3);
    gInterval i4("chr1", 31, 38, "temp4");
    setofIntervals.push_back(&i4);

    std::string output = "\nL:\nL:\nIntervals overlapping 10.000000\n#chr1:1-10,temp1";
    output += "\nIntervals overlapping 24.000000\n#chr1:15-24,temp2\nIntervals overlapping ";
    output += "29.000000\n#chr1:17-29,temp3\nR:\nL:\nIntervals overlapping 38.000000\n";
    output += "#chr1:31-38,temp4\nIntervals overlapping 40.000000\n#chr1:30-40,temp5";

    // Act: call methods on SUT, capture output
    CITree sut(setofIntervals);

    // std::cout << sut.write_Full_Tree() << std::endl;

    // Assert: Verify the outcome
    EXPECT_THAT(output, sut.write_Full_Tree());
}
    */

// Need to test a known result, an edge case (edge of interval), and not in tree.
TEST(ITree, treeSearchcontains)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *>setofIntervals;
    gInterval i1("chr1", 1, 10, "temp1");
    setofIntervals.push_back(&i1);
    gInterval i2("chr1", 15, 24, "temp2");
    setofIntervals.push_back(&i2);
    gInterval i5("chr1", 30, 40, "temp5");
    setofIntervals.push_back(&i5);
    gInterval i3("chr1", 17, 29, "temp3");
    setofIntervals.push_back(&i3);
    gInterval i4("chr1", 31, 38, "temp4");
    setofIntervals.push_back(&i4);

    CITree sut(setofIntervals);
  
    // Act: call methods on SUT, capture output
    std::vector<gInterval *> results;
    results = sut.searchPoint((double)33);

    //std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
     //   std::cout << (*it)->write_out() << std::endl;
    //}

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST(ITree, treeSearchedge)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *>setofIntervals;
    gInterval i1("chr1", 1, 10, "temp1");
    setofIntervals.push_back(&i1);
    gInterval i2("chr1", 15, 24, "temp2");
    setofIntervals.push_back(&i2);
    gInterval i5("chr1", 30, 40, "temp5");
    setofIntervals.push_back(&i5);
    gInterval i3("chr1", 17, 29, "temp3");
    setofIntervals.push_back(&i3);
    gInterval i4("chr1", 31, 38, "temp4");
    setofIntervals.push_back(&i4);

    CITree sut(setofIntervals);
  
    // Act: call methods on SUT, capture output
    std::vector<gInterval *> results;
    results = sut.searchPoint((double)17);

    std::vector<gInterval *>::iterator it;
    for (it = results.begin(); it != results.end(); it++) {
        std::cout << (*it)->write_out() << std::endl;
    }

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 2);
}

TEST(ITree, treeSearchmissing)
{
    // Arrange: bring SUT to desired state
    std::vector<gInterval *>setofIntervals;
    gInterval i1("chr1", 1, 10, "temp1");
    setofIntervals.push_back(&i1);
    gInterval i2("chr1", 15, 24, "temp2");
    setofIntervals.push_back(&i2);
    gInterval i5("chr1", 30, 40, "temp5");
    setofIntervals.push_back(&i5);
    gInterval i3("chr1", 17, 29, "temp3");
    setofIntervals.push_back(&i3);
    gInterval i4("chr1", 31, 38, "temp4");
    setofIntervals.push_back(&i4);

    // Act: call methods on SUT, capture output
    CITree sut(setofIntervals);
  
    std::vector<gInterval *> results;
    results = sut.searchPoint((double)12);

    // std::vector<gInterval *>::iterator it;
    // for (it = results.begin(); it != results.end(); it++) {
      //  std::cout << (*it)->write_out() << std::endl;
    //}

    // Assert: Verify the outcome
    EXPECT_EQ(results.size(), 0);
}
