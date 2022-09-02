/**
 * @file test_Bedfile.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile and Bedgraph classes: \ref Bedfile Bedgraph
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"

#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>

#include "Intervals.h"
#include "ITree.h"
#include "Bedfile.h"

class BedfileIOTest: public :: testing::Test {
  protected:
  void SetUp() override { }
  void RunIO(std::string filename) {
    // Read in the file (both to bedfile object and string stream)
    sut.load_file(filename);
    std::ifstream input(filename.c_str()); 
    original << input.rdbuf(); 

    // Now write the file back out and read it that to a second string buffer. 
    outfn = filename + "_repro";
    sut.write_file(outfn);
    std::ifstream second(outfn.c_str()); 
    duplicate << second.rdbuf(); 
  }
  void TearDown() override {
    std::remove(outfn.c_str()); // Delete the temp file made.
  }
  Bedfile sut;
  std::ostringstream original; 
  std::string outfn;
  std::ostringstream duplicate; 

};

/**
 * @brief Do I write back out the same file I read in.
 * Caveats: doesn't keep comments and writes in chromosome sorted order,
 * so the input file MUST behave under these conditions!
 */
TEST_F(BedfileIOTest, ReadWriteBed4) 
{
    // Arrange
    std::string filename = "../examples/multiple.nocomments.bed";    //BED4
    // Act
    RunIO(filename);
    //Assert
    EXPECT_EQ(original.str(), duplicate.str());
}

TEST_F(BedfileIOTest, ReadWriteBed6) 
{
    // Arrange
    std::string filename = "../examples/multiple.bed6";    //BED6
    // Act
    RunIO(filename);
    //Assert
    EXPECT_EQ(original.str(), duplicate.str());
}

TEST_F(BedfileIOTest, ReadWriteBed12) 
{
    // Arrange
    std::string filename = "../examples/multiple.bed12";    //BED6
    // Act
    RunIO(filename);
    //Assert
    EXPECT_EQ(original.str(), duplicate.str());
}
