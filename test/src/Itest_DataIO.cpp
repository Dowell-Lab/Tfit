/**
 * @file Itest_DataIO.cpp 
 * @author Robin Dowell
 * @brief Integration of loading bed and bedGraphs with 
 * searching regions for specific intervals.
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"

#include "Bed.h"
#include "ITree.h"
#include "Bedfile.h"

// Bedfile:  I/O for bed4, bed6, bed12

TEST(ITest_DataIO, multiChromCheck) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   
   // Act: call methods on SUT, capture output
   sut.load_file(file_name);

   int chrcount = 0;
   std::map<int, CITree *>::iterator it;
   for (it = sut.setRegions.searchable.begin(); it != sut.setRegions.searchable.end(); it++) {
       chrcount++;
   }

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(sut.setRegions.chr_names.num_elements, chrcount);

   sut.setRegions.clearTrees();

}

TEST(ITest_DataIO, checkSearch_noChr) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<bed12 *> results;
   bed12 query("chr51", 19731000, 19768000, "example");

   // Act: call methods on SUT, capture output
   results = sut.setRegions.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.empty(), 1);   // Expects empty
   sut.setRegions.clearTrees();
}

TEST(ITest_DataIO, checkSearch_hasoverlap) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<bed12 *> results;
   bed12 query("chr22", 19731000, 19768000, "example");

   // std::cout << sut.reportBedfileContents() << std::endl;
   // Act: call methods on SUT, capture output
   results = sut.setRegions.findOverlapIntervals(&query);

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 2);
   sut.setRegions.clearTrees();
}

TEST(ITest_DataIO, checkSearch_noOverlap) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<bed12 *> results;
   bed12 query("chr21", 19731000, 19768000, "example");

   // Act: call methods on SUT, capture output
   results = sut.setRegions.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 0);
   sut.setRegions.clearTrees();
}


