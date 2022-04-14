/**
 * @file test_Bedfile.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile class (reading in a bedfile) 
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"

#include "Intervals.h"
#include "ITree.h"
#include "Bedfile.h"

TEST(Bedfile, multiChromCheck) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   
   // Act: call methods on SUT, capture output
   sut.load_file(file_name);

   int chrcount = 0;
   std::map<int, CITree *>::iterator it;
   for (it = sut.intervals.begin(); it != sut.intervals.end(); it++) {
       chrcount++;
   }

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(sut.chr_names.num_elements, chrcount);
}

TEST(Bedfile, checkSearch_noChr) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<gInterval *> results;
   gInterval query("chr51", 19731000, 19768000, "example");

   // Act: call methods on SUT, capture output
   results = sut.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.empty(), 1);   // Expects empty
}

TEST(Bedfile, checkSearch_hasoverlap) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<gInterval *> results;
   gInterval query("chr22", 19731000, 19768000, "example");

   // std::cout << sut.reportBedfileContents() << std::endl;
   // Act: call methods on SUT, capture output
   results = sut.findOverlapIntervals(&query);

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 2);
}

TEST(Bedfile, checkSearch_noOverlap) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut;
   sut.load_file(file_name);

   std::vector<gInterval *> results;
   gInterval query("chr21", 19731000, 19768000, "example");

   // Act: call methods on SUT, capture output
   results = sut.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 0);
}