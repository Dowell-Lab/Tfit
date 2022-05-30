/**
 * @file test_Bedfile.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile class: \ref Bedfile 
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
   for (it = sut.setRegions.searchable.begin(); it != sut.setRegions.searchable.end(); it++) {
       chrcount++;
   }

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(sut.setRegions.chr_names.num_elements, chrcount);

   sut.setRegions.clearTrees();

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
   results = sut.setRegions.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.empty(), 1);   // Expects empty
   sut.setRegions.clearTrees();
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
   results = sut.setRegions.findOverlapIntervals(&query);

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 2);
   sut.setRegions.clearTrees();
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
   results = sut.setRegions.findOverlapIntervals(&query);

//   std::cout << sut.reportBedfileContents() << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(results.size(), 0);
   sut.setRegions.clearTrees();
}

TEST(Bedfile, load_Bedgraph) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/typical_region.bg"; 
   Bedgraph bg;

   // Act: call methods on SUT, capture output
   bg.load_file(file_name,0);   // Load without exiting bed

   // std::cout << bg.reportBedfileContents() << std::endl;
   std::vector<gInterval*> roi = bg.setRegions.regions[0];
   gInterval *output = roi[0];

   // Assert: Verify the outcome
   EXPECT_EQ(output->chromosome, "chr21");
   ASSERT_NE(output->data, nullptr);
   ASSERT_NE(output->data->cdata, nullptr);
   EXPECT_EQ(output->data->cdata->delta, 10);
}


