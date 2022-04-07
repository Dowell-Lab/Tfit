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

TEST(Bedfile, fromBedFile) 
{
    // Arrange: bring SUT to desired state
   std::string file_name = "../examples/multiple.bed"; 
   Bedfile sut(file_name);
   
   // Act: call methods on SUT, capture output
   sut.load_file();

   // std:: cout << sut.print_tree_at_chromosome((std::string)"chr22") << std::endl;

   // Assert: Verify the outcome
   EXPECT_EQ(sut.chr_names.num_elements, 2);
}
