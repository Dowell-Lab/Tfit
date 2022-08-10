/**
 * @file test_Bedfile.cpp 
 * @author Robin Dowell
 * @brief Testing the Bedfile and Bedgraph classes: \ref Bedfile Bedgraph
 * @date 2022-03-18
 * 
 */
#include "gmock/gmock.h"

#include "Intervals.h"
#include "ITree.h"
#include "Bedfile.h"

// Bedfile:  Test I/O for bed3, bed4, bed6, bed12 ?

/*
TEST(Bedfile, ReadWriteEquivalent) 
{

}
*/

TEST(Bedgraph, ReadErrorCheck) 
{
    std::string filename = "../examples/typical_region.bg";    //chr21 33401693 33407411	
    Bedgraph sut;

    sut.load_file(filename, false);

    sut.reportBedfileContents();

    EXPECT_TRUE(true);
}