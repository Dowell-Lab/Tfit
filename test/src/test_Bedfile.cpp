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
#include <sstream>

#include "Intervals.h"
#include "ITree.h"
#include "Bedfile.h"

// Bedfile:  Test I/O for bed3, bed4, bed6, bed12 ?

/*
TEST(Bedfile, ReadWriteEquivalent) 
{

}
*/

TEST(Bedfile, ReadWriteCheck) 
{
    std::string filename = "../examples/multiple.bed";    //chr21 33401693 33407411	
    Bedfile sut;
    sut.load_file(filename);

    // Read the original file into a buffer called original
    std::ostringstream original; 
    std::ifstream input(filename.c_str()); 
    original << input.rdbuf(); 

    std::string outfn = "../examples/multiple_repro.bed";    //chr21 33401693 33407411	
    sut.write_file(outfn);

    // Read the wrote file into a buffer called duplicate 
    std::ostringstream duplicate; 
    std::ifstream second(outfn.c_str()); 
    duplicate << second.rdbuf(); 

    EXPECT_EQ(original.str(), duplicate.str());
}