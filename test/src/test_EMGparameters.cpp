/**
 * @file test_EMGparameters.cpp
 * @author Robin Dowell
 * @brief Unit Testing: testing Tfit/src/single_model.cpp \ref EMGparameters
 * @version 0.1
 * @date 2022-02-10
 * 
 */
#include "gmock/gmock.h"
#include "single_model.h"


TEST(EMGparameters, ReadWrite)
{
    // Arrange
    EMGparameters sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3);
    EMGparameters input;

    // Act
    std::string outString = sut.write();
    input.read(outString);

    // Assert
    EXPECT_EQ(sut.mu, (double)35000.0); // Checks basic I/O
    EXPECT_EQ(sut.mu,input.mu);  
    EXPECT_EQ(sut.sigma,input.sigma);  
    EXPECT_EQ(sut.lambda,input.lambda);  
    EXPECT_EQ(sut.pi,input.pi);  
    EXPECT_EQ(sut.footprint,input.footprint);  
    EXPECT_EQ(sut.omega,input.omega);  


}

TEST(EMGparameters, Start_Stop)
{
    // Arrange
    EMGparameters sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3);

    // Assert
    EXPECT_EQ(sut.getStart(), (double)34605);  
    EXPECT_EQ(sut.getEnd(), (double)35395);  
}

TEST(EMGparameters, fetch_asStrings)
{
    // Arrange
    EMGparameters sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3);

    // Act
    std::vector<std::string> params = sut.fetch_as_strings();

    // Assert
    EXPECT_EQ(params.size(), 6);

    /* Note: to_string with doubles can be a problem in 
    terms of things matching afterwards.  I side step that here
    by using the same to_string call as the fetch_as_strings. */
    EXPECT_EQ(params[0],std::to_string(35000.0));  
    EXPECT_EQ(params[1],std::to_string(350.0));  
    EXPECT_EQ(params[2],std::to_string(45.0));  
    EXPECT_EQ(params[3],std::to_string(0.8));  
    EXPECT_EQ(params[4],std::to_string(30.0));  
    EXPECT_EQ(params[5],std::to_string(0.3));  

    // Teardown
    params.clear();
}
