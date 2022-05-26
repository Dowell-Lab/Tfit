/**
 * @file reproduce_Bidirectional.cpp 
 * @author Robin Dowell
 * @brief Testing the bidirectional model class reproduces the old code
 * @date 2022-05-26
 * 
 */
#include "gmock/gmock.h"
#include <fstream>
#include "Models.h"
#include "model.h"

TEST(BidirEQ, reproducePDF)
{
    // Arrange: bring SUT to desired state
    Bidirectional sut(35., 0.5, 0.25, 0.5, 2.);
    EMG oldmethod(35., 0.5, 0.25, 1, 0.5);

    // Act: call methods on SUT, capture output
    double newresult = sut.pdf(37.,'+');
    double oldresult = oldmethod.pdf(37., 1);

    std::cout << to_string(oldresult) + " is now " 
          + to_string(newresult) << std::endl;

    // Assert: Verify the outcome
    ASSERT_LE(abs(oldresult - newresult), 0.0001);  
}
