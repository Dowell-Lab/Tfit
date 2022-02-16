/**
 * @file test_SetEMGparam.cpp
 * @author Robin Dowell
 * @brief Unit Testing: testing Tfit/src/single_model.cpp
 * @version 0.1
 * @date 2022-02-16
 * 
 */
#include "gmock/gmock.h"
#include "single_model.h"

TEST(Set_EMGparam, read_from_K_models)
{
    // Arrange
    std::string line  = (std::string)"~2,-18405.714702" + "\t" + "57970016,57974849" 
        + "\t" + "69.438995,37.975287" 
        + "\t" + "371.734252,408.579080" + "\t"	+ "0.486291,0.586718"
        + "\t" + "121.206766,250.000000" + "\t" + "0.617766,0.069621,0.000284|0.310201,0.000341,0.000348"
        + "\t" + "57976184.000000,57976184.000000" + "\t" + "57965508.000000,57965508.000000";
    Set_EMGparameters sut;

    // Act
    sut.read_from_K_models(line);

    // Assert
    EXPECT_EQ(sut.collection[0]->mu, (double)57970016);

    // Teardown
}
