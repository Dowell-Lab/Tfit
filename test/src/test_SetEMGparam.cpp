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
    EXPECT_EQ(sut.collection[1]->mu, (double)57974849);
    EXPECT_EQ(sut.collection[0]->omega, (double)0.617766);
    EXPECT_EQ(sut.collection[1]->omega, (double)0.310201);

    // Teardown
}

TEST(Set_EMGparam, write)
{
    // Arrange
    Set_EMGparameters sut(2);
    sut.collection[0] = new EMGparameters(45333,65,44,0.345,32,0.5);
    sut.collection[1] = new EMGparameters(11111,35,22,0.081,11,0.3);

    // Act
    std::string Toutput = sut.write();
    std::string output;
    output  = (std::string) "~2,0.000000\t" + "45333.000000,11111.000000\t" 
        + "65.000000,35.000000\t" + "44.000000,22.000000\t" "0.345000,0.081000\t"
        + "32.000000,11.000000\t" + "0.500000,0.300000";

    // Assert
    EXPECT_TRUE(Toutput == output);

    // Teardown
}
