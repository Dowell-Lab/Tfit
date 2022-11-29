/**
 * @file test_ModelParamSet.cpp
 * @author Robin Dowell
 * @brief Unit Testing: Both the set of parameters associated with a 
 * single model fit and a container of them. 
 * @version 0.1
 * @date 2022-02-10
 * 
 */
#include "gmock/gmock.h"
#include "ModelParamSet.h"

TEST(ModelParams, ReadWrite)
{
    // Arrange
    ModelParams sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3, 0.2, 0.1, 34000.0, 38000.0);
    ModelParams input;

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
    EXPECT_EQ(sut.omega[0],input.omega[0]);  
    EXPECT_EQ(sut.omega[1],input.omega[1]);  
    EXPECT_EQ(sut.omega[2],input.omega[2]);  
}

TEST(ModelParams, ReadWriteJSON)
{
    // Arrange
    ModelParams sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3, 0.2, 0.1, 34000.0, 38000.0);
    ModelParams input;

    // Act
    std::string outString = sut.writeAsJSON();
    // std::cout << outString + "\n" << std::endl;
    input.readFromJSON(outString);

    // Assert
    EXPECT_EQ(sut.mu, (double)35000.0); // Checks basic I/O
    EXPECT_EQ(sut.mu,input.mu);  
    EXPECT_EQ(sut.sigma,input.sigma);  
    EXPECT_EQ(sut.lambda,input.lambda);  
    EXPECT_EQ(sut.pi,input.pi);  
    EXPECT_EQ(sut.footprint,input.footprint);  
    EXPECT_EQ(sut.omega[0],input.omega[0]);  
    EXPECT_EQ(sut.omega[1],input.omega[1]);  
    EXPECT_EQ(sut.omega[2],input.omega[2]);  
}

TEST(ModelParams, Start_Stop)
{
    // Arrange
    ModelParams sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3, 0.2, 0.1, 34000.0, 37000.0);

    // Assert
    EXPECT_EQ(sut.getBedStart(), (double)34605);  
    EXPECT_EQ(sut.getBedEnd(), (double)35395);  
}

TEST(ModelParams, fetch_asStrings)
{
    // Arrange
    ModelParams sut(35000.0, 350.0, 45.0, 0.8, 30.0, 0.3, 0.2, 0.1, 34000.0, 39000.0);

    // Act
    std::vector<std::string> params = sut.fetch_as_strings();

    // Assert
    EXPECT_EQ(params.size(), 10);

    /* Note: to_string with doubles can be a problem in 
    terms of things matching afterwards.  I side step that here
    by using the same to_string call as the fetch_as_strings. */
    EXPECT_EQ(params[0],std::to_string(35000.0));  
    EXPECT_EQ(params[1],std::to_string(350.0));  
    EXPECT_EQ(params[2],std::to_string(45.0));  
    EXPECT_EQ(params[3],std::to_string(0.8));  
    EXPECT_EQ(params[4],std::to_string(30.0));  
    EXPECT_EQ(params[5],std::to_string(0.3));  
    EXPECT_EQ(params[6],std::to_string(0.2));  
    EXPECT_EQ(params[7],std::to_string(0.1));  

    // Teardown
    params.clear();
}

TEST(ModelParamSet, Set_read_from_K_models)
{
    // Arrange
    std::string line  = (std::string)"~2,-18405.714702" + "\t" + "57970016,57974849" 
        + "\t" + "69.438995,37.975287" 
        + "\t" + "371.734252,408.579080" + "\t"	+ "0.486291,0.586718"
        + "\t" + "121.206766,250.000000" + "\t" + "0.617766,0.069621,0.000284|0.310201,0.000341,0.000348"
        + "\t" + "57976184.000000,57976184.000000" + "\t" + "57965508.000000,57965508.000000";
    ModelParamSet sut;

    // Act
    sut.readFromKmodels(line);

    // Assert
    EXPECT_EQ(sut.collection[0]->mu, (double)57970016);
    EXPECT_EQ(sut.collection[1]->mu, (double)57974849);
    EXPECT_EQ(sut.collection[0]->omega[0], (double)0.617766);
    EXPECT_EQ(sut.collection[0]->omega[1], (double)0.069621);
    EXPECT_EQ(sut.collection[1]->omega[0], (double)0.310201);
    EXPECT_EQ(sut.collection[1]->omega[1], (double)0.000341);

    // Teardown
}

TEST(ModelParamSet, Set_write)
{
    // Arrange
    ModelParamSet sut(2);
    sut.collection[0] = new ModelParams(45333,65,44,0.345,32,0.5, 0.2, 0.1, 44000, 48000);
    sut.collection[1] = new ModelParams(11111,35,22,0.081,11,0.3, 0.2, 0.1, 10000, 14000);

    // Act
    std::string Toutput = sut.writeAsKmodels();
    std::cout << Toutput << std::endl;
    std::string output;


    output  = (std::string) "~2,0.000000\t45333.000000,11111.000000\t"
        + "65.000000,35.000000\t44.000000,22.000000\t0.345000,0.081000\t"
        + "32.000000,11.000000\t0.500000,0.200000,0.100000|0.300000,0.200000,0.100000\t"
        + "44000.000000,10000.000000\t48000.000000,14000.000000";

    // Assert
    EXPECT_TRUE(Toutput == output);

    // Teardown
}
