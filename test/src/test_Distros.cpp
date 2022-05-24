/**
 * @file test_Distros.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-05-24
 * 
 */
#include "gmock/gmock.h"
#include <fstream>
#include "Distro.h"

TEST(Distros, standardNormal)
{
    // Arrange: bring SUT to desired state
    Normal standard;
    int num_samples = 100000;

    // Act: call methods on SUT, capture output
    std::vector<double> output = standard.generate_data(num_samples);

    // Calculate sample mean
    double total = 0;
    for (int i = 0; i < num_samples; i++) {
       total += output[i];
    }
    double mean = total/num_samples;

    // Calculate sample standard deviation 
    total = 0;
    for (int i = 0; i < num_samples; i++) {
      total += pow((output[i] - mean),2);
    }
    double std = sqrt(total /(num_samples-1));

    // Assert: Verify the outcome
    EXPECT_EQ(standard.write_out(), "N(0.00,1.00)");
    EXPECT_LE(abs(mean - 0.00), 0.01);  
    EXPECT_LE(abs(std - 1.0), 0.01);  
}

TEST(Distros, Normal)
{
    // Arrange: bring SUT to desired state
    Normal norm(5,2);
    int num_samples = 100000;

    // Act: call methods on SUT, capture output
    std::vector<double> output = norm.generate_data(num_samples);

    // Calculate sample mean
    double total = 0;
    for (int i = 0; i < num_samples; i++) {
       total += output[i];
    }
    double mean = total/num_samples;

    // Calculate sample standard deviation 
    total = 0;
    for (int i = 0; i < num_samples; i++) {
      total += pow((output[i] - mean),2);
    }
    double std = sqrt(total /(num_samples-1));

    // Assert: Verify the outcome
    EXPECT_EQ(norm.write_out(), "N(5.00,2.00)");
    EXPECT_LE(abs(mean - 5.00), 0.01);  
    EXPECT_LE(abs(std - 2.0), 0.01);  
    EXPECT_LE(abs(norm.pdf(3) - 0.1209854), 0.0001);  // R: dnorm(3,mean=5,sd=3)
    EXPECT_LE(abs(norm.cdf(3) - 0.1586553), 0.0001);  // R: pnorm(3,mean=5,sd=3)
}

TEST(Distros, Exponential)
{
    // Arrange: bring SUT to desired state
    Exponential exp(5);
    int num_samples = 100000;

    // Act: call methods on SUT, capture output
    std::vector<double> output = exp.generate_data(num_samples);

    // Calculate sample mean
    double total = 0;
    for (int i = 0; i < num_samples; i++) {
       total += output[i];
    }
    double mean = total/num_samples;

    // Assert: Verify the outcome
    EXPECT_EQ(exp.write_out(), "E(5.00)");
    EXPECT_LE(abs(mean - (1/exp.lambda)), 0.01);  
    EXPECT_LE(abs(exp.pdf(0.2) - 1.839397), 0.0001);  // R: dexp(0.2,rate=5);
    EXPECT_LE(abs(exp.cdf(0.2) - 0.6321206), 0.0001);  // R: pexp(0.2,rate=5)
}

TEST(Distros, Uniform)
{
    // Arrange: bring SUT to desired state
    Uniform uni(5,10);
    int num_samples = 100000;

    // Act: call methods on SUT, capture output
    std::vector<double> output = uni.generate_data(num_samples);

    // Calculate sample mean
    double total = 0;
    for (int i = 0; i < num_samples; i++) {
       total += output[i];
    }
    double mean = total/num_samples;

    // Assert: Verify the outcome
    EXPECT_EQ(uni.write_out(), "U(5.00,10.00)");
    EXPECT_LE(abs(mean - 0.5*(uni.lower+uni.upper)), 0.01);  
    EXPECT_LE(abs(uni.pdf(6) - 0.2), 0.0001);  // R: dunif(6,min=5,max=10);
    EXPECT_LE(abs(uni.cdf(7) - 0.4), 0.0001);  // R: punif(7,min=5,max=10);
}
