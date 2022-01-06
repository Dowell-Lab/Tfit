/**
 * @file test_split.cpp
 * @author Robin Dowell
 * @brief Unit Testing: testing Tfit/src/split.cpp
 * @version 0.1
 * @date 2022-01-06
 * 
 */
#include "gmock/gmock.h"
#include "split.h"

using namespace std;

TEST(Split, TestCaseOne)
{
    // Arrange
    std::string testString;
    testString = "Blanket\tPillow";
    vector<string> answer;
    // Act
    answer = string_split(testString, '\t');
    // Assert
    EXPECT_THAT(answer.back(), "Pillow");
    EXPECT_EQ(answer.size(),2);
}
