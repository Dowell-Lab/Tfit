/**
 * @file test_main.cpp
 * @author Robin Dowell
 * @brief Unit Testing Main 
 * @version 0.1
 * @date 2022-01-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "gmock/gmock.h"

int main(int argc, char** argv) {
	::testing::InitGoogleMock(&argc, argv);
	return RUN_ALL_TESTS();
}

