# Description : Definition of all source and header files of one package.
#               To add this package to a project, include this file in the corresponding CMakeLists.txt file.
#      Author : Robin Dowell
# ================================================================================================
cmake_minimum_required (VERSION 3.7)
project(TestTfit)

# add definitions
add_definitions(-fprofile-arcs -ftest-coverage)

find_package(MPI REQUIRED)

# Include diretories
include_directories(/usr/include)
link_directories(/usr/lib)
#Bring the headers, such as Student.h into the project
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Set the compiler options
##  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -O0 -Wall -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -std=c++11 -fopenmp -g -D_LGIBCXX_USE_CXX1_ABI=0 -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)

# Create OBJECT_DIR variable
set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/TestTfit.dir/src)
message("-- Object files will be output to: ${OBJECT_DIR}")

# Find required packages
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Set the sources
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB sources "../src/*.cpp")
# Remove other main functions
list(REMOVE_ITEM sources "/home/robindowell/src/Tfit/Dowell/Tfit/test/../src/arg_main.cpp")
list(REMOVE_ITEM sources "/home/robindowell/src/Tfit/Dowell/Tfit/test/../src/main.cpp")
list(REMOVE_ITEM sources "/home/robindowell/src/Tfit/Dowell/Tfit/test/../src/MainTesting.cpp")
list(APPEND sources
                src/explore_coreEM.cpp
                src/reproduce_Bidirectional.cpp
                src/reproduce_loadBedGraph.cpp
                src/reproduce_Nested.cpp
                src/reproduce_segment.cpp
                src/reproduce_templatematching.cpp
                src/test_main.cpp
                src/test_Bedfile.cpp
                src/test_Data.cpp
                src/test_Distros.cpp
                src/test_helper.cpp
                src/test_Intervals.cpp
                src/test_ITree.cpp
                src/test_ModelParamSet.cpp
                src/test_Models.cpp
                src/test_ModelSupport.cpp
                src/test_Regions.cpp
                src/test_split.cpp
                )

# Set Include directories
include_directories(
                src/
                ../src/
                )

# Create the executable
add_executable(TestTfit ${sources})

# Create the gcov target. Run coverage tests with 'make gcov'
add_custom_target(gcov
    COMMAND mkdir -p gcoverage
    COMMAND ${CMAKE_MAKE_PROGRAM} test
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
add_custom_command(TARGET gcov
    COMMAND echo "=================== GCOV ===================="
    COMMAND gcov -b ${CMAKE_SOURCE_DIR}/src/*.cpp -o ${OBJECT_DIR}
    COMMAND echo "-- Source diretorie: ${CMAKE_SOURCE_DIR}/src/"
    COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/gcoverage"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gcoverage
    )
add_dependencies(gcov TestTfit)

# Create the lcov target. Run result tests with 'make lcov'
add_custom_target(lcov
    COMMAND mkdir -p lcoverage
    )
add_custom_command(TARGET lcov
    COMMAND echo "=================== LCOV ===================="
    COMMAND echo "-- Passing lcov tool under code coverage"
    COMMAND lcov --capture --directory ../ --output-file lcoverage/main_coverage.info
    COMMAND echo "-- Generating HTML output files"
    COMMAND genhtml lcoverage/main_coverage.info --output-directory lcoverage
    )

# Make sure to clean up the coverage folder
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gcoverage)

# Create the gcov-clean target. This cleans the build as well as generated 
# .gcda and .gcno files.
add_custom_target(init
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND rm -f ${OBJECT_DIR}/*.gcno
    COMMAND rm -f ${OBJECT_DIR}/*.gcda
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

# Testing
enable_testing()

add_test(output_test ${CMAKE_CURRENT_BINARY_DIR}/TestTfit)

# link Boost libraries to the new target
target_link_libraries(TestTfit gtest)
target_link_libraries(TestTfit gmock)
target_link_libraries(TestTfit pthread)
target_link_libraries(TestTfit -fprofile-arcs)
target_link_libraries(TestTfit ${MPI_CXX_LIBRARIES})
#target_link_libraries(TestTfit gcov)

