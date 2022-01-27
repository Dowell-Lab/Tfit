cmake_minimum_required(VERSION 3.12)
project(Tfit)

option(BUILD_TEST "Enable Testing code" ON)

if(BUILD_TEST) 
    add_subdirectory(test)
endif()

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -std=c++11 -fopenmp -g -D_LGIBCXX_USE_CXX1_ABI=0 ")

#Bring the headers, such as Student.h into the project
include_directories(src)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")

# This is a temporary testing stub.
# message(${SOURCES})
list(REMOVE_ITEM SOURCES "/home/robindowell/src/Tfit/Dowell/Tfit/src/arg_main.cpp")

add_executable(Tfit ${SOURCES})

target_link_libraries(Tfit ${MPI_CXX_LIBRARIES})
