cmake_minimum_required(VERSION 3.16)
project(ranger_test)

## ======================================================================================##
## Check for C++11. For GCC this is >=4.7
## ======================================================================================##
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)

## ======================================================================================##
## Compiler flags
## ======================================================================================##
add_compile_options(-Wall)

## ======================================================================================##
## Dependencies
## ======================================================================================##
add_subdirectory(../ ranger_build)

## ======================================================================================##
## Google test
## ======================================================================================##
enable_testing()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# for CTest auto registration (gtest_add_tests)
include(GoogleTest)

## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(runUnitTests test_all.cpp utility_test.cpp)
target_link_libraries(runUnitTests
        PRIVATE
        ranger-lib
        GTest::gtest
        GTest::gtest_main)
gtest_add_tests(TARGET runUnitTests)

