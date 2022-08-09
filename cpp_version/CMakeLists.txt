cmake_minimum_required(VERSION 3.16)
project(ranger VERSION 0.14.1)

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
## Generate Version Header
## ======================================================================================##
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/global/version.h @ONLY)

## ======================================================================================##
## Default to Release
## ======================================================================================##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

## ======================================================================================##
# Dependencies
## ======================================================================================##
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

## ======================================================================================##
## Workaround for clang on windows
## ======================================================================================##
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(WIN32_CLANG TRUE CACHE BOOL "Clang on windows does not ship getopt.h, so we must provide it")
  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/getopt/getopt.h")
    file(DOWNLOAD "https://raw.githubusercontent.com/alex85k/wingetopt/master/src/getopt.h"
            "${CMAKE_BINARY_DIR}/getopt/getopt.h"
            TLS_VERIFY ON
            EXPECTED_HASH SHA256=87534c9d7ace0230c021f1154cfb8806f70ece100b7f3dd9e0913292fdcd8531)
  endif()
  if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/getopt/getopt.c")
    file(DOWNLOAD "https://raw.githubusercontent.com/alex85k/wingetopt/master/src/getopt.c"
            "${CMAKE_BINARY_DIR}/getopt/getopt.c"
            TLS_VERIFY ON
            EXPECTED_HASH SHA256=2d0f6ef53db31f0bfa88c949ec57b4e869110b38c89ff7661927dd6e36a9ff71)
  endif()

  add_library(getopt STATIC "${CMAKE_BINARY_DIR}/getopt/getopt.c")
  target_include_directories(getopt PUBLIC
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/getopt>
          )
endif()

## ======================================================================================##
## Library
## ======================================================================================##
add_library(ranger-lib STATIC
        src/Tree/TreeRegression.cpp
        src/Tree/TreeClassification.cpp
        src/Tree/Tree.cpp
        src/Tree/TreeProbability.cpp
        src/Tree/TreeSurvival.cpp
        src/utility/utility.cpp
        src/utility/Data.cpp
        src/utility/ArgumentHandler.cpp
        src/Forest/ForestSurvival.cpp
        src/Forest/ForestRegression.cpp
        src/Forest/Forest.cpp
        src/Forest/ForestProbability.cpp
        src/Forest/ForestClassification.cpp
        )

target_include_directories(ranger-lib PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/global>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Forest>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Tree>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/utility>
        )

target_link_libraries(ranger-lib
        PRIVATE
        Threads::Threads
        $<$<BOOL:${WIN32_CLANG}>:getopt> # getopt.h is not available on Windows with clang.
        )

## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(ranger src/main.cpp)

target_link_libraries(ranger
        PRIVATE
        ranger-lib)

if(WIN32 AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(ranger PRIVATE -static)
endif()
