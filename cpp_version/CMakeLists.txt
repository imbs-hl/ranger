cmake_minimum_required(VERSION 3.1)
project(ranger)

## ======================================================================================##
## Check for C++14
## ======================================================================================##
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED TRUE)

## ======================================================================================##
## Compiler flags
## ======================================================================================##
add_compile_options(-Wall)

## ======================================================================================##
## Threads
## ======================================================================================##
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

## ======================================================================================##
## Subdirectories and source files
## ======================================================================================##
include_directories(src src/utility src/Forest src/Tree)
file(GLOB_RECURSE SOURCES src/*.cpp)

## ======================================================================================##
## Debug and release targets
## ======================================================================================##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )
ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(ranger ${SOURCES})
target_link_libraries(ranger Threads::Threads)

## ======================================================================================##
## Test
## ======================================================================================##
enable_testing()
add_test(NAME help COMMAND ranger --help)
