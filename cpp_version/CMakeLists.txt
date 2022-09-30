cmake_minimum_required(VERSION 3.2)

project(ranger)

## ======================================================================================##
## Check for C++14. For GCC this is >=4.7
## ======================================================================================##
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
  message("Compiler with C++14 support found.")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler, i.e. gcc >= 5 or Clang >= 3.4.")
endif()

## ======================================================================================##
## Compiler flags
## ======================================================================================##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")

## ======================================================================================##
## In Clang phtread flag only for compiler, not for linker. For
## windows use static linking
## ======================================================================================##
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS -pthread)
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -static")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

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

