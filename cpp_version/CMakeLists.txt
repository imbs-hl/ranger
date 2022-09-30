cmake_minimum_required(VERSION 3.2)

project(ranger)

include (GNUInstallDirs)

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
find_package (Threads)

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
## Library
## ======================================================================================##
file(GLOB_RECURSE SOURCES src/utility/*.cpp src/Forest/*.cpp src/Tree/*.cpp)
if (MSVC)
  list(APPEND SOURCES src/getopt/getopt.c)
endif ()
option(BUILD_SHARED_LIBS "shared/static lib" OFF)
add_library(ranger ${SOURCES})
if (MSVC)
  target_include_directories(ranger PUBLIC src/getopt)
endif ()
if (Threads_FOUND)
  target_link_libraries(ranger ${CMAKE_THREAD_LIBS_INIT})
endif ()
target_include_directories(ranger PUBLIC src/utility src/Forest src/Tree src)
install(TARGETS ranger RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
file(GLOB_RECURSE HEADERS src/*.h)
install (FILES ${HEADERS} DESTINATION include/ranger)

## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(ranger-cli src/main.cpp)
set_target_properties(ranger-cli PROPERTIES OUTPUT_NAME ranger)
target_link_libraries(ranger-cli ranger)
install(TARGETS ranger-cli DESTINATION ${CMAKE_INSTALL_BINDIR})

## ======================================================================================##
## Test
## ======================================================================================##
enable_testing()
add_test(NAME cli_help COMMAND ranger-cli --help)
