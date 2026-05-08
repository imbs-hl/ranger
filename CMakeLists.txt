cmake_minimum_required(VERSION 3.23)

project(ranger
    DESCRIPTION "A Fast Implementation of Random Forests"
    HOMEPAGE_URL "https://github.com/imbs-hl/ranger"
    LANGUAGES C CXX
    VERSION 0.18.0)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_library(ranger)
add_library(Ranger::ranger ALIAS ranger)

target_compile_features(ranger PUBLIC cxx_std_14)

target_sources(ranger
  PRIVATE
    src/Data.cpp
    src/Forest.cpp
    src/ForestClassification.cpp
    src/ForestProbability.cpp
    src/ForestRegression.cpp
    src/ForestSurvival.cpp
    src/Tree.cpp
    src/TreeClassification.cpp
    src/TreeProbability.cpp
    src/TreeRegression.cpp
    src/TreeSurvival.cpp
    src/utility.cpp
    src/utility/ArgumentHandler.cpp
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS include
    FILES
      include/ranger/ArgumentHandler.h
      include/ranger/Data.h
      include/ranger/DataChar.h
      include/ranger/DataDouble.h
      include/ranger/DataFloat.h
      include/ranger/Forest.h
      include/ranger/ForestClassification.h
      include/ranger/ForestProbability.h
      include/ranger/ForestRegression.h
      include/ranger/ForestSurvival.h
      include/ranger/Tree.h
      include/ranger/TreeClassification.h
      include/ranger/TreeProbability.h
      include/ranger/TreeRegression.h
      include/ranger/TreeSurvival.h
      include/ranger/getopt.h
      include/ranger/globals.h
      include/ranger/utility.h
      include/ranger/version.h
)

target_include_directories(ranger
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ranger
)

target_link_libraries(ranger PUBLIC Threads::Threads)

add_executable(ranger-cli src/cli/main.cpp)
set_target_properties(ranger-cli PROPERTIES OUTPUT_NAME ranger)
target_include_directories(ranger-cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/ranger)
target_link_libraries(ranger-cli PRIVATE Ranger::ranger)

install(TARGETS ranger ranger-cli
  EXPORT RangerTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT RangerTargets
  FILE RangerTargets.cmake
  NAMESPACE Ranger::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Ranger
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/RangerConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RangerConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/RangerConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Ranger
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/RangerConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/RangerConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Ranger
)
