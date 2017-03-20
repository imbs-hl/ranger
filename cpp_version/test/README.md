To test the utility functions in C++, first install the google test framework (in this directory):

    wget https://github.com/google/googletest/archive/release-1.7.0.tar.gz
    tar xf release-1.7.0.tar.gz
    mv googletest-release-1.7.0 gtest-1.7.0

Then run the usual cmake commands:

    mkdir build
    cd build
    cmake ..
    make
    ./runUnitTests
