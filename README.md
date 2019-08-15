[![Travis Build Status](https://travis-ci.org/imbs-hl/ranger.svg?branch=master)](https://travis-ci.org/imbs-hl/ranger)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/imbs-hl/ranger?branch=master&svg=true)](https://ci.appveyor.com/project/mnwright/ranger)
[![Coverage Status](https://coveralls.io/repos/github/imbs-hl/ranger/badge.svg?branch=master)](https://coveralls.io/github/imbs-hl/ranger?branch=master)
![CRAN Downloads month](http://cranlogs.r-pkg.org/badges/ranger?color=brightgreen)
![CRAN Downloads overall](http://cranlogs.r-pkg.org/badges/grand-total/ranger?color=brightgreen)
## ranger: A Fast Implementation of Random Forests
Marvin N. Wright

### Introduction
ranger is a fast implementation of random forests (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, and survival forests are supported. Classification and regression forests are implemented as in the original Random Forest (Breiman 2001), survival forests as in Random Survival Forests (Ishwaran et al. 2008). Includes implementations of extremely randomized trees (Geurts et al. 2006) and quantile regression forests (Meinshausen 2006).

ranger is written in C++, but a version for R is available, too. We recommend to use the R version. It is easy to install and use and the results are readily available for further analysis. The R version is as fast as the standalone C++ version.

### Installation
#### R version
To install the ranger R package from CRAN, just run

```R
install.packages("ranger")
```

R version >= 3.1 is required. With recent R versions, multithreading on Windows platforms should just work. If you compile yourself, the new RTools toolchain is required.

To install the development version from GitHub using `devtools`, run

```R
devtools::install_github("imbs-hl/ranger")
```

#### Standalone C++ version
To install the C++ version of ranger in Linux or Mac OS X you will need a compiler supporting C++11 (i.e. gcc >= 4.7 or Clang >= 3.0) and Cmake. To build start a terminal from the ranger main directory and run the following commands

```bash
cd cpp_version
mkdir build
cd build
cmake ..
make
```

After compilation there should be an executable called "ranger" in the build directory. 

To run the C++ version in Microsoft Windows please cross compile or ask for a binary.

### Usage
#### R version
For usage of the R version see ?ranger in R. Most importantly, see the Examples section. As a first example you could try 

```R  
ranger(Species ~ ., data = iris)
```

#### Standalone C++ version
In the C++ version type 

```bash
./ranger --help 
```

for a list of commands. First you need a training dataset in a file. This file should contain one header line with variable names and one line with variable values per sample (numeric only). Variable names must not contain any whitespace, comma or semicolon. Values can be seperated by whitespace, comma or semicolon but can not be mixed in one file. A typical call of ranger would be for example

```bash
./ranger --verbose --file data.dat --depvarname Species --treetype 1 --ntree 1000 --nthreads 4
```

If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. 

Please cite our paper if you use ranger.

### References
* Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. https://doi.org/10.18637/jss.v077.i01.
* Schmid, M., Wright, M. N. & Ziegler, A. (2016). On the use of Harrell's C for clinical risk prediction via random survival forests. Expert Syst Appl 63:450-459. https://doi.org/10.1016/j.eswa.2016.07.018.
* Wright, M. N., Dankowski, T. & Ziegler, A. (2017). Unbiased split variable selection for random survival forests using maximally selected rank statistics. Stat Med 36:1272-1284. https://doi.org/10.1002/sim.7212.
* Nembrini, S., König, I. R. & Wright, M. N. (2018). The revival of the Gini Importance? Bioinformatics. https://doi.org/10.1093/bioinformatics/bty373.
* Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. https://doi.org/10.1023/A:1010933404324.
* Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. Ann Appl Stat 2:841-860. https://doi.org/10.1097/JTO.0b013e318233d835.
* Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods Inf Med 51:74-81. https://doi.org/10.3414/ME00-01-0052.
* Hastie, T., Tibshirani, R., Friedman, J. (2009). The Elements of Statistical Learning. Springer, New York. 2nd edition.
* Geurts, P., Ernst, D., Wehenkel, L. (2006). Extremely randomized trees. Mach Learn 63:3-42. https://doi.org/10.1007/s10994-006-6226-1.
* Meinshausen (2006). Quantile Regression Forests. J Mach Learn Res 7:983-999. http://www.jmlr.org/papers/v7/meinshausen06a.html.
* Sandri, M. & Zuccolotto, P. (2008). A bias correction algorithm for the Gini variable importance measure in classification trees. J Comput Graph Stat, 17:611-628. https://doi.org/10.1198/106186008X344522.
* Coppersmith D., Hong S. J., Hosking J. R. (1999). Partitioning nominal attributes in decision trees. Data Min Knowl Discov 3:197-217. https://doi.org/10.1023/A:1009869804967.
