# Rangerts, a ranger package for time series

![Random forest](www/random-forest.jpg?raw=true "For fun")
credit: https://thinkr.fr/premiers-pas-en-machine-learning-avec-r-volume-4-random-forest/

## The original package
*ranger* (https://github.com/imbs-hl/ranger) is an open source R package on github, created and maintained by Marvin N. Wright, with a clear explanation in the article below :    
* Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Software 77:1-17. https://doi.org/10.18637/jss.v077.i01.

## What does rangerts do
In the *rangerts* package, we try to get out of the comfort zone of random forest applications. We try to predict a time series, with exogenous variables.     

The main idea of this modified version is to test the random forest algorithm using the block bootstrapping (help take time dependency of the data into account) during the period of tree growing, instead of standard resampling mode. To see whether this could help improve model's accuracy.       

In order to benefit from the efficient implementation of the *ranger* package, we based on its c++ codes and we added 4 different kinds of block bootstrapping: non-overlapping blocks, moving blocks, stationary blocks, and circular blocks.        

## New parameters
All functions are the same as the *ranger* package (even the main function name :stuck_out_tongue:).    
We added 3 parameters in the ranger function: bootstrap.ts, by.end, block.size. All these parameters are to be used **with caution** together with the parameters already included in the *ranger* original function.   

* bootstrap.ts: string parameter, empty (= NULL) for iid bootstrap, or takes its value in **nonoverlapping**, **moving**, **stationary**, **circular**, and **seasonal**, by default = NULL. Some research works have demonstrated that moving or stationary might be more beneficial.    
* by.end: boolean, by default = TRUE, build block from the end to the start of time series.     
* block.size: the number of observations per block, by default = 10. In the **stationary** block bootstrapping mode, this parameter define the geometric law with $p = 1/block.size$.     
* period: the number of observations per period, if **seasonal** bootstrap is selected.

## Installation
To install the development version from GitHub using `devtools`, run
```R
# quiet = TRUE to mask c++ compilation messages, optional
devtools::install_github("hyanworkspace/rangerts", quiet = T)
```


## Key parameters
### Bootstrap mode selection and sample fraction
The default bootstrapping method is the **i.i.d mode, with replacement** in the *ranger* package.    

Variantes exist when changing the parameter `replace = FALSE` or give a weight vector `case.weights` over the training observations to modify the probabilities that some observations will have more chance to be selected in the bag.      

Fraction of observations to sample is 1 by default for sampling with replacement and 0.632 ( ` = (exp(1)-1)/exp(1)` ) for sampling without replacement. This could be changed manuelly by the parameter `sample.fraction` in the `ranger` function.     

Among all the block bootstrapping we implemented, by the nature of their design,  except *nonoverlapping* and *seasonal* (under condition), all others have replacement in the bootstrapped sample inbag.    

The parameter `replace = TRUE` or `FALSE` makes no difference if you use these three block bootstrapping modes, `bootstrap.ts = "moving"`, `bootstrap.ts = "stationary"`, `bootstrap.ts = "circular"`, and `replace = FALSE` with the `bootstrap.ts = "seasonal"` can be dangerous if block.size and period are not properly chosen.     

### Block size
Our experiments have shown that this is the key parameters to be tuned. As we are treating with time series and weekly data, candidate values can be 4-5 (almost a month), or 52 (a year). Small values might be beneficial to be tested too.    

We suggest train a standard ranger model then study the autocorrelation of the residuals, to get some hint on what values to take for the block size, or maybe directly study the autocorrelation of the target variable.    

## Examples
We provide an open source dataset of French weekly electricity consumption, along with several features :     

* Time : observation index
* Day : day of month
* Month : month of year
* Year : year
* NumWeek : This feature goes from 0 to 1, from 1st January to 31th December, and increases linearly
* Load : electricity consumption in MW, target variable to predict
* Load1 : lag 1 of load
* Temp : temperature
* Temp1 : lag 1 of temperature
* IPI : industrial index

``` R
library(rangerts)
# to check the function ranger function helper
?rangerts::rangerts

# load consumption data in the package ----
data <- rangerts::elec_data

# feature engineering
data$Time2 <- data$Time^2
data$TempTime <- data$Time*data$Temp
data$TempChauf <- pmin(0, data$Temp - 15)
data$TempChaufTime <- pmin(0, data$Temp - 15) * data$Time

noel <- which(abs(data$Day - 24) <= 3 & data$Month == 12)
consoNoel = vector("numeric", length(data$Time))
consoNoel[noel] = 1
data$consoNoel <- consoNoel
data$MonthF <- as.factor(data$Month)

# split train and test
df_train <- data %>%
  dplyr::filter(Test == 0) %>%
  dplyr::select(- Test)

df_test <- data %>%
  dplyr::filter(Test == 1) %>%
  dplyr::select(- Test)

# set general parameters
nb_trees <- 1000
mtry <- floor(sqrt(ncol(df_train)))
block_size <- 52

# Use case 1
# the default ranger with bootstrap i.i.d and with replacement
# thus the sample fraction is the default value = 1
rf_iid_rep <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T,
                 seed = 1, # for reproductibility
                 keep.inbag = T) # to keep trace of in-bag samples
# 679 observations in total
nrow(df_train)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_iid_rep$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 429.52
# total number of inbag samples
purrr::map_dbl(rf_iid_rep$inbag.counts, sum) %>%
  mean()
# 679

# Use case 2
# the default ranger with bootstrap i.i.d and with replacement
# thus the sample fraction = 0.632
rf_iid <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = F,
                 seed = 1,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_iid$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 429
# total number of inbag samples
purrr::map_dbl(rf_iid$inbag.counts, sum) %>%
  mean()
# 429

# Use case 3
# the nonoverlapping mode with replacement
# thus the sample fraction is the default value = 1
rf_no_rep <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T, # default = T too
                 seed = 1,
                 bootstrap.ts = "nonoverlapping",
                 block.size = block_size,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_no_rep$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 439.266
# total number of inbag samples
purrr::map_dbl(rf_no_rep$inbag.counts, sum) %>%
  mean()
# 679

# Use case 4
# the nonoverlapping mode with replacement
# thus the sample fraction is the default value = 1
rf_no <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = F, # in this case, every sample in-bag is taken only once
                 seed = 1,
                 bootstrap.ts = "nonoverlapping",
                 block.size = block_size,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_no$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 429
# total number of inbag samples
purrr::map_dbl(rf_no$inbag.counts, sum) %>%
  mean()
# 429


# Use case 5
# the moving mode with replacement
# thus the sample fraction is the default value = 1
rf_mv <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T, # default = T too
                 seed = 1,
                 bootstrap.ts = "moving",
                 block.size = block_size,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_mv$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 430.157
# total number of inbag samples
purrr::map_dbl(rf_mv$inbag.counts, sum) %>%
  mean()
# 679

# Use case 6
# the stationary mode with replacement
# thus the sample fraction is the default value = 1
rf_st <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T, # default = T too
                 seed = 1,
                 bootstrap.ts = "stationary",
                 block.size = block_size,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_st$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 446.788

# Use case 7
# the circular mode with replacement
# thus the sample fraction is the default value = 1
rf_cr <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T, # default = T too
                 seed = 1,
                 bootstrap.ts = "circular",
                 block.size = block_size,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_cr$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 438.124

# Use case 8
# the seasonal mode with replacement
# thus the sample fraction is the default value = 1
rf_se <- rangerts::rangerts(Load ~ ., data = df_train,
                 num.trees = nb_trees,
                 mtry = mtry,
                 replace = T, # default = T too
                 seed = 1,
                 bootstrap.ts = "seasonal",
                 block.size = block_size,
                 period = 4,
                 keep.inbag = T)
# the average number of different observations
# that are at least taken in-bag once in the trees
purrr::map_int(rf_se$inbag.counts,
               ~ length(which(.x != 0))) %>%
  mean()
# 438.124

# final model list
model_list <- list(rf_iid, rf_iid_rep,
                   rf_no, rf_no_rep,
                   rf_mv, rf_st,
                   rf_cr, rf_se)

# compare rmse & mape
algo_spec <- c("iid without replacement",
               "iid with replacement",
               "nonoverlapping without replacement",
               "nonoverlapping with replacement",
               "moving with replacement",
               "stationary with replacement",
               "circular with replacement",
               "seasonal with replacement"
               )
rmse <- purrr::map_dbl(model_list,
               ~ yardstick::rmse_vec(df_test$Load,
                                   predict(.x, df_test)$predictions))
cbind(algo_spec, round(rmse, 2))

mape <- purrr::map_dbl(model_list,
               ~ yardstick::mape_vec(df_test$Load,
                                   predict(.x, df_test)$predictions))
cbind(algo_spec, round(mape, 2))
```
