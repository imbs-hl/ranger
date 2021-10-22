library(Mcomp)
library(M4comp2018)
library(tidyverse)
library(rangerts)
M3

niveau_acf <- 0.9

# # yearly ----
# ts <- subset(M3, "yearly")
# purrr::map_int(ts,  ~.x$n) %>% max()
# # 41
# 
# ts_M4 <- Filter(function(l) l$period == "Yearly", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 835
# 
# # others ----
# ts <- subset(M3, "other")
# purrr::map_int(ts,  ~.x$n) %>% max()
# # 96
# 
# ts_M4 <- Filter(function(l) l$period == "Weekly", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 2597
# 
# ts_M4 <- Filter(function(l) l$period == "Daily", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 9919
# 
# ts_M4 <- Filter(function(l) l$period == "Hourly", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 960
# 
# 
# # monthly data ----
# ts <- subset(M3, "monthly")
# purrr::map_int(ts,  ~.x$n) %>% max()
# # 126
# 
# 
# ts_M4 <- Filter(function(l) l$period == "Monthly", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 2794
# 
# 
# # quarterly ----
# ts <- subset(M3, "quarterly")
# purrr::map_int(ts,  ~.x$n) %>% max()
# # 64
# 
# ts_M4 <- Filter(function(l) l$period == "Quarterly", M4)
# purrr::map_int(ts_M4,  ~.x$n) %>% max()
# # 866


# m3 ----
ts_q <- subset(M3, "quarterly")
ts_m <- subset(M3, "monthly")

ts <- append(ts_q, ts_m)

res_m3 <- map_dfr(seq_along(ts), function(i) {
  
  this_ts <- ts[[i]]
  
  score <- NULL
  if (this_ts$period == "QUARTERLY") {
    df <- tibble(obs = c(this_ts$x, this_ts$xx)) %>% 
      mutate(quarter = rep(1:4, length.out = n()),
             time = 1:n())
    
    df_train <- df %>% slice(1:this_ts$n)
    df_test <- df %>% slice((this_ts$n+1):nrow(df))
    
    rf_iid <- rangerts(obs ~ ., df_train,
                       seed = i)
    prev_iid <- predict(rf_iid, df_test)$predictions
    
    acf_res <- acf(df_train$obs, plot = F)
    block <- last(which(acf_res$acf >= niveau_acf)) - 1
    
    rf_nono <- rangerts(obs ~ ., df_train, 
                      bootstrap.ts = "nonoverlapping", 
                      by.end = FALSE,
                      block.size = max(2, block),
                      seed = i)
    prev_nono <- predict(rf_nono, df_test)$predictions
    
    rf_mov <- rangerts(obs ~ ., df_train, 
                     bootstrap.ts = "moving", 
                     by.end = FALSE,
                     block.size = max(2, block),
                     seed = i)
    prev_mov <- predict(rf_mov, df_test)$predictions
    
    score <- tibble(iid_rmse = yardstick::rmse_vec(df_test$obs, prev_iid),
                    iid_mape = yardstick::mape_vec(df_test$obs, prev_iid),
                    nono_rmse = yardstick::rmse_vec(df_test$obs, prev_nono),
                    nono_mape = yardstick::mape_vec(df_test$obs, prev_nono),
                    mov_rmse = yardstick::rmse_vec(df_test$obs, prev_mov),
                    mov_mape = yardstick::mape_vec(df_test$obs, prev_mov)) %>% 
      mutate(freq = "quarterly")
    
    cat("MAPE: Nono:", score$nono_mape, "Mov:", score$mov_mape, "\n")
    
    
  } else if (this_ts$period == "MONTHLY") {
    
    df <- tibble(obs = c(this_ts$x, this_ts$xx)) %>% 
      mutate(month = rep(1:12, length.out = n()),
             time = 1:n())
    
    df_train <- df %>% slice(1:this_ts$n)
    df_test <- df %>% slice((this_ts$n+1):nrow(df))
    
    rf_iid <- rangerts(obs ~ ., df_train,
                       seed = i)
    prev_iid <- predict(rf_iid, df_test)$predictions
    
    acf_res <- acf(df_train$obs - predict(rf_iid, df_train)$predictions, 
                   plot = F)
    block <- last(which(acf_res$acf >= niveau_acf)) - 1
    
    rf_nono <- rangerts(obs ~ ., df_train, 
                      bootstrap.ts = "nonoverlapping", 
                      by.end = FALSE,
                      block.size = max(2, block),
                      seed = i)
    prev_nono <- predict(rf_nono, df_test)$predictions
    
    rf_mov <- rangerts(obs ~ ., df_train, 
                     bootstrap.ts = "moving", 
                     by.end = FALSE,
                     block.size = max(2, block),
                     seed = i)
    prev_mov <- predict(rf_mov, df_test)$predictions
    
    score <- tibble(iid_rmse = yardstick::rmse_vec(df_test$obs, prev_iid),
                    iid_mape = yardstick::mape_vec(df_test$obs, prev_iid),
                    nono_rmse = yardstick::rmse_vec(df_test$obs, prev_nono),
                    nono_mape = yardstick::mape_vec(df_test$obs, prev_nono),
                    mov_rmse = yardstick::rmse_vec(df_test$obs, prev_mov),
                    mov_mape = yardstick::mape_vec(df_test$obs, prev_mov)) %>% 
      mutate(freq = "monthly")
    
    cat("MAPE: Nono:", score$nono_mape, "Mov:", score$mov_mape, "\n")
    
  } 
  score
  
})

write_rds(res_m3, paste0("results/m3_", niveau_acf, "_bis.rds"))
summary(res_m3)
boxplot(res_m3$mov_mape)


# m4 ----
ts_q <- Filter(function(l) l$period == "Quarterly", M4)
ts_m <- Filter(function(l) l$period == "Monthly", M4)
ts_w <- Filter(function(l) l$period == "Weekly", M4)
ts_d <- Filter(function(l) l$period == "Daily", M4)
ts_h <- Filter(function(l) l$period == "Hourly", M4)

ts <- append(ts_q, ts_m)
# %>% 
#   append(ts_w) %>% 
#   append(ts_d) %>% 
#   append(ts_h)

res_m4 <- map_dfr(seq_along(ts), function(i) {
  
  this_ts <- ts[[i]]
  
  if (this_ts$period == "Quarterly") {
    df <- tibble(obs = c(this_ts$x, this_ts$xx)) %>% 
      mutate(quarter = rep(1:4, length.out = n()),
             time = 1:n())
    
    df_train <- df %>% slice(1:this_ts$n)
    df_test <- df %>% slice((this_ts$n+1):nrow(df))
    
    rf_iid <- rangerts(obs ~ ., df_train, seed = i)
    prev_iid <- predict(rf_iid, df_test)$predictions
    
    acf_res <- acf(df_train$obs, plot = F)
    block <- last(which(acf_res$acf >= niveau_acf)) - 1
    
    rf_nono <- rangerts(obs ~ ., df_train, 
                      bootstrap.ts = "nonoverlapping", 
                      by.end = FALSE,
                      block.size = max(2, block),
                      seed = i)
    prev_nono <- predict(rf_nono, df_test)$predictions
    
    rf_mov <- rangerts(obs ~ ., df_train, 
                     bootstrap.ts = "moving", 
                     by.end = FALSE,
                     block.size = max(2, block),
                     seed = i)
    prev_mov <- predict(rf_mov, df_test)$predictions
    
    score <- tibble(iid_rmse = yardstick::rmse_vec(df_test$obs, prev_iid),
                    iid_mape = yardstick::mape_vec(df_test$obs, prev_iid),
                    nono_rmse = yardstick::rmse_vec(df_test$obs, prev_nono),
                    nono_mape = yardstick::mape_vec(df_test$obs, prev_nono),
                    mov_rmse = yardstick::rmse_vec(df_test$obs, prev_mov),
                    mov_mape = yardstick::mape_vec(df_test$obs, prev_mov)) %>% 
      mutate(freq = "quarterly")
    
    cat("MAPE: Nono:", score$nono_mape, "Mov:", score$mov_mape, "\n")
    
    
  } else if (this_ts$period == "Monthly") {
    
    df <- tibble(obs = c(this_ts$x, this_ts$xx)) %>% 
      mutate(month = rep(1:12, length.out = n()),
             time = 1:n())
    
    df_train <- df %>% slice(1:this_ts$n)
    df_test <- df %>% slice((this_ts$n+1):nrow(df))
    
    rf_iid <- rangerts(obs ~ ., df_train,
                       seed = i)
    prev_iid <- predict(rf_iid, df_test)$predictions
    
    acf_res <- acf(df_train$obs, plot = F)
    block <- last(which(acf_res$acf >= niveau_acf)) - 1
    
    rf_nono <- rangerts(obs ~ ., df_train, 
                      bootstrap.ts = "nonoverlapping", 
                      by.end = FALSE,
                      block.size = max(2, block),
                      seed = i)
    prev_nono <- predict(rf_nono, df_test)$predictions
    
    rf_mov <- rangerts(obs ~ ., df_train, 
                     bootstrap.ts = "moving", 
                     by.end = FALSE,
                     block.size = max(2, block),
                     seed = i)
    prev_mov <- predict(rf_mov, df_test)$predictions
    
    score <- tibble(iid_rmse = yardstick::rmse_vec(df_test$obs, prev_iid),
                    iid_mape = yardstick::mape_vec(df_test$obs, prev_iid),
                    nono_rmse = yardstick::rmse_vec(df_test$obs, prev_nono),
                    nono_mape = yardstick::mape_vec(df_test$obs, prev_nono),
                    mov_rmse = yardstick::rmse_vec(df_test$obs, prev_mov),
                    mov_mape = yardstick::mape_vec(df_test$obs, prev_mov)) %>% 
      mutate(freq = "monthly")
    
    cat("MAPE: Nono:", score$nono_mape, "Mov:", score$mov_mape, "\n")
    
  }
  score
  
})

summary(res_m4)
write_rds(res_m4, paste0("results/m4_", niveau_acf, ".rds"))


# beepr::beep()
