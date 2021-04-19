# result analysis
library(tidyverse)

# m3 res
res_files <- list.files("results/")
m3_res <- res_files[str_detect(res_files, "m3")]
# m3_res <- m3_res[!str_detect(m3_res, "bis")]
m3_res <- m3_res[!str_detect(m3_res, "block")]

res <- map_dfr(m3_res,  ~read_rds(paste0("results/", .x)) %>% 
      rename(nono_rmse = nono, mov_rmse = mov) %>% 
      mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))


res %>% 
  nest(data = !niveau_acf) %>% 
  mutate(summary = lapply(data, summary)) %>% 
  select(summary) %>% 
  .[[1]]

p <- ggplot(res, aes(x = freq, y = mov_mape, color = niveau_acf)) +
  geom_boxplot()

library(plotly)
ggplotly(p)


# m4
m4_res <- res_files[str_detect(res_files, "m4")]
m4_res <- m4_res[!str_detect(m4_res, "bis")]

res <- map_dfr(m4_res,  ~read_rds(paste0("results/", .x)) %>% 
                 rename(nono_rmse = nono, mov_rmse = mov) %>% 
                 mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))


res %>% 
  nest(data = !niveau_acf) %>% 
  mutate(summary = lapply(data, summary)) %>% 
  select(summary) %>% 
  .[[1]]

ggplot(res, aes(x = freq, y = nono_mape, color = niveau_acf)) +
  geom_boxplot()

ggplot(res, aes(x = freq, y = nono_rmse, color = niveau_acf)) +
  geom_boxplot()

ggplot(res, aes(x = freq, y = mov_mape, color = niveau_acf)) +
  geom_boxplot()

ggplot(res, aes(x = freq, y = mov_rmse, color = niveau_acf)) +
  geom_boxplot()

library(plotly)
ggplotly(p)



# m3 ----
res <- map_dfr(m3_res,  ~read_rds(paste0("results/", .x)) %>% 
                 rename(nono_rmse = nono, mov_rmse = mov) %>% 
                 mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))


# percentage
res %>% 
  group_by(niveau_acf) %>% 
  summarise(percentage = sum(as.numeric(mov_mape > 0)) / n())

# which ones are better ?
res_with_idx <- res %>%
  nest(data = !freq) %>% 
  mutate(data_t = lapply(data, function(x) {
    x %>% mutate(idx = rep(1:(nrow(x) / 5), times = 5))
  })) %>% 
  unnest(data_t) %>% 
  select(- data)
# 
# ts_q <- subset(M3, "quarterly")
# ts_m <- subset(M3, "monthly")
# 
# ts <- append(ts_q, ts_m)

# blocksize_m3 <- map_dfr(seq(0.5, 0.9, by = 0.1), function(niveau_acf) {
#   
#   map_dfr(seq_along(ts), function(id) {
#     
#     this_ts <- ts[[id]]
#     
#     df <- tibble(obs = c(this_ts$x, this_ts$xx))
#     
#     df_train <- df %>% slice(1:this_ts$n)
#     acf_res <- acf(df_train$obs, 
#                    plot = F)
#     block <- last(which(acf_res$acf >= niveau_acf)) - 1
#     
#     tibble(block_size = max(2, block), 
#            idx = id)
#     
#     
#   }) %>% 
#     mutate(niveau_acf = niveau_acf)
# })
#   
# write_rds(blocksize_m3, "results/blocksize_m3.rds")
blocksize_m3 <- read_rds("results/blocksize_m3.rds")

res_with_idx <- res_with_idx %>% 
  left_join(blocksize_m3 %>% 
              mutate(niveau_acf = as.factor(niveau_acf)), 
            by = c("niveau_acf", "idx"))

res_with_idx %>% 
  mutate(better_with_mov = as.numeric(mov_mape > 0),
         better_with_nono = as.numeric(nono_mape > 0)) %>% 
  group_by(niveau_acf, better_with_mov) %>% 
  summarise(mean_block_size = mean(block_size))

res_with_idx %>% 
  mutate(better_with_mov = as.numeric(mov_mape > 0),
         better_with_nono = as.numeric(nono_mape > 0)) %>% 
  group_by(niveau_acf, better_with_nono) %>% 
  summarise(mean_block_size = mean(block_size))



# m4 ----
# ts_q <- Filter(function(l) l$period == "Quarterly", M4)
# ts_m <- Filter(function(l) l$period == "Monthly", M4)
# 
# ts <- append(ts_q, ts_m)

# blocksize_m4 <- map_dfr(seq(0.5, 0.9, by = 0.1), function(niveau_acf) {
#   
#   map_dfr(seq_along(ts), function(id) {
#     
#     this_ts <- ts[[id]]
#     
#     df <- tibble(obs = c(this_ts$x, this_ts$xx))
#     
#     df_train <- df %>% slice(1:this_ts$n)
#     acf_res <- acf(df_train$obs, 
#                    plot = F)
#     block <- last(which(acf_res$acf >= niveau_acf)) - 1
#     
#     tibble(block_size = max(2, block), 
#            idx = id)
#     
#     
#   }) %>% 
#     mutate(niveau_acf = niveau_acf)
# })
# 
# write_rds(blocksize_m4, "results/blocksize_m4.rds")


m4_res <- res_files[str_detect(res_files, "m4")]
m4_res <- m4_res[!str_detect(m4_res, "bis")]
m4_res <- m4_res[!str_detect(m4_res, "blocksize")]

res <- map_dfr(m4_res,  ~read_rds(paste0("results/", .x)) %>% 
                 rename(nono_rmse = nono, mov_rmse = mov) %>% 
                 mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))

res_with_idx <- res %>%
  nest(data = !freq) %>% 
  mutate(data_t = lapply(data, function(x) {
    x %>% mutate(idx = rep(1:(nrow(x) / 5), times = 5))
  })) %>% 
  unnest(data_t) %>% 
  select(- data)

blocksize_m4 <- read_rds("results/blocksize_m4.rds")

res_with_idx <- res_with_idx %>% 
  left_join(blocksize_m4 %>% 
              mutate(niveau_acf = as.factor(niveau_acf)), 
            by = c("niveau_acf", "idx"))

res_with_idx %>% 
  mutate(better_with_mov = as.numeric(mov_mape > 0),
         better_with_nono = as.numeric(nono_mape > 0)) %>% 
  group_by(niveau_acf, better_with_mov) %>% 
  summarise(mean_block_size = mean(block_size))

res_with_idx %>% 
  mutate(better_with_mov = as.numeric(mov_mape > 0),
         better_with_nono = as.numeric(nono_mape > 0)) %>% 
  group_by(niveau_acf, better_with_nono) %>% 
  summarise(mean_block_size = mean(block_size))

res_with_idx %>% 
  mutate(better_with_mov = as.numeric(mov_mape > 0),
         better_with_nono = as.numeric(nono_mape > 0)) %>% 
ggplot(aes(x = freq, y = block_size, color = niveau_acf)) +
  geom_boxplot() + 
  facet_grid(~better_with_mov)


# test M4 ----

res_files <- list.files("results/")
m4_res <- res_files[str_detect(res_files, "m4")]
m4_res <- m4_res[!str_detect(m4_res, "bis")]
m4_res <- m4_res[!str_detect(m4_res, "stat")]
m4_res <- m4_res[!str_detect(m4_res, "blocksize")]
stat_m4 <- read_rds("results/stat_m4.rds")

res <- map_dfr(m4_res,  ~read_rds(paste0("results/", .x)) %>% 
                 mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))  %>% 
  bind_cols(bind_rows(stat_m4, stat_m4, stat_m4, stat_m4, stat_m4)) %>% 
  mutate(nono_nmape = (iid_mape - nono_mape) / iid_mape,
         mov_nmape = (iid_mape - mov_mape) / iid_mape,
         iid_nrmse_std = iid_rmse / std,
         iid_nrmse_mean = iid_rmse / mean,
         iid_nrmse_med = iid_rmse / median,
         iid_rmse_q = iid_rmse / q1_q3,
         nono_nrmse_std = nono_rmse / std,
         nono_nrmse_mean = nono_rmse / mean,
         nono_nrmse_med = nono_rmse / median,
         nono_rmse_q = nono_rmse / q1_q3,
         mov_nrmse_std = mov_rmse / std,
         mov_nrmse_mean = mov_rmse / mean,
         mov_nrmse_med = mov_rmse / median,
         mov_rmse_q = mov_rmse / q1_q3)

# Normalized RMSE ----
niveaux_acf <- seq(0.5, 0.9, by = 0.1)

# std with all test (freq and niveau_acf) ----
nrmse_std <- res %>% 
  group_by(freq, niveau_acf) %>%
  nest() %>% 
  mutate(metric = lapply(data, function(df) {
    df %>% 
      select(iid_nrmse_std, nono_nrmse_std, mov_nrmse_std) %>% 
      pivot_longer(cols = c(iid_nrmse_std, nono_nrmse_std, mov_nrmse_std)) %>% 
      mutate(name = str_remove(name, fixed("_nrmse_std"))) %>% 
      filter(!is.infinite(value)) %>% 
      mutate(value = as.numeric(value)) %>% 
      filter(value < 20)
  }))
  
# sapply(nrmse_std$metric, function(df) {
#   
#   ggplot(df, aes(x = value, color = name)) +
#       geom_density()
#   
# })

# p <- ggplot(nrmse_std, aes(x = value, color = name)) +
#   geom_density()
# ggplotly(p)


res_nrmse_std <- res %>% 
  group_by(freq, niveau_acf) %>%
  nest() %>% 
  mutate(by_group = lapply(data, function(df) {
    filter(df, !is.infinite(iid_nrmse_std))
  }))
  

walk(seq_along(res_nrmse_std$by_group), function(i) {
  
  cat("freq: ", res_nrmse_std$freq[i], "\n")
  cat("niveau_acf: ", niveaux_acf[res_nrmse_std$niveau_acf[i]], "\n")
  df <- res_nrmse_std$by_group[[i]]
  print(wilcox.test(x = df$iid_nrmse_std,
                    y = df$mov_nrmse_std, 
                    paired = T, alternative = "greater"))
})


for (acf in niveaux_acf) {
  
  for (frequen in unique(res$freq)) {
    
    cat("Test for results with niveau_acf =", acf, frequen, "data", "\n")
    res_nrmse_std <- res %>% 
      filter(!is.infinite(iid_nrmse_std)) %>% 
      filter(niveau_acf == acf, 
             freq == frequen)
    
    cat("------------------------- Moving -------------------------\n")
    print(wilcox.test(x = res_nrmse_std$iid_nrmse_std,
                      y = res_nrmse_std$mov_nrmse_std, 
                      paired = T, alternative = "greater"))
    
    cat("------------------------- Non-over -----------------------\n")
    print(wilcox.test(x = res_nrmse_std$iid_nrmse_std,
                      y = res_nrmse_std$nono_nrmse_std, 
                      paired = T, alternative = "greater"))
    
    cat("----------------------------------------------------------\n")
  }
  
}


# median
nrmse_med <- res %>% 
  select(iid_nrmse_med, nono_nrmse_med, mov_nrmse_med) %>% 
  pivot_longer(cols = c(iid_nrmse_med, nono_nrmse_med, mov_nrmse_med)) %>% 
  mutate(name = str_remove(name, fixed("_nrmse_med"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 2)

p <- ggplot(nrmse_med, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_med <- res %>% 
  filter(!is.infinite(iid_nrmse_med))

wilcox.test(x = res_nrmse_med$iid_nrmse_med,
            y = res_nrmse_med$mov_nrmse_med, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_med$iid_nrmse_med,
            y = res_nrmse_med$nono_nrmse_med, 
            paired = T, alternative = "greater")

# mean
nrmse_mean <- res %>% 
  select(iid_nrmse_mean, nono_nrmse_mean, mov_nrmse_mean) %>% 
  pivot_longer(cols = c(iid_nrmse_mean, nono_nrmse_mean, mov_nrmse_mean)) %>% 
  mutate(name = str_remove(name, fixed("_nrmse_mean"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 2)

p <- ggplot(nrmse_mean, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_mean <- res %>% 
  filter(!is.infinite(iid_nrmse_mean))

wilcox.test(x = res_nrmse_mean$iid_nrmse_mean,
            y = res_nrmse_mean$mov_nrmse_mean, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_mean$iid_nrmse_mean,
            y = res_nrmse_mean$nono_nrmse_mean, 
            paired = T, alternative = "greater")

# IQR
nrmse_q <- res %>% 
  select(iid_rmse_q, nono_rmse_q, mov_rmse_q) %>% 
  pivot_longer(cols = c(iid_rmse_q, nono_rmse_q, mov_rmse_q)) %>% 
  mutate(name = str_remove(name, fixed("_rmse_q"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 2)

p <- ggplot(nrmse_q, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_q <- res %>% 
  filter(!is.infinite(iid_rmse_q))

wilcox.test(x = res_nrmse_q$iid_rmse_q,
            y = res_nrmse_q$mov_rmse_q, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_q$iid_rmse_q,
            y = res_nrmse_q$nono_rmse_q, 
            paired = T, alternative = "greater")

# Normalized MAPE ----
nmape <- res %>% 
  select(nono_nmape, mov_nmape) %>% 
  pivot_longer(cols = c(nono_nmape, mov_nmape)) %>% 
  mutate(name = str_remove(name, fixed("_nmape"))) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(abs(value) < 1)

p <- ggplot(nmape, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nmape <- res %>% 
  filter(!is.infinite(mov_nmape))

wilcox.test(x = res_nmape$mov_nmape, alternative = "greater")
wilcox.test(x = res_nmape$nono_nmape, alternative = "greater")



# RMSE et MAPE ----
rmse <- res %>% 
  select(iid_rmse, nono_rmse, mov_rmse) %>% 
  pivot_longer(cols = c(iid_rmse, nono_rmse, mov_rmse)) %>% 
  mutate(name = str_remove(name, fixed("_rmse"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 1e4)


p <- ggplot(rmse, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_rmse <- res %>% 
  filter(!is.infinite(iid_rmse))

wilcox.test(x = res_rmse$iid_rmse,
            y = res_rmse$mov_rmse, 
            paired = T, alternative = "greater")

wilcox.test(x = res_rmse$iid_rmse,
            y = res_rmse$nono_rmse, 
            paired = T, alternative = "greater")


mape <- res %>% 
  select(iid_mape, nono_mape, mov_mape) %>% 
  pivot_longer(cols = c(iid_mape, nono_mape, mov_mape)) %>% 
  mutate(name = str_remove(name, fixed("_mape"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 1e3)


p <- ggplot(mape, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_mape <- res %>% 
  filter(!is.infinite(iid_mape))

wilcox.test(x = res_mape$iid_mape,
            y = res_mape$mov_mape, 
            paired = T, alternative = "greater")

wilcox.test(x = res_mape$iid_mape,
            y = res_mape$nono_mape, 
            paired = T, alternative = "greater")

# test M3 ----
res_files <- list.files("results/")
m3_res <- res_files[str_detect(res_files, "m3")]
m3_res <- m3_res[!str_detect(m3_res, "bis")]
m3_res <- m3_res[!str_detect(m3_res, "stat")]
m3_res <- m3_res[!str_detect(m3_res, "blocksize")]
stat_m3 <- read_rds("results/stat_m3.rds")

res <- map_dfr(m3_res,  ~read_rds(paste0("results/", .x)) %>% 
                 mutate(niveau_acf = str_remove(.x, ".rds"))) %>% 
  mutate(niveau_acf = 
           as.numeric(str_sub(niveau_acf, start = 4)),
         niveau_acf = as.factor(niveau_acf))  %>% 
  bind_cols(bind_rows(stat_m3, stat_m3, stat_m3, stat_m3, stat_m3)) %>% 
  mutate(nono_nmape = (iid_mape - nono_mape) / iid_mape,
         mov_nmape = (iid_mape - mov_mape) / iid_mape,
         iid_nrmse_std = iid_rmse / std,
         iid_nrmse_mean = iid_rmse / mean,
         iid_nrmse_med = iid_rmse / median,
         iid_rmse_q = iid_rmse / q1_q3,
         nono_nrmse_std = nono_rmse / std,
         nono_nrmse_mean = nono_rmse / mean,
         nono_nrmse_med = nono_rmse / median,
         nono_rmse_q = nono_rmse / q1_q3,
         mov_nrmse_std = mov_rmse / std,
         mov_nrmse_mean = mov_rmse / mean,
         mov_nrmse_med = mov_rmse / median,
         mov_rmse_q = mov_rmse / q1_q3)

# Normalized RMSE ----

# std
nrmse_std <- res %>% 
  select(iid_nrmse_std, nono_nrmse_std, mov_nrmse_std) %>% 
  pivot_longer(cols = c(iid_nrmse_std, nono_nrmse_std, mov_nrmse_std)) %>% 
  mutate(name = str_remove(name, fixed("_nrmse_std"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 20)

p <- ggplot(nrmse_std, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_std <- res %>% 
  filter(!is.infinite(iid_nrmse_std))

wilcox.test(x = res_nrmse_std$iid_nrmse_std,
            y = res_nrmse_std$mov_nrmse_std, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_std$iid_nrmse_std,
            y = res_nrmse_std$nono_nrmse_std, 
            paired = T, alternative = "greater")

# median
nrmse_med <- res %>% 
  select(iid_nrmse_med, nono_nrmse_med, mov_nrmse_med) %>% 
  pivot_longer(cols = c(iid_nrmse_med, nono_nrmse_med, mov_nrmse_med)) %>% 
  mutate(name = str_remove(name, fixed("_nrmse_med"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 2)

p <- ggplot(nrmse_med, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_med <- res %>% 
  filter(!is.infinite(iid_nrmse_med))

wilcox.test(x = res_nrmse_med$iid_nrmse_med,
            y = res_nrmse_med$mov_nrmse_med, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_med$iid_nrmse_med,
            y = res_nrmse_med$nono_nrmse_med, 
            paired = T, alternative = "greater")

# mean
nrmse_mean <- res %>% 
  select(iid_nrmse_mean, nono_nrmse_mean, mov_nrmse_mean) %>% 
  pivot_longer(cols = c(iid_nrmse_mean, nono_nrmse_mean, mov_nrmse_mean)) %>% 
  mutate(name = str_remove(name, fixed("_nrmse_mean"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 2)

p <- ggplot(nrmse_mean, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_mean <- res %>% 
  filter(!is.infinite(iid_nrmse_mean))

wilcox.test(x = res_nrmse_mean$iid_nrmse_mean,
            y = res_nrmse_mean$mov_nrmse_mean, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_mean$iid_nrmse_mean,
            y = res_nrmse_mean$nono_nrmse_mean, 
            paired = T, alternative = "greater")

# IQR
nrmse_q <- res %>% 
  select(iid_rmse_q, nono_rmse_q, mov_rmse_q) %>% 
  pivot_longer(cols = c(iid_rmse_q, nono_rmse_q, mov_rmse_q)) %>% 
  mutate(name = str_remove(name, fixed("_rmse_q"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 20)

p <- ggplot(nrmse_q, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nrmse_q <- res %>% 
  filter(!is.infinite(iid_rmse_q))

wilcox.test(x = res_nrmse_q$iid_rmse_q,
            y = res_nrmse_q$mov_rmse_q, 
            paired = T, alternative = "greater")

wilcox.test(x = res_nrmse_q$iid_rmse_q,
            y = res_nrmse_q$nono_rmse_q, 
            paired = T, alternative = "greater")

# Normalized MAPE ----
nmape <- res %>% 
  select(nono_nmape, mov_nmape) %>% 
  pivot_longer(cols = c(nono_nmape, mov_nmape)) %>% 
  mutate(name = str_remove(name, fixed("_nmape"))) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(abs(value) < 1)

p <- ggplot(nmape, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_nmape <- res %>% 
  filter(!is.infinite(mov_nmape))

wilcox.test(x = res_nmape$mov_nmape, alternative = "greater")
wilcox.test(x = res_nmape$nono_nmape, alternative = "greater")



# RMSE et MAPE ----
rmse <- res %>% 
  select(iid_rmse, nono_rmse, mov_rmse) %>% 
  pivot_longer(cols = c(iid_rmse, nono_rmse, mov_rmse)) %>% 
  mutate(name = str_remove(name, fixed("_rmse"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 1e4)


p <- ggplot(rmse, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_rmse <- res %>% 
  filter(!is.infinite(iid_rmse))

wilcox.test(x = res_rmse$iid_rmse,
            y = res_rmse$mov_rmse, 
            paired = T, alternative = "greater")

# t.test(x = res_rmse$iid_rmse,
#        y = res_rmse$mov_rmse, 
#        paired = TRUE, alternative = "greater")


wilcox.test(x = res_rmse$iid_rmse,
            y = res_rmse$nono_rmse, 
            paired = T, alternative = "greater")


mape <- res %>% 
  select(iid_mape, nono_mape, mov_mape) %>% 
  pivot_longer(cols = c(iid_mape, nono_mape, mov_mape)) %>% 
  mutate(name = str_remove(name, fixed("_mape"))) %>% 
  filter(!is.infinite(value)) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value < 1e3)


p <- ggplot(mape, aes(x = value, color = name)) +
  geom_density()
ggplotly(p)


res_mape <- res %>% 
  filter(!is.infinite(iid_mape))

wilcox.test(x = res_mape$iid_mape,
            y = res_mape$mov_mape, 
            paired = T)

wilcox.test(x = res_mape$iid_mape,
            y = res_mape$nono_mape, 
            paired = T)
t.test(x = res_mape$iid_mape,
       y = res_mape$nono_mape, 
       paired = TRUE, alternative = "greater")
