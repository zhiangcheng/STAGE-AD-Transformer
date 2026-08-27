# Use jack knife to test whether two genetic correlations are different 


args = commandArgs(trailingOnly=TRUE)
infile_rg1 = args[1]
infile_rg2 = args[2]
infile_c1 = args[3]
infile_c2 = args[4]
infile_hsq1 = args[5]
infile_hsq2 = args[6]
infile_hsq3 = args[7]
infile_hsq4 = args[8]
trait = args[9]
directory = args[10]


#########################################################################


### Packages ###

library(tidyverse)


### Read In data ###

# Global genetic correlation 1
rg1 <- read.table(file = infile_rg1, header = T)

# Delete file for 2nd genetic covariance
rg2 <- read.table(file = infile_rg2, header = T)



# Delete file for first genetic covariance
delete_c1 <- read.table(file = infile_c1, header = F, col.names = "gen_cov1")

# Delete file for second genetic covariance
delete_c2 <- read.table(file = infile_c2, header = F, col.names = "gen_cov2")



# Delete file heritability of trait 1 from correlation 1
delete_hsq1 <- read.table(file = infile_hsq1, header = F, col.names = "hsq1")

# Delete file heritability of trait 2 from correlation 1
delete_hsq2 <- read.table(file = infile_hsq2, header = F, col.names = "hsq2")

# Delete file heritability of trait 1 from correlation 2
delete_hsq3 <- read.table(file = infile_hsq3, header = F, col.names = "hsq3")

# Delete file heritability of trait 2 from correlation 2
delete_hsq4 <- read.table(file = infile_hsq4, header = F, col.names = "hsq4")



### Global rg values ###

global_rg1 <- rg1 %>% 
  filter(p2 == trait) %>% 
  pull(rg)

global_rg2 <- rg2 %>% 
  filter(p2 == trait) %>% 
  pull(rg)


### Delete values (estimated by systematically excluding 200 blocks) ###

delete1 <- delete_c1 %>% 
  bind_cols(delete_hsq1) %>% 
  bind_cols(delete_hsq2) %>% 
  mutate(rg1 = gen_cov1 / sqrt(hsq1 * hsq2))

delete2 <- delete_c2 %>% 
  bind_cols(delete_hsq3) %>% 
  bind_cols(delete_hsq4) %>% 
  mutate(rg2 = gen_cov2 / sqrt(hsq3 * hsq4))

delete <- delete1 %>% 
  bind_cols(delete2) %>% 
  mutate(rg_diff = rg1 - rg2)


### Jackknife pseudovalues ###

n_blocks <- nrow(delete)

jackknife <- delete %>% 
  mutate(pseudovalues = n_blocks * (global_rg1 - global_rg2) - (n_blocks - 1) * rg_diff)


### Z-score and p-value for rg1 and rg2 being different ###

jackknife_estimate <- jackknife %>% 
  summarise(mean_pseudovalues = mean(pseudovalues)) %>% 
  pull(mean_pseudovalues)

jackknife_variance <- jackknife %>% 
  summarise(var_pseudovalues = var(pseudovalues)) %>% 
  pull(var_pseudovalues)

jackknife_se <- sqrt(jackknife_variance/n_blocks)

z_score <- jackknife_estimate / jackknife_se
p_value <- 2 * (1 - pnorm(abs(z_score)))

result <- tibble(
  trait = trait,
  Z = z_score,
  P = p_value)

write.table(result, 
            file = paste0(directory, "Comparison_of_female_and_male_rg_with_", trait, ".txt", sep = ""),
            quote = F, sep = "\t", col.names = T, row.names = F)
