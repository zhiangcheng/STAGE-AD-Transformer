args = commandArgs(trailingOnly=TRUE)
infile1 = args[1]
infile2 = args[2]
output = args[3]

#########################################################################

### Packages ###
library(data.table)   # To load data
library(tidyverse)

### Load Data ###
sumstats1 <- read.table(infile1, header = T, stringsAsFactors = F)
sumstats2 <- read.table(infile2, header = T, stringsAsFactors = F)

### Correlation ###

cor <- cor(sumstats1$Effect, sumstats2$Effect)

# Print results to file

output_correlation <- paste0(output, "_correlation.txt")

cat("Correlation between two traits =", cor, file = output_correlation)

