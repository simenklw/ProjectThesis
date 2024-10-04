library(arrow)
library(adagenet)
library(heirfast)


data <- read_feather("data/processed/massBV.feather")

ind <- as.character(data$ringnr)
population <- as.character(data$hatchisland)

locus <- data[,10:100]

Mydata <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
Mydata








