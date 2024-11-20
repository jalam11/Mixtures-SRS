
# load packages

## Reading and saving files
library(readr); library(writexl) # for read_csv(); write_xlsx()  
library(here)

# Cleaning
library(plyr) # mapvalues()
library(MASS)
library(dplyr)
library(tidyr) # for gather()
select <- dplyr::select

## Visualizing
library(ggplot2)
library(ggrepel); library(ggpubr); library(scales); library(RColorBrewer) # extra ggplot features
library(signs) # for signs_format(), uses proper minus sign instead of hyphen
library(RColorBrewer) # For selecting the colour palettefor the figures

library(flextable)
library(gtsummary)
library(gt)
library(labelled)
library(officer)

# Descriptive statistics
library(EnvStats) # For calculating geoMean, geoSD
library(reshape2) # For melt(), used in correlation matrix

# Analysis
library(gWQS)
library(qgcomp)
library(qgcompint)

library(mice)
library(mitools)






