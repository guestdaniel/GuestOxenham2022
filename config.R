# Configure directories
root_directory = '/home/daniel/GuestOxenham2021'

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set plotting parameters
n_breaks = 10
breaks = c(seq(2, 10, 10/n_breaks) %o% 10^(-5:3))
labels = as.character(breaks)
labels[!(log10(breaks)%%1==0) & breaks != 2000] = ''
ticksizes = rep(.25, length(breaks))
ticksizes[log10(breaks)%%1==0] = 1
font_scale = 8
size_point=1.0
size_smooth=0.5

