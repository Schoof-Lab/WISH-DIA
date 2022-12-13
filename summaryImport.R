#Import the data into the plotting notebook

library(tidyverse)
library(ggprism)
source('helperFunctions.R')

# DIA windows size for different input load

window_size <- read_SN_summary("data/window_stats.csv") %>% 
            mutate(Load = factor(str_match(X, pattern = "45_(.+)_.\\.raw")[,2], levels = rev(c("1ng", "5ng", "10ng", "100ng"))),
                   window_size = factor(str_match(X, pattern = "w_(.+?)_")[,2]), 
                   window_number = factor(str_match(X, pattern = "DIA_(.+?)_")[,2]),
                   resolution = factor(str_match(X, pattern = "mz_(.+?)_")[,2]),
                   injectionTime = factor(str_match(X, pattern = "_(.{1,6})_1CV")[,2]), # need to limit the amount of characters as lazy matching does not seem to work from end 
                   MS2_type = factor(str_match(X, pattern ="40SPD_(.+?)-")[,2])) %>% 
                   group_by(Load,window_size,window_number,resolution,injectionTime, MS2_type) %>% # some further modification of variables for nicer plots 
                   mutate(window_size = gsub("var-mz", "100", window_size),
                          window_size = factor(gsub("mz", "", window_size), levels = c("10", "20", "40", "80", "100")))
                         


# HRMS1-aquisition method for 1ng input 

hrms1 <- read_SN_summary('data/DIA-R6_stats.csv', thousand_sep = "\\.") %>% 
           mutate(Load = "1ng",
                  resolution = factor(str_match(X, pattern = '120k-(.+?)-')[,2]),
                  window_size = gsub("mz", "", factor(str_match(X, pattern = '-(.mz?|..mz?)-')[,2])),
                  method = case_when(window_size %in% c("8", "5") ~ 'HRMS1-DIA', window_size %in% c("10", "15") ~ "DIA"),
                  MS2_type = case_when(resolution == "Normal" ~ "LIT", resolution != "Normal" ~ 'OT')) %>% 
                group_by(resolution, method, MS2_type)

#Some manual editing due to lost trailing zeros by gsub...
hrms1$Precursors[8] <- 5130;hrms1$Peptides[7] <- 4640

#HRMS1 vs DIA 1ng input with different window sizes
hrms1_win <- read_SN_summary('data/R7_stats.csv', thousand_sep = "\\.") %>%
           mutate(Load = "1ng",
                  ms1_resolution = factor(c(str_match(X, pattern = 'DIA_(.+?)_')[1:12,2],str_match(X, pattern = 'HRMS1_(.+?)_')[13:24,2])),
                  ms2_resolution = factor(str_match(X, pattern = 'k_(.+?)_')[,2], levels = c('15k','30k', '60k', '120k', '240k')),
                  window_size = factor(gsub("mz", "", factor(str_match(X, pattern = 'w_(.+?)_')[,2])), levels = c('10', '15', '20', '30', '40', '60', '100', '120')),
                  method = str_match(X, pattern = '40SPD_(.+?)_')[,2]) %>% 
                  group_by(method, window_size, ms2_resolution) 

hrms1_win$Precursors[c(7,8,12,16,18)] <- c(3130,3130, 5550, 4740, 4750)
hrms1_win$Peptides[c(9)] <- c(3040)
hrms1_win$Protein.Groups[c(10,18)] <- c(1440,1290)


#calculate summary statistics for plots 
window_size_summary <- calculate_summary(window_size)
hrms1_summary <- calculate_summary(hrms1)
hrms1_win_summary <- calculate_summary(hrms1_win)