#Import the data into the plotting notebook
library(tidyverse)
library(ggprism)
source('helperFunctions.R')




#uPAC DIA method test data. Missing the full data 
uPAC_singlecol <- read_SN_summary("data/SN17/R_ERWIN_stats.csv", thousand_sep = "\\.", col_sep = "\t") %>% 
                filter(!grepl("LIT-DIA", X)) %>% #removing all the LIT-DIA 
                mutate(Load = factor(str_match(X, pattern =  "mz_(.+?)_Pierce")[,2], levels = c('10ng', '5ng', '1ng', '250pg')),
                        resolution = factor(str_match(X, pattern =  "HRMS1_(.+?)_")[,2]),
                        window_size = factor(str_match(X, pattern = "win_(.+?)_")[,2]),
                        gradient = factor(str_match(X, pattern = "DIA_(.+?)_")[,2]),
                        scheme = "Single-column") %>% 
                        group_by(Load, gradient)
uPAC_singlecol$Peptides[c(2,10,19,33)] <- c(8660,17490, 24950,8080)
uPAC_singlecol$Protein.Groups[c(2,19)] <- c(2280,4680)
uPAC_singlecol$Load[c(17,18,19)] <- '10ng'


#uPAC with precol set-up

uPAC_precol <- read_SN_summary("data/SN17/R2_stats.csv", thousand_sep = "\\.", col_sep = "\t") %>% 
            mutate(Load = factor(str_match(X, pattern = "mz_(.+)_Pierce")[,2], , levels = c('10ng', '5ng', '1ng', '250pg')),
                   Gradient = str_match(X, pattern = "DIA_(.+)_200nl")[,2], scheme = "Pre-column") %>% group_by(Load, Gradient)
#Some manual edit as the terminal zero keeps getting dropped. Will figure this out later. Probably will just rewrite the summary python script so I do not have to deal with the dot comma separator 
uPAC_precol$Protein.Groups[c(1,17)] <- c(4800, 4360)
uPAC_precol$Peptides[c(5)] <- c(13460)



#Shorter gradients for 1ng 
short_1ng <- read_SN_summary('data/SN17/R3_stats.csv', col_sep = "\t", thousand_sep = "\\.") %>%
            mutate(Gradient = str_match(X, pattern = "DIA_(.+?)_200")[,2], Load = "1ng") %>% 
            group_by(Gradient)

short_1ng$Protein.Groups[6] <- c(2690)

#Shorter gradients for 250pg
parallel_1ng <- read_SN_summary('data/SN17/R4_stats.csv', col_sep ="\t", thousand_sep = "\\.") %>% 
                mutate(Gradient = str_match(X, pattern = "precol_(.+?)_200")[,2], Load = "250g") 


#Single-cell 96well, 120k vs 240k 

#Summaries for initial single-cell runs 
sc_96 <- read_SN_summary('data/SN17/R8_stats.csv', col_sep ="\t", thousand_sep = "\\.") %>% 
             mutate(resolution = ifelse(grepl("240k", X), "240k", "120k"),
                      injection = factor(ifelse(grepl("direct", X), "Direct", "Transfered"), levels = c("Transfered", "Direct")),
                      gradient = str_match(X, pattern = "precol_(.+?)_")[,2]) %>% group_by(resolution, gradient) 

sc_96$Peptides[c(5,17, 24)] <- c(1510,1710,2680)
sc_96$Protein.Groups[c(24)] <- c(1020)

sc_384 <- read_SN_summary('data/SN17/R10_stats.csv', col_sep ="\t", thousand_sep = "\\.") %>%
                mutate(resolution = ifelse(grepl("240k", X), "240k", "120k"),
                      injection = factor(ifelse(grepl("direct", X), "Direct", "Transfered"), levels = c("Transfered", "Direct")),
                      gradient = str_match(X, pattern = "precol_(.+?)_")[,2]) 
sc_384$Peptides[c(6,9)] <- c(3940,3960)
sc_test <- bind_rows(sc_96,sc_384) %>% group_by(resolution, gradient, injection)
                
#calculate summary statistics for plots 
uPAC_singlecol_summary <- calculate_summary(uPAC_singlecol)
uPAC_precol_summary <- calculate_summary(uPAC_precol)
short_1ng_summary <- calculate_summary(short_1ng)
parallel_1ng_summary <- calculate_summary(parallel_1ng)
sc_96_summary <- calculate_summary(sc_96)
sc_test_summary <- calculate_summary(sc_test)
