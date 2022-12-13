# Function for Spectronaut summary data transformation 

#Importing data 

read_SN_summary <- function(file_path, thousand_sep = "\\,", col_sep = ","){
       # file_path - the path to the file to read
       temp <- read.csv(file_path, sep = col_sep) %>%
            mutate(Precursors = as.numeric(gsub(thousand_sep, "", Precursors)),
                   Peptides = as.numeric(gsub(thousand_sep, "", Peptides)), 
                   Protein.Groups = as.numeric(gsub(thousand_sep, "", Protein.Groups))) 
       as_tibble(temp) # returns the read data 
}


#Calculating summaries 

calculate_summary <- function(df){
       #Takes the output generate by read_SN_summary and calculate the average number of detected precursors, peptides and proteins in replicates.
       #df has to be pregrouped for the variables to calculate the right summaries e.g. Gradient used or sample Load
       df %>%  summarise(Precursors = round(mean(Precursors)),
                        Peptides = round(mean(Peptides)),
                        Protein.Groups = round(mean(Protein.Groups)))
}


# SILAC analysis functions 