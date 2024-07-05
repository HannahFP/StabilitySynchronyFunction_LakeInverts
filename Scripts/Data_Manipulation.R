# Load necessary libraries
library(ade4)
library(dplyr)

# Set working directory
setwd("\\\\home.ansatt.ntnu.no/Hannafri/Documents/IS_Cleaned")

# Load data
IVrepl <- read.csv("Data/Replicate_Abun_Final.csv", header = TRUE, sep = ';', check.names = FALSE, encoding = "UTF-8")
Traits <- read.csv2("Data/Traits_DR_Final.csv", header = TRUE, sep = ';')

# Initialize variables
n.lakes <- 109

# Create empty lists to store manipulated data for each lake
obs_N_list <- vector("list", n.lakes)
AQUA_Et1_list <- vector("list", n.lakes)
nr_replicates_list <- vector("list", n.lakes)
Abund_list <- vector("list", n.lakes)
Lake_list <- vector("list", n.lakes)

# Function to process each lake's data
process_lake_data <- function(i) {
  tryCatch({
    # Subset the dataframe for the current lake
    LakeX <- IVrepl[IVrepl$LakeID == i, 2:(ncol(IVrepl))]
    
    # Remove species never observed by filtering columns where the sum of observations is greater than 0
    LakeX <- LakeX[, c(1:3, (4:ncol(LakeX))[colSums(LakeX[, 4:ncol(LakeX)]) > 0])]
    
    # Aggregate data by year
    tr <- aggregate(LakeX[, -c(1:3)], by = list(Year = LakeX$Year), sum)
    
    # Calculate the number of years each species was observed
    num_year_obs <- apply(tr[, -1] > 0, 2, sum)
    
    # Select species observed in more than 10 years
    sel_Species <- num_year_obs[num_year_obs > 10]
    
    # Define the columns to keep
    morecol <- c("Lake", "Year", "Rep_num")
    cols <- c(morecol, names(sel_Species))
    
    # Subset the dataframe to include only the specified columns
    LakeX <- LakeX[, cols]
    
    # Get the unique years present in the dataframe
    unYear <- unique(LakeX$Year)
    
    # Identify the missing years within the range of the minimum and maximum years
    missing_years <- (min(unYear):max(unYear))[!(min(unYear):max(unYear) %in% unYear)]

    # Check if there are missing years
    if (length(missing_years) > 0) {
      # Create a data frame for the missing years with the appropriate structure and data types
      tmp <- data.frame(
        Lake = rep(unique(LakeX$Lake), length(missing_years)),
        Year = missing_years,
        Rep_num = NA,
        matrix(NA, nrow = length(missing_years), ncol = ncol(LakeX) - 3)
      )
      
      # Ensure the column names are the same
      colnames(tmp) <- colnames(LakeX)
      
      # Bind the original dataframe with the new rows
      LakeX <- rbind(LakeX, tmp)
    } else {
      # Create an empty dataframe with the same structure as LakeX
      tmp <- data.frame(matrix(ncol = ncol(LakeX), nrow = 0))
      colnames(tmp) <- colnames(LakeX)
      LakeX <- rbind(LakeX, tmp)
    }
    
    # Summarize the abundance data by year
    Abund <- LakeX[, c(2, 4:ncol(LakeX))] %>% 
      group_by(Year) %>% 
      summarise(across(everything(), ~sum(.x, na.rm = TRUE)))
    
    # Select columns that match the taxa in the Traits dataframe
    col.num <- which(colnames(Abund) %in% Traits$Taxa)
    Abund <- Abund[, col.num]
    
    # Subset the data to include only observed taxa
    obs_N <- t(LakeX[, -(1:2)])
    obs_N <- subset(obs_N, rownames(obs_N) %in% Traits$Taxa)
    obs_N <- obs_N[order(row.names(obs_N)), ]
    
    # Ensure obs_N is numeric while maintaining row names
    rownames_obs_N <- rownames(obs_N)
    obs_N <- as.data.frame(obs_N, as.numeric)
    rownames(obs_N) <- rownames_obs_N
    
    obs_N <- as.matrix(obs_N)
    
    # Store the number of replicates per year
    nr_replicates <- as.data.frame(table(LakeX$Year))
    colnames(nr_replicates) <- c("Year", "Freq")
    
    # Proceed only if there are at least 10 taxa
    if (nrow(obs_N) >= 10) {
      obs_N_list[[i]] <<- obs_N
      nr_replicates_list[[i]] <<- nr_replicates
      Abund_list[[i]] <<- Abund
      Lake_list[[i]] <<- unique(LakeX$Lake)

      # Process trait data for functional similarity calculation
      Traitstop <- subset(Traits, Traits$Taxa %in% rownames(obs_N))
      rownames(Traitstop) <- Traitstop[, 1]
      Traitstop <- Traitstop[, -1]
      Traitstop <- Traitstop[order(row.names(Traitstop)), ]
      Traitstop <- prep.fuzzy.var(Traitstop, c(7, 2, 3, 4, 8, 4, 5, 8, 9, 8, 9, 3, 3, 5, 6))
      AQUA_Et1 <- dist.ktab(ktab.list.df(list(Traitstop)), type = c("F"), option = c("noscale"))
      
      AQUA_Et1_list[[i]] <<- AQUA_Et1
      
    }
  }, error = function(e) { 
    cat("ERROR:", conditionMessage(e), "\n") 
  })
}

# Apply the process_lake_data function to each lake
for (i in 1:n.lakes) {
  process_lake_data(i)
}

# Save the manipulated data to file
saveRDS(list(
  obs_N_list = obs_N_list,
  AQUA_Et1_list = AQUA_Et1_list,
  nr_replicates_list = nr_replicates_list,
  Abund_list = Abund_list,
  Lake_list = Lake_list
), file = "Data/lake_processed_data.RDS")
