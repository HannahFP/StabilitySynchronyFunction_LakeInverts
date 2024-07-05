# Load necessary libraries
library(ade4)
library(arm)
library(FD)
library(TMB)
library(dplyr)

# Set working directory
setwd("\\\\home.ansatt.ntnu.no/Hannafri/Documents/IS_Cleaned")

# Source utility functions
source("Scripts/MAR_UtilityFunctions.R")

# Compile and manage TMB model
compile("Scripts/tmb_MAR1_comm_lv_replicates.cpp")
file.copy(from = "Scripts/tmb_MAR1_comm_lv_replicates.dll", to = "tmb_MAR1_comm_lv_replicates.dll")

# Load manipulated data
processed_data <- readRDS("Data/lake_processed_data.RDS")
obs_N_list <- processed_data$obs_N_list
AQUA_Et1_list <- processed_data$AQUA_Et1_list
nr_replicates_list <- processed_data$nr_replicates_list
Abund_list <- processed_data$Abund_list

# Initialize variables
n.lakes <- length(obs_N_list)
MeanSYFSCVFD <- matrix(NA, n.lakes, 7)
DC_contributions <- numeric(n.lakes)
DC_CI <- matrix(NA, n.lakes, 2)
rho_e_contributions <- numeric(n.lakes)
rho_e_CI <- matrix(NA, n.lakes, 2)
rho_x_contributions <- numeric(n.lakes)
rho_x_CI <- matrix(NA, n.lakes, 2)
corlist <- list()

# Create empty lists to store results
my_list1 <- vector("list", n.lakes)
my_list2 <- vector("list", n.lakes)
my_list3 <- vector("list", n.lakes)
my_list4 <- vector("list", n.lakes)
my_list5 <- vector("list", n.lakes)
my_list6 <- vector("list", n.lakes)
my_list7 <- vector("list", n.lakes)
my_list8 <- vector("list", n.lakes)
my_list9 <- vector("list", n.lakes)
my_list10 <- vector("list", n.lakes)

pdf(file = "Plots/all_lakes.pdf")

# Define the function to run models and extract estimates for each lake
process_lake_model <- function(i) {
  tryCatch({
    # Load the manipulated data for the current lake
    obs_N <- obs_N_list[[i]]
    AQUA_Et1 <- AQUA_Et1_list[[i]]
    nr_replicates <- nr_replicates_list[[i]]
    Abund <- Abund_list[[i]]
    LakeName <- Lake_list[[i]]
    
    if (!is.null(obs_N) && !is.null(AQUA_Et1) && !is.null(nr_replicates)) {
      # Convert nr_replicates back to a table
      nr_replicates <- as.table(setNames(nr_replicates$Freq, nr_replicates$Year))
      
      # Calculate the total abundance per year
      AbundYr <- rowSums(Abund)
      
      # Calculate and store the stability as the reciprocal of the coefficient of variation
      CoV <- sd(AbundYr, na.rm = TRUE) / mean(AbundYr, na.rm = TRUE)
      stability <- 1 / CoV
      MeanSYFSCVFD[i,2] <<- stability
      
      ## Calculate correlation matrix to see how it compares to rho_x
      cor_mat <- cor(obs_N)
      cor_mat[upper.tri(cor_mat)] <- NA
      diag(cor_mat) <- NA
      
      my_list10[[i]] <<- cor_mat
      
      # Proceed only if there are at least 10 taxa
      if (nrow(obs_N) >= 10) {
        my_list5[[i]] <<- rownames(obs_N)
        
        # Estimate model parameters
        f11 <- estMARcomm_lv_overdisp(Y = obs_N, nr_replicates = nr_replicates, nu = 1, n_lv = 3, family = "Poisson", silent = TRUE)
        my_list4[[i]] <<- f11$out$model$message
        
        # Extract estimates and calculate correlations
        sumf1 <- extractEstimates(f11)
        pars <- sumf1$AR1_pars
        lnk <- pars[, 1] / (1 - pars[, 2])
        cor <- cor(log(rowMeans(obs_N)), lnk)
        corlist[[i]] <<- cor
        
        DC_values <- sumf1$DC
        rho_e_values <- sumf1$rho_e
        rho_x_values <- sumf1$rho_x
        
        rho_e_contributions[i] <<- mean(rho_e_values, na.rm = TRUE)
        rho_e_CI[i,] <<- quantile(rho_e_values, c(0.025, 0.975))
        
        rho_x_contributions[i] <<- mean(rho_x_values, na.rm = TRUE)
        rho_x_CI[i,] <<- quantile(rho_x_values, c(0.025, 0.975))
        
        DC_contributions[i] <<- mean(1-DC_values, na.rm = TRUE)
        DC_CI[i,] <<- quantile(1-DC_values, c(0.025, 0.975))
        
        # Store variance of parameters
        my_list6[[i]] <<- var(pars[, 1])
        my_list7[[i]] <<- var(pars[, 2])
        
        # Calculate temporal synchrony matrix
        SYRhoe <- sumf1$rho_e_mat
        SYRhox <- sumf1$rho_x_mat
        
        my_list1[[i]] <<- SYRhoe
        my_list9[[i]] <<- SYRhox
        MeanSYFSCVFD[i,1] <<- paste(LakeName)
        MeanSYFSCVFD[i,3] <<- mean(SYRhoe, na.rm = TRUE)
        MeanSYFSCVFD[i,4] <<- mean(SYRhox, na.rm = TRUE)
        
        # Calculate functional similarity matrix
        FS <- as.matrix(AQUA_Et1)
        FS <- 1 - FS
        my_list2[[i]] <<- FS
        MeanSYFSCVFD[i,5] <<- mean(FS[upper.tri(FS)], na.rm = TRUE)
        
        # Calculate functional diversity and store results
        filtered_abundance_matrix <- Abund[rowSums(Abund) != 0, ]
        FD <- dbFD(AQUA_Et1, filtered_abundance_matrix)
        MeanSYFSCVFD[i,6] <<- mean(FD$FRic, na.rm = TRUE)
        MeanSYFSCVFD[i,7] <<- mean(FD$FDis)
        my_list8[[i]] <<- FD$FRic
        
        # Store pairwise synchrony and functional similarity data
        my_list3[[i]] <<- data.frame(Lake = i, Pair = 1:length(SYRhoe[lower.tri(SYRhoe, diag = FALSE)]), SYRhoe = SYRhoe[lower.tri(SYRhoe, diag = FALSE)], FS = FS[lower.tri(FS, diag = FALSE)])
        
        # Plot functional similarity vs temporal synchrony
        plot(FS[lower.tri(FS, diag = FALSE)], SYRhoe[lower.tri(SYRhoe, diag = FALSE)], 
             xlab = "Functional similarity", ylab = "Temporal synchrony", main = paste("Lake", i), pch = 19)
        abline(lm(SYRhoe[lower.tri(SYRhoe, diag = FALSE)] ~ FS[lower.tri(FS, diag = FALSE)]), lwd = 3)
      }
    }
  }, error = function(e) { 
    cat("ERROR:", conditionMessage(e), "\n") 
  })
}

# Apply the process_lake_model function to each lake
for (i in 1:n.lakes) {
  process_lake_model(i)
}

dev.off()

# Save the results to file
saveRDS(list(
  MeanSYFSCVFD = MeanSYFSCVFD,
  DC_contributions = DC_contributions,
  DC_CI = DC_CI,
  rho_e_contributions = rho_e_contributions,
  rho_e_CI = rho_e_CI,
  rho_x_contributions = rho_x_contributions,
  rho_x_CI = rho_x_CI,
  corlist = corlist,
  my_list1 = my_list1,
  my_list2 = my_list2,
  my_list3 = my_list3,
  my_list4 = my_list4,
  my_list5 = my_list5,
  my_list6 = my_list6,
  my_list7 = my_list7,
  my_list8 = my_list8,
  my_list9 = my_list9,
  my_list10 = my_list10
), file = "Data/results.RDS")
