##Can only be run after big main analysis is completed and loaded

library(reshape2)
library(ggplot2)
library(dplyr)

setwd("\\\\home.ansatt.ntnu.no/Hannafri/Documents/IS_Cleaned")

colnames(MeanSYFSCVFD) <- c("Lake", "Stability", "SYRhoe", "SYRhox", "FS", "FRic", "FDis")
MeanSYFSCVFD2 <- as.data.frame(MeanSYFSCVFD[complete.cases(MeanSYFSCVFD),])
MeanSYFSCVFD2[ , c(2:7)] <- apply(MeanSYFSCVFD2[ , c(2:7)], 2,           
                                  function(x) as.numeric(x))
MeanSYFSCVFD2 <- data.frame(MeanSYFSCVFD2[complete.cases(MeanSYFSCVFD2),])

##Code to produce Fig. 2
my_list1[sapply(my_list1, is.null)] <- NULL

pdf(file="Plots/SYBoxPlots.pdf",  width=8, height=5)
boxplot(my_list1, ylab="Inter-taxa synchrony", xlab="Lakes")
abline(h=0, lwd=2, lty=2)
dev.off()

##Code to produce Fig. 4
# Load necessary libraries
library(lme4)
library(arm)

# Construct the data frame from the list
df <- do.call(rbind, Map(data.frame, my_list3))

# Calculate FS_Mean and FS_Deviation
df$FS_Mean <- rep(aggregate(df$FS, list(df$Lake), FUN=mean)[[2]], table(df$Lake))
df$FS_Deviation <- df$FS - df$FS_Mean

# Create a numeric LakeID
df$LakeID <- as.numeric(as.character(factor(df$Lake, levels = unique(df$Lake), labels = seq_along(unique(df$Lake)))))

# Fit the first linear mixed model for bloops
#bloop <- lmer(FS ~ 1 + (1 | Lake), data = df)
#bloops <- coefficients(bloop)
#df$FS_Bloop <- rep(bloops$Lake$`(Intercept)`, table(df$Lake))
#df$Bloop_Deviation <- df$FS - df$FS_Bloop

# Fit the second linear mixed model for Bloop_mod
#Bloop_mod <- lmer(SYRhoe ~ FS_Bloop + Bloop_Deviation + (Bloop_Deviation | Lake), data = df)
#summary(Bloop_mod)

# Fit the third linear mixed model for mod
mod <- lmer(SYRhoe ~ FS_Mean + FS_Deviation + (FS_Deviation | LakeID), data = df)
summary(mod)
smod <- sim(mod, 1000)

# Plot the coefficients of the model
#plot(coefficients(mod)$LakeID, pch = 19)

# Calculate confidence intervals
CI <- apply(smod@ranef$LakeID[,,2] + smod@fixef[,3], 2, quantile, c(0.025, 0.975))
CImean <- quantile(smod@fixef[,3], c(0.025, 0.975))

# Define the mean vector
mean <- data.frame(intercept = fixef(mod)[1], slope = fixef(mod)[3], lCI = quantile(smod@fixef[,3], 0.025), uCI = quantile(smod@fixef[,3], 0.975), col = "red")

# Create the data frame with intercept, slope, and confidence intervals
df2 <- data.frame(intercept = coefficients(mod)$LakeID[,1],
                  slope = coefficients(mod)$LakeID[,3],
                  lCI = CI[1,],
                  uCI = CI[2,],
                  col = "black")

# Combine the data frames and sort by slope
df2 <- rbind(df2, mean)
df2 <- df2[order(as.numeric(df2$slope)),]

# Display the first few rows of df2
head(df2)

# Plotting Fig. 3
pdf(file = "Plots/Slopes1.pdf", width = 6, height = 8)
plot(df2$slope, 1:nrow(df2), xlim = c(-0.5, 0.8), pch = 19, xlab = "Slope SY~FS", 
     ylab = "Lake ID", col = df2$col)
segments(df2$lCI, 1:nrow(df2), df2$uCI, 1:nrow(df2), col = df2$col)
abline(v = 0, lty = 2)
dev.off()

##Code to produce Fig. 4
pdf(file="Plots/CoV_SY_FS.pdf", width = 3, height = 8)
par(mfrow=c(3,1))
plot(MeanSYFSCVFD2$Stability~MeanSYFSCVFD2$SYRhoe, pch=19, main="Synchrony and Stability",
     xlab="Community Synchrony", ylab="Stability")
abline(lm(MeanSYFSCVFD2$Stability~MeanSYFSCVFD2$SYRhoe), lwd=2)


plot(MeanSYFSCVFD2$Stability~MeanSYFSCVFD2$FDis, pch=19, main="Functional Diversity and Stability",
     xlab="Mean FD", ylab="Stability")
abline(lm(MeanSYFSCVFD2$Stability~MeanSYFSCVFD2$FDis), lwd=2)


plot(MeanSYFSCVFD2$SYRhoe~MeanSYFSCVFD2$FDis, pch=19, main="Functional Diversity and Synchrony",
     xlab="Mean FD", ylab="Community Synchrony")
abline(lm(MeanSYFSCVFD2$SYRhoe~MeanSYFSCVFD2$FDis), lwd=2, lty=2)

dev.off()


##MODEL IS rho_x(k) = rho_e(k)*DC(k)

# Create a dataframe for plotting
contribution_data <- data.frame(
  LakeID = 1:n.lakes,
  EC = rho_e_contributions,
  ECCI_lower = rho_e_CI[,1],
  ECCI_upper = rho_e_CI[,2],
  DC = DC_contributions,
  DCCI_lower = DC_CI[,1],
  DCCI_upper = DC_CI[,2]
)

# Filter out rows with any missing values
contribution_data <- contribution_data[complete.cases(contribution_data), ]

# Adjust the values: add a sign column to help position the bars
contribution_data <- contribution_data %>%
  mutate(LakeID = factor(LakeID, levels = rev(unique(LakeID)))) # Reverse the order of LakeID

# Define a color palette
palette <- c("DC" = "#1f78b4", "EC" = "#33a02c")

# Create a separate dataframe for the error bars
error_data <- data.frame(
  LakeID = rep(contribution_data$LakeID, 2),
  Contribution = c(contribution_data$DC, contribution_data$EC),
  CI_lower = c(contribution_data$DCCI_lower, contribution_data$ECCI_lower),
  CI_upper = c(contribution_data$DCCI_upper, contribution_data$ECCI_upper),
  Type = rep(c("DC", "EC"), each = nrow(contribution_data))
)



# Plot the contributions

plot <- ggplot() +
  # Add points for the means
  geom_point(data = contribution_data, aes(x = DC, y = LakeID, color = "DC"), size = 3) +
  geom_point(data = contribution_data, aes(x = EC, y = LakeID, color = "EC"), size = 3) +
  # Add error bars for the confidence intervals
  geom_errorbarh(data = error_data, aes(xmin = CI_lower, xmax = CI_upper, y = LakeID, color = Type), height = 0.2) +
  scale_color_manual(values = palette) +
  labs(title = "Demographic and Environmental\nComponents of Synchrony",
       x = "Mean values",
       y = "Lake ID",
       color = NULL) + # Remove legend title
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1), # Keep y-axis labels horizontal
        plot.title = element_text(hjust = 0.5, face = "bold"),         # Center and bold the title
        legend.position = "top")                                      # Move legend to the top

ggsave("Plots/DC_EC_ForestPlot.jpeg", plot = plot, width = 4, height = 10)

# Assuming you have the necessary data:
# rho_e_contributions, DC_contributions, and rho_x values

# Calculate the relative contributions
relative_contributions <- data.frame(
  LakeID = 1:n.lakes,
  rho_x = rho_x_contributions, # Include rho_x values
  rho_e = rho_e_contributions / (rho_e_contributions + DC_contributions),
  DC = DC_contributions / (rho_e_contributions + DC_contributions)
)


# Filter out rows with any missing values
relative_contributions <- relative_contributions[complete.cases(relative_contributions), ]

# Melt the dataframe for easier plotting
relative_contributions_melt <- melt(relative_contributions, id.vars = c("LakeID", "rho_x"), variable.name = "Type", value.name = "Contribution")

# Flip the order of LakeID
relative_contributions_melt <- relative_contributions_melt %>%
  mutate(LakeID = factor(LakeID, levels = rev(1:n.lakes)))

# Define a color palette
palette <- c("DC" = "#1f78b4", "rho_e" = "#33a02c")

# Plot the relative contributions

plot <- ggplot(relative_contributions_melt, aes(y = LakeID, x = Contribution * rho_x, fill = Type)) + # Scale contribution by rho_x
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = palette) +
  labs(title = "Relative Contributions of rho_e\nand DC to rho_x",
       x = "Contribution to rho_x",
       y = "Lake ID",
       fill = NULL) + # Remove legend title
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1), # Keep y-axis labels horizontal
        plot.title = element_text(hjust = 0.5, face = "bold"),          # Center and bold the title
        legend.position = "top")                                        # Move legend to the top


ggsave("Plots/DC_EC_Contributions.jpeg", plot = plot, width = 4, height = 10)
