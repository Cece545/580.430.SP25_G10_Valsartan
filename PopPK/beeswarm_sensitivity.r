library(R.matlab)
library(ggbeeswarm)
library(ggplot2)
library(tidyr)
library(gridExtra)

# Set working directory to the location of the .mat file
setwd("/Users/ceceliazhang/Desktop/JHU/SP25/Systems_Pharmacology_Personalized_Medicine/Project")

# Load the .mat file
mat_data <- readMat("valsartan_relative_sensitivity.mat")

# Extract data
data <- list(
  AUC_receptor = as.vector(mat_data$S.AUC.receptor),
  AUC_weight = as.vector(mat_data$S.AUC.weight),
  Cmax_receptor = as.vector(mat_data$S.Cmax.receptor),
  Cmax_weight = as.vector(mat_data$S.Cmax.weight),
  Cmin_receptor = as.vector(mat_data$S.Cmin.receptor),
  Cmin_weight = as.vector(mat_data$S.Cmin.weight)
)

# Convert to data frame
data_df <- do.call(cbind, data)
colnames(data_df) <- names(data)
data_df <- as.data.frame(data_df)

# Reshape data for plotting
data_long <- pivot_longer(data_df, cols = everything(), names_to = "Parameter", values_to = "Value")

# Function to create beeswarm plot without 'Parameter' text but keeping x-axis labels
generate_beeswarm <- function(data, parameter) {
  ggplot(data[data$Parameter == parameter, ], aes(x = Parameter, y = Value)) +
    geom_beeswarm(aes(color = ifelse(grepl('weight', Parameter), 'Weight', 'Receptor'))) +
    scale_color_manual(values = c('Weight' = '#56B4E9', 'Receptor' = '#E69F00')) +
    theme_minimal(base_size = 15) +
    ylim(-1.25, 0.5) +
    labs(x = NULL, y = "Value") +  # Remove x-axis label text
    theme(
      legend.position = 'none',
      plot.title = element_blank(),
      axis.title.x = element_blank(),  # Remove 'Parameter' text label
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1)
    )
}

# Define the receptor and weight parameter sets
receptor_params <- c('AUC_receptor', 'Cmax_receptor', 'Cmin_receptor')
weight_params <- c('AUC_weight', 'Cmax_weight', 'Cmin_weight')

# Generate subplots for receptor parameters
png('beeswarm_receptor_group.png', width = 12, height = 6, units = 'in', res = 300)
receptor_plots <- lapply(receptor_params, function(param) generate_beeswarm(data_long, param))
grid.arrange(grobs = receptor_plots, ncol = 3)
dev.off()

# Generate subplots for weight parameters
png('beeswarm_weight_group.png', width = 12, height = 6, units = 'in', res = 300)
weight_plots <- lapply(weight_params, function(param) generate_beeswarm(data_long, param))
grid.arrange(grobs = weight_plots, ncol = 3)
dev.off()
