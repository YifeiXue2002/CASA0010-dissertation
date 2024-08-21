setwd("~/Desktop/Dissertation")
library("sf")
library("ggplot2")
library("spdep")
library("sp")
library("readxl")
library("spgwr")
library("car")
library(sjPlot)
library(ggspatial)

# load files
data <- read.csv("Data/imd_england_lsoa.csv")
data2<- read.csv("Data/LSOA2011 AvPTAI2015.csv")

# Merge the datasets on 'lsoa_code'
merged_df <- merge(data, data2, by = "lsoa_code", all = FALSE)
print(merged_df)

lsoashp <- read_sf("Data/statistical-gis-boundaries-london/ESRI/LSOA_2011_London_gen_MHW.shp")

# Merge the spatial data with the merged_df
spatialdatafile <- merge(lsoashp, merged_df, by.x = "LSOA11CD", by.y = "lsoa_code")

# Linear model
model1 <- lm(log10(AvPTAI2015) ~ IMD_decile + Income_decile + Education_decile, data = spatialdatafile)
tab_model(model1)

# Summary and VIF
options(scipen = 7)
summary(model1)
vif(model1)

# Extract residuals
spatialdatafile$RESIDUALS <- model1$residuals

# Map residuals
ggplot(data = spatialdatafile) +
  geom_sf(aes(fill = RESIDUALS), color = "black", alpha = 0.7) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, 
                       breaks = quantile(spatialdatafile$RESIDUALS, probs = seq(0, 1, 0.2), na.rm = TRUE),
                       labels = scales::percent_format()) +
  theme_minimal() +
  theme(legend.position = c(0.7, 0.85),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  labs(fill = "Residuals") +
  coord_sf()

# Spatial autocorrelation
spatialdatafile$ROWNUM <- 1:nrow(spatialdatafile)
spatialdatafile_2.0 <- as(spatialdatafile, "Spatial")
Weights <- poly2nb(spatialdatafile_2.0, row.names = spatialdatafile_2.0$ROWNUM)
WeightsMatrix <- nb2mat(Weights, style='B', zero.policy=TRUE)
Residual_WeightMatrix <- mat2listw(WeightsMatrix, style='W', zero.policy=TRUE)
lm.morantest(model1, Residual_WeightMatrix, alternative="two.sided")

# Transform and calculate centroids
spatialdatafile <- st_transform(spatialdatafile, 27700)
spatialdatafile <- st_centroid(spatialdatafile)
spatialdatafile <- cbind(spatialdatafile, st_coordinates(spatialdatafile))

# GWR bandwidth selection
BwG <- gwr.sel(log10(AvPTAI2015) ~ IMD_decile + Income_decile + Education_decile, 
               data = spatialdatafile, 
               coords = cbind(spatialdatafile$X, spatialdatafile$Y), 
               adapt = TRUE)

# GWR model
start.timer <- proc.time()
gwr.model <- gwr(log10(AvPTAI2015) ~ IMD_decile + Income_decile + Education_decile, 
                 data = spatialdatafile, 
                 coords = cbind(spatialdatafile$X, spatialdatafile$Y), 
                 adapt = BwG, 
                 hatmatrix = TRUE, 
                 se.fit = TRUE)
end.timer <- proc.time() - start.timer
end.timer

gwr.model

# GWR results
gwr.data <- as.data.frame(gwr.model$SDF)

# Create spatial data frame for results
lsoa_result <- lsoashp

# Insert coefficients and standard errors
lsoa_result$CoefIMD <- gwr.data[,"IMD_decile"]
lsoa_result$CoefIncome <- gwr.data[,"Income_decile"]
lsoa_result$CoefEducation <- gwr.data[,"Education_decile"]

lsoa_result$SEIMD <- gwr.data[,"IMD_decile_se"]
lsoa_result$SEIncome <- gwr.data[,"Income_decile_se"]
lsoa_result$SEEducation <- gwr.data[,"Education_decile_se"]

# Insert localR2 estimates
lsoa_result$localR2 <- gwr.data[,"localR2"]

# Function to create a map
create_map <- function(data, var, title, midpoint = 0, palette = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = midpoint)) {
  ggplot() +
    geom_sf(data = data, aes(fill = {{var}}), color = "black", size = 0.1) +
    palette +
    theme_minimal() +
    theme(
      legend.position = c(0.9, 0.93),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title, x = "", y = "") +
    coord_sf()
}

# Visualize coefficients
create_map(lsoa_result, CoefIncome, "Coefficient: Income Decile")
summary(lsoa_result$CoefIncome)

create_map(lsoa_result, CoefEducation, "Coefficient: Education Decile")
summary(lsoa_result$CoefEducation)

create_map(lsoa_result, CoefIMD, "Coefficient: IMD Decile")
summary(lsoa_result$CoefIMD)

# Arrange plots
grid.arrange(p4, p5, p6, ncol = 3)

# Compute t-score statistics
lsoa_result$tstatIncome <- lsoa_result$CoefIncome / lsoa_result$SEIncome
lsoa_result$tstatEducation <- lsoa_result$CoefEducation / lsoa_result$SEEducation
lsoa_result$tstatIMD <- lsoa_result$CoefIMD / lsoa_result$SEIMD

# Create significance columns
lsoa_result$significantIncome <- cut(lsoa_result$tstatIncome,
                                     breaks = c(-Inf, -2, 2, Inf),
                                     labels = c("Reduction: Significant", "Not Significant", "Increase: Significant"))

lsoa_result$significantEducation <- cut(lsoa_result$tstatEducation,
                                        breaks = c(-Inf, -2, 2, Inf),
                                        labels = c("Reduction: Significant", "Not Significant", "Increase: Significant"))

lsoa_result$significantIMD <- cut(lsoa_result$tstatIMD,
                                  breaks = c(-Inf, -2, 2, Inf),
                                  labels = c("Reduction: Significant", "Not Significant", "Increase: Significant"))
# Plot for Significance: Income
ggplot(data = lsoa_result) +
  geom_sf(aes(fill = significantIncome), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("blue", "white", "red"), name = "Significance: Income") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.93),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5)) +
  labs(fill = "Significance: Income") +
  coord_sf() +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.5)

# Plot for Significance: Education
ggplot(data = lsoa_result) +
  geom_sf(aes(fill = significantEducation), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("blue", "white", "red"), name = "Significance: Education") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.93),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5)) +
  labs(fill = "Significance: Education") +
  coord_sf() +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.5)

# Plot for Significance: IMD
ggplot(data = lsoa_result) +
  geom_sf(aes(fill = significantIMD), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("blue", "white", "red"), name = "Significance: IMD") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.93),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5)) +
  labs(fill = "Significance: IMD") +
  coord_sf() +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.5)

# Arrange plots
grid.arrange(p7, p8, p9, ncol = 3)

#local R2
ggplot() +
  geom_sf(data = lsoashp, fill = NA, color = "black", size = 0.1) +
  geom_sf(data = lsoa_result, aes(fill = localR2), color = NA) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    name = "Adaptive: Local R2"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.9, 0.9),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Adaptive: Local R2", x = "", y = "") +
  coord_sf()
