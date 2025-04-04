####CODE for OLS and SARs for Bioregions
#Used in: "The latitudinal variation in amphibian speciation rates revisited"
#Author: Adrian Garcia Rodriguez
#Last update: March 24, 2025
rm(list=ls())
gc()

##Load packages
PKGs <- c("raster","rlang", "dplyr","rgdal","maps", "maptools","sp","spocc","plyr", "tidyverse", "pbapply", "sf", "ggplot2", "rgeos", "R.utils","data.table", "here",  "RColorBrewer", "rasterVis", "beepr", "stringr", "corrplot", "AER", "sjPlot", "sjmisc", "sjlabelled", "ggpubr", "blmeco","partR2", "magrittr", "glmm.hp", "GLMMadaptive", "gridExtra", "lme4", "furrr", "spdep", "ncf","spdep", "flextable", "stargazer")
sapply(PKGs, require, character.only = TRUE)

setwd("C:/Users/Lenovo/Dropbox/PANDEMIC PRIORITIES/5.LGD_Speciation/THIRD TRY/Docs_Submission_NEE/Codes/")
setwd("C:/Users/garciaa/Dropbox/PANDEMIC PRIORITIES/5.LGD_Speciation/THIRD TRY/Docs_Submission_NEE/Codes/")
setwd("C:/Users/garciaa/Dropbox/MANUSCRIPTS/2025/LGSpeciation/Codes/GitHub")

setwd("C:/Users/Administrator/Dropbox/MANUSCRIPTS/2025/LGSpeciation/")

##Load environment
#load("Codes/SARs_preds_bioregs.RData")

#Load data 
BioregData<- read.csv("Data/DataBioregs_Spec_Preds_FINAL.csv", sep=",")
BioregData2<-BioregData[ ,c(2:5,7,9,10, 11:14,16,18)]
names(BioregData2)
lm_spec_rich<- lm(Richn1dg ~ Mean.Spec.1dg, data=BioregData2)
summary(lm_spec_rich)

##Test for multicolinearity among predictors
preds<-BioregData2[, c(6:13)] ##Choose numeric predictors
df_new <- preds[, order(colnames(preds))]

df_new<-round(df_new,3)
df_new<-na.omit(df_new)

##Pairwise correlations
p.mat <- cor(df_new)

##Plot
corrplot(p.mat,title = "Correlations among Predictors", method ="color", type = "lower", diag = FALSE,addCoef.col = "grey20",#
         tl.col = "grey20",  order = "hclust",hclust.method = "ward.D", number.cex = 1,tl.cex = 1, cl.cex=1, outline = F, mar=c(0,0,0,8), sig.level = 0.01)


########################################################################################
########################################################################################
###Scale data
BioregData_sc<-cbind(BioregData2[,1:5], scale(BioregData2[,6:13]))
hist(BioregData_sc$Mean.Spec.1dg) ##response with normal distribution

####OLS: Spec Rate as Response and predictors scaled
olsfull<-lm(Mean.Spec.1dg ~  Time.Area + Area.Product + Temp  + Rough1dg + NRI1dg + ClimV1dg, data=BioregData_sc)
summary(olsfull)
res.ols <- residuals(olsfull)
fitted.ols<-fitted(olsfull)

###TEST OF LEVERAGE (considering the regression Spec~TopoComplex)
ols<-lm(Mean.Spec.1dg ~  Rough1dg, data=BioregData_sc)
summary(ols)
res.ols <- residuals(ols)
fitted.ols<-fitted(ols)

# Extract leverage values using hatvalues()
hats <- hatvalues(ols)

# Convert the leverage values into a data frame
hats_df <- data.frame(Bioregion = BioregData$Bioregion, Leverage = hats)
# Order the Bioregion column alphabetically
hats_df$Bioregion <- factor(hats_df$Bioregion, levels = sort(unique(hats_df$Bioregion)))

# Calculate the cutoff value
avglv_low <- 2 * (7 / 32)
avglv_high <- 3 * (7 / 32)

# Create a horizontal bar plot
library(ggplot2)

ggplot(hats_df, aes(x = Leverage, y = Bioregion)) + 
  geom_bar(stat = "identity", fill = "cadetblue4") +  # Bars for leverage values
  geom_vline(xintercept = avglv_low, linetype = "dashed", color = "orange", linewidth = .25) +  # Vertical dashed line at cutoff
  geom_vline(xintercept = avglv_high, linetype = "dashed", color = "red", linewidth = .25) +  # Vertical dashed line at cutoff
  geom_text(aes(x = avglv_low + 0.005, y = max(as.numeric(Bioregion)) - 1), 
            label =round(avglv_low, 3), 
            color = "black", hjust = 0, , size = 3) +  # Add text label next to dashed line
  geom_text(aes(x = avglv_high + 0.005, y = max(as.numeric(Bioregion)) - 1), 
            label = round(avglv_high, 3), 
            color = "black", hjust = 0, size = 3) +  # Add text label next to dashed line
  labs(title = "Leverage Values with low and high cutoffs", x = "Leverage Value", y = "Bioregion", size=3) +
  theme_pubr() +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 7)) +
  xlim(0, 0.8) # Adjust y-axis text size for better visibility
###Plot saved as PNG

### Defined neighborhood structure
nb <- spdep::dnearneigh(x=BioregData_sc[,2:3], d1 = 0, d2 = 50, longlat = FALSE)# this create a spatial neighbor objects based on distance thresholds
plot(nb, BioregData_sc[,2:3])#visualize neighbor relationships
str(nb, list.len=5, give.attr = F)
lw <-nb2listw(nb, style = "B", zero.policy = TRUE)#convert the nb object into spatial weights matrix
str(lw$weights, list.len=5, give.attr = F)

#####MORAN TEST (to detect spatial autocorrelation)
MORAN<-lm.morantest(ols, lw, zero.policy = F)
MORAN$p.value ###p<5.60e-16 seems to be spatially autocorrelated

#####RAO TEST (Rao’s score test allow a distinction between spatial error models and spatial lag models)
lm.RStests(ols, listw=lw, test="all")
#
##Results indicate that the spatial structure in the error terms is highly significant

#lag.listw() is used to compute the spatial lag of a variable based on a spatial weights matrix (listw). The spatial lag represents a weighted average of the neighboring values for each observation, where the weights are defined by the spatial relationships specified in the weights matrix.What is a Spatial Lag?A spatial lag is a transformation that allows you to include the influence of neighboring observations in your model. It provides a way to quantify how the values of a variable at neighboring locations affect the value at a particular location

### MORAN Plot
moran.plot(res.ols, lw, zero.policy = FALSE, labels = TRUE)##shows a positive spatial autocorrelation

###Perform a Monte Carlo simulation to calculate a randomization for the p-value
moran <- moran.mc(res.ols,lw, nsim = 9999, zero.policy = TRUE)
moran #confirms significance
plot(moran)

###RUN LAG AND ERROR MODELS AND COMPARE THEM
# define spatial neighborhoods - nb object
library(pbapply)
###Function to run the SARs with varying distances defining the neighbor structure
SARsVarDists<- pblapply(10:100, function(x){
  spec.nb <- dnearneigh(BioregData_sc[,2:3], d1 = 0, d2 =  x, longlat = FALSE)
  
  # define spatial weights for neighborhoods - listw object
  spec.lw <- nb2listw(spec.nb, style = "B", zero.policy = TRUE)
  
  # fit spatial lag model
  spec.slm <- spatialreg::lagsarlm(Mean.Spec.1dg ~ Time.Area +  Area.Product + Temp  +  NRI1dg + Rough1dg + ClimV1dg, data = BioregData_sc, listw = spec.lw, zero.policy = TRUE)
  #summary(spec.slm)
  
  # fit spatial error model
  spec.sem <- spatialreg::errorsarlm(Mean.Spec.1dg ~ Time.Area +  Area.Product + Temp  +  NRI1dg + Rough1dg + ClimV1dg, data=BioregData_sc, listw = spec.lw, zero.policy = TRUE)
  #summary(spec.sem, Nagelkerke=TRUE)
  
  AICs<-AIC(spec.sem, spec.slm) 
  AICs$dist<-x
  return(AICs) ###Chosen model error model d2=99 (longlat= FALSE)
  
})

SARsVarDistsDF<- do.call(rbind,SARsVarDists)

# Select the model with lowest AIC
mod.lowestAIC <- SARsVarDistsDF[which.min(SARsVarDistsDF$AIC),]
optimal_dist<-mod.lowestAIC$dist

####RUN MODEL WITH THE OPTIMAL DISTANCE
spec.nb <- dnearneigh(BioregData_sc[,2:3], d1 = 0, d2 =  optimal_dist, longlat = FALSE)
plot(spec.nb, BioregData_sc[,2:3])#visualize neighbor relationships

# define spatial weights for neighborhoods - listw object
spec.lw <- nb2listw(spec.nb, style = "B", zero.policy = TRUE)

spec.semFINAL <- spatialreg::errorsarlm(Mean.Spec.1dg ~ Time.Area +  Product + Temp  +  NRI1dg + Rough1dg + ClimV1dg, data=BioregData_sc, listw = spec.lw, zero.policy = TRUE)
summary(spec.semFINAL, Nagelkerke=TRUE)

coefficients_sar <- summary(spec.semFINAL)$coefficients

# Assuming predictors in your model are named predictor1, predictor2, etc.
predictors <- names(coefficients_sar)[-1]  # Exclude the intercept

library(ggplot2)
# For each predictor, create a plot showing the regression line
for (predictor in predictors) {
  
  # Create a data frame with the predictor and dependent variable
  plot_data <- BioregData %>%
    dplyr::select(Mean.Spec.1dg, all_of(predictor))
  
  # Calculate the predicted values for the SAR model using the extracted coefficient
  plot_data$predicted <- coefficients_sar[predictor] * plot_data[[predictor]] + coefficients_sar["(Intercept)"]
  
  # Plot the regression with ggplot2
  p <- ggplot(plot_data, aes_string(x = predictor, y = "Mean.Spec.1dg")) +
    geom_point(color="dodgerblue4", size=4, alpha=.5) +  # Scatter plot of the original data
    geom_line(aes_string(y = "predicted"), color = "red", size = 1) +  # Regression line
    labs(
      title = paste("Regression of", predictor, "from SAR model"),
      x = predictor,
      y = "Mean DR"
    ) +
    theme_bw()
  
  # Print the plot for each predictor
  print(p)
}


####PRESENT RESULTS AS FOREST PLOT
# Get the summary of the model
model_summary <- summary(spec.semFINAL, Nagelkerke=TRUE)
str(model_summary)
model_summary$Coef

# Get the coefficients (point estimates)
estimates <- model_summary$Coef[-1, 1] * 1000 # Point estimates (excluding the intercept)

# Get the standard errors of the coefficients
std_errors <- model_summary$Coef[-1, 2] * 1000  # Standard errors (excluding the intercept)

# Calculate the 95% confidence intervals (Estimate ± 1.96 * SE)
lower_ci <- estimates -  std_errors  # Lower bound of CI
upper_ci <- estimates +  std_errors  # Upper bound of CI

# Prepare data for the forest plot
plot_data <- data.frame(
  variable = names(estimates),
  estimate = estimates,
  lower = lower_ci,
  upper = upper_ci
)

# Create a data frame for ggplot
gg_data <- data.frame(
  variable = plot_data$variable,
  estimate = plot_data$estimate,
  lower = lower_ci,
  upper = upper_ci
)

##Add proper names to the variables
gg_data$variable<- c("Time-Area", "Productivity", "Temperature", "Net relatedness index", "Topographic complexity", "Climatic velocity")

# Add a color column, default to gray, then change specific variable colors
gg_data$color <- ifelse(gg_data$variable %in% c("Climatic velocity", "Topographic complexity"), "dodgerblue4", "gray70")
gg_data <- gg_data[order(gg_data$estimate), ]

# Create the plot with ggplot
library(ggpubr)

num_rows <- nrow(gg_data)

# Create the plot with alternating gray and white bands
forest <- ggplot(gg_data, aes(x = estimate, y = variable)) +
  # Add alternating background bands with more control
  geom_rect(
    aes(
      ymin = as.numeric(reorder(variable, estimate)) - 0.5,
      ymax = as.numeric(reorder(variable, estimate)) + 0.5,
      xmin = -Inf, 
      xmax = Inf,
      fill = ifelse(as.numeric(reorder(variable, estimate)) %% 2 == 0, "even", "odd")
    ),
    alpha = 0.3
  ) +
  scale_fill_manual(
    values = c("odd" = "gray85", "even" = "white"),
    guide = "none"
  ) +
  
  # Error bars
  geom_errorbarh(
    aes(xmin = lower, xmax = upper), 
    color = "gray50", 
    height = 0, 
    linewidth = 0.75
  ) +  
  
  # Points with proper color mapping
  geom_point(
    aes(color = color), 
    size = 3.5
  ) +
  scale_color_identity() +  # Use colors as-is from the data
  
  # Reference line
  geom_vline(
    xintercept = 0, 
    linetype = "dashed", 
    color = "gray30", 
    linewidth = 0.75
  ) +
  
  # Axis scaling
  scale_x_continuous(limits = c(-7, 7)) +
  
  # Theme
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 11),    
    axis.text.y = element_text(
      size = 12,
      margin = margin(r = 20)  # Add space between y-axis labels and plot
    ), 
    axis.text.x = element_text(size = 9),
    plot.title = element_text(size = 18,  hjust = 0.5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.border = element_rect(color = "gray80", fill = NA)  # Lighter panel border
  ) +
  
  labs(
    x = "Effect sizes", 
    y = NULL   # More elegant than empty string
  )

# Display the plot
print(forest)

library(grid)


predictor= 'ClimV1dg'
# Create a data frame with the predictor and dependent variable
plot_data <- BioregData_sc %>%
  dplyr::select(Mean.Spec.1dg, all_of(predictor))

# Calculate the predicted values for the SAR model using the extracted coefficient
plot_data$predicted <- coefficients_sar[predictor] * plot_data[[predictor]] + coefficients_sar["(Intercept)"]

# Plot the regression with ggplot2
p.clim <- ggplot(plot_data, aes_string(x = predictor, y = "Mean.Spec.1dg")) +
  geom_point(color="dodgerblue4", size=4, alpha=.5) +  # Scatter plot of the original data
  geom_line(aes_string(y = "predicted"), color = "red", size = 1) +  # Regression line
  labs(
    title = paste("Regression of", predictor, "from SAR model"),
    x = predictor,
    y = "Mean DR"
  ) +
  theme_bw()

predictor= 'Rough1dg'
# Create a data frame with the predictor and dependent variable
plot_data <- BioregData_sc %>%
  dplyr::select(Mean.Spec.1dg, all_of(predictor))

# Calculate the predicted values for the SAR model using the extracted coefficient
plot_data$predicted <- coefficients_sar[predictor] * plot_data[[predictor]] + coefficients_sar["(Intercept)"]

# Plot the regression with ggplot2
p.rough <- ggplot(plot_data, aes_string(x = predictor, y = "Mean.Spec.1dg")) +
  geom_point(color="dodgerblue4", size=4, alpha=.5) +  # Scatter plot of the original data
  geom_line(aes_string(y = "predicted"), color = "red", size = 1) +  # Regression line
  labs(
    title = paste("Regression of", predictor, "from SAR model"),
    x = predictor,
    y = "Mean DR"
  ) +
  theme_bw()

ggarrange(forest,p.clim,p.rough,nrow = 3, ncol = 1)


