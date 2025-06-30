# Modern GLM Analysis for Species Distribution Modeling
# Comprehensive version integrating stepwise selection and advanced plotting
# Updated to address deprecated functions and provide coherent analysis

# Load required libraries
library(biomod2)
library(MASS)
library(ggplot2)
library(gridExtra)  # For arranging multiple plots
library(dplyr)      # For data manipulation
library(viridis)    # For better color palettes
library(lattice)    # For additional plotting options

# Load the dataset
#mammals_data <- read.csv("data/tabular/species/mammals_and_bioclim_table.csv", row.names = 1)
ddia_adjsd_bpi <- read.csv("data/ne_canyons_dscs_btm_20250630.csv", row.names = 1)

# Data exploration and summary
# cat("=== Dataset Summary ===\n")
# cat("Total observations:", nrow(mammals_data), "\n")
# cat("Red fox presence records:", sum(mammals_data$VulpesVulpes == 1), "\n")
# cat("Red fox absence records:", sum(mammals_data$VulpesVulpes == 0), "\n")
# cat("Prevalence:", round(mean(mammals_data$VulpesVulpes), 3), "\n\n")

# Data exploration and summary
cat("=== Dataset Summary ===\n")
cat("Total observations:", nrow(ddia_adjsd_bpi), "\n")
cat("Desmophyllum dianthus presence records:", sum(ddia_adjsd_bpi$Desmophyllum.dianthus == 1), "\n")
cat("Desmophylum dianthus absence records:", sum(ddia_adjsd_bpi$Desmophyllum.dianthus == 0), "\n")
cat("Prevalence:", round(mean(ddia_adjsd_bpi$Desmophyllum.dianthus), 3), "\n\n")

# Check for missing values
# missing_summary <- colSums(is.na(mammals_data))
# if(any(missing_summary > 0)) {
#   cat("Missing values found:\n")
#   print(missing_summary[missing_summary > 0])
# } else {
#   cat("No missing values detected.\n")
# }

# Check for missing values
missing_summary <- colSums(is.na(ddia_adjsd_bpi))
if(any(missing_summary > 0)) {
  cat("Missing values found:\n")
  print(missing_summary[missing_summary > 0])
} else {
  cat("No missing values detected.\n")
}

### 10.3	Generalized Additive Models


```{r code_10.3_Generalized_Additive_Models_b1}
if(is.element("package:mgcv", search())) detach("package:mgcv") ## make sure the mgcv package is not loaded to avoid conflicts
library(gam)
gam1 = gam(VulpesVulpes ~ s(bio3,2) + s(bio7,2) + s(bio11,2) + s(bio12,2), data=mammals_data, family="binomial")
gam2 = gam(VulpesVulpes ~ s(bio3,4) + s(bio7,4) + s(bio11,4) + s(bio12,4), data=mammals_data, family="binomial")
```


```{r GAMb2 10.6, opts.label = 'fig_half_page', fig.cap = 'Figure 10. 6: Response curves of model gam1 expressed in logit scale (function plot.gam() from the gam package).'}
par(mfrow=c(2,2))
plot(gam1, se=T)
```


```{r GAMb3 10.7, message=FALSE,warning=FALSE, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10. 7: Response curves of the gam1 (degree of smoothing = 2) and gam2 (degree of smoothing = 4) models.'}
rp <- response.plot2(models = c('gam1', 'gam2'),
                     Data = mammals_data[,c("bio3", "bio7", "bio11", "bio12")],
                     show.variables = c("bio3",  "bio7", "bio11", "bio12"),
                     fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

gg.rp <- ggplot(rp, aes(x = expl.val, y = pred.val, lty = pred.name)) +
  geom_line() + ylab("prob of occ") + xlab("") + 
  rp.gg.theme + 
  facet_grid(~ expl.name, scales = 'free_x')
print(gg.rp)
```

```{r code_10.3_Generalized_Additive_Models_b4}
gamStart <- gam(VulpesVulpes~1, data=mammals_data, family=binomial)
gamModAIC <- step.gam(gamStart, biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 4), trace=F, direction = "both")        
```

```{r code_10.3_Generalized_Additive_Models_b4b, echo=FALSE, eval = FALSE}
## test of the multiple smoother case
biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 2:4)
gamModAIC.test <- step.gam(gamStart, biomod2:::.scope(mammals_data[1:3,c("bio3",  "bio7", "bio11", "bio12")], "s", 2:4), trace=F, direction = "both")
```


```{r GAMb5 10.8, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10. 8. Observed (black=presence, light gray= absence) and potential distribution of Vulpes vulpes extracted from gamModAIC. The gray scale of prediction illustratesshows habitat suitability values between 0 (light, unsuitable) and 1 (dark, highly suitable)'}
par(mfrow=c(1,2))
level.plot(mammals_data$VulpesVulpes, XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")
level.plot(fitted(gamModAIC), XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="Stepwise GAM with AIC")
```



```{r code_10.3_Generalized_Additive_Models_b6, fig.keep = FALSE}
if(is.element("package:gam", search())) detach("package:gam") ## make sure the gam package is not loaded to avoid conflicts
library(mgcv)
gam_mgcv <- gam(VulpesVulpes~s(bio3)+s(bio7)+s(bio11)+s(bio12),data = mammals_data, family = "binomial")
## see a range of summary statistics
summary(gam_mgcv)
gam.check(gam_mgcv)
```


```{r GAMb7 10.9, opts.label = 'fig_half_page', fig.cap = 'Figure 10.9. Response curves of model gam_mgcv plotted using the internal function of mgcv().'}
plot(gam_mgcv,pages=1, seWithMean=TRUE)  
```


```{r GAMb8 10.10, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10.10. Response curves from the model calibrated with the mgcv package (gam_mgcv).'}
rp <- response.plot2(models = c('gam_mgcv'),
                     Data = mammals_data[,c("bio3", "bio7", "bio11", "bio12")],
                     show.variables = c("bio3",  "bio7", "bio11", "bio12"),
                     fixed.var.metric = 'mean', plot = FALSE, use.formal.names = TRUE)

gg.rp <- ggplot(rp, aes(x = expl.val, y = pred.val, lty = pred.name)) +
  geom_line() + ylab("prob of occ") + xlab("") + 
  rp.gg.theme + 
  facet_grid(~ expl.name, scales = 'free_x')
print(gg.rp)
```

```{r GAMb9 10.11, opts.label = 'fig_quarter_page', fig.cap = 'Figure 10.11. Observed (black=presence, light gray= absence) and potential distribution of Vulpes vulpes extracted from the gam_mgcv object. The gray scale of predictions illustrates shows habitat suitability values between 0 (light, unsuitable) and 1 (dark, highly suitable).'}
par(mfrow=c(1,2))
level.plot(mammals_data$VulpesVulpes, XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3,level.range=c(0,1), show.scale=F, title="Original data")
level.plot(fitted(gam_mgcv), XY=mammals_data[,c("X_WGS84", "Y_WGS84")], color.gradient = "grey", cex=0.3, level.range=c(0,1), show.scale=F, title="GAM with mgcv")
par(mfrow=c(1,1))
```