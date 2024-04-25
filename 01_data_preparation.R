################################################################################
# IMPORTANT: THE MORTALITY DATA provided are simulated FROM MODEL 4 (SB-DLNM 
# WITH TIME-SERIES DESIGN) IN THE ORIGINAL SB-DLNM PAPER (https://doi.org/10.1093/ije/dyae061) 
# due to confidentiality restrictions of the original dataset. Original data 
# may be requested to the data owner (The Barcelona Public Health Agency) who 
# should share them in similar terms than those applying to this study. The 
# supplied mortality, with exactly the same structure than the original dataset, 
# allows reproducing the code provided so that interested readers may run it as 
# an example of use. However, the results from the original data are 
# additionally supplied (input/result_inla folder). We provide here the 
# implementation of B-DLNMs and SB-DLNMS in R-INLA, a significantly faster and 
# user-friendly software compared to WinBUGS.
################################################################################

################################################################################
# In this R project we implement Bayesian and spatial Bayesian distributed 
# lag non-Linear models (B-DLNM and SB-DLNM) for the case study of short-term 
# associations between temperature and mortality in the city of Barcelona.
################################################################################

################################################################################
# CODE 1: DATA PREPARATION
# Data preparation: Transforming time-series temperature and mortality datasets
# for B-DLNMs and SB-DLNMs.
################################################################################

# Load libraries

library(sf)
library(spdep)
library(lubridate)
library(dlnm)

# Load data

load("input/daily_data.RData")
shapefile_bcn <- read_sf("input/shapefile_bcn.shp")

# Subset data to only summer months of 2007 to 2016

data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

# Generate a list of spatial structure from the shapefile for use in INLA.
list_neig <- nb2mat(poly2nb(shapefile_bcn), style = "B")
save(list_neig, file = "output/list_neighbours.RData")

# Set variables defining the dlnm model

dlnm_var <- list(
  var_prc = c(0.50, 0.90),
  var_fun = "ns",
  max_lag = 8,
  lagnk = 2,
  n_reg = 73,
  n_coef = 12)

# Set variables for trend and seasonality

df_seas <- 4
df_trend_10years <- 1 # 1 df every 10 years to control for long-term trends
df_trend <- round(length(min(year(data$date)):max(year(data$date))) / df_trend_10years / 10) # Here we assume the time period for all regions is the same
rm(df_trend_10years)

dlnm_var$df_seas <- df_seas
dlnm_var$df_trend <- df_trend

rm(df_seas, df_trend)

save(dlnm_var, file = "output/dlnm_configuration.RData")

# Create crossbasis for each region

# Ensure that the data is ordered by region to maintain alignment between 
# crossbasis and regions
if(is.unsorted(data$region)) {
  stop("data needs to be ordered by region for the next loop")}

list_cb <- lapply(formatC(1:dlnm_var$n_reg, width = 2, flag = "0"), 
    function(i_reg) {
      
      temp <- subset(data, region == i_reg, 
                     select = c("temp", paste0("lag", 1:dlnm_var$max_lag)))
      
      temp_knots <- quantile(temp$temp, dlnm_var$var_prc, na.rm = TRUE)
      temp_boundary <- range(temp, na.rm = TRUE)
      
      cb <- crossbasis(temp,
                       argvar = list(fun = dlnm_var$var_fun,
                                     knots = temp_knots,
                                     Boundary.knots = temp_boundary),
                       arglag = list(fun = "ns",
                                     knots = logknots(dlnm_var$max_lag, 
                                                      dlnm_var$lagnk),
                                     intercept = TRUE))
      
      cb
      
})
cb <- do.call(rbind, list_cb)
rm(list_cb)

# Create different objects for the case-crossover and time-series designs

data_cco <- data; cb_cco <- cb
data_ts <- data; cb_ts <- cb

rm(data, cb)

# PREPARE DATA FOR THE CASE-CROSSOVER DESIGN

# Create strata for the case-crossover
# (neighborhood - year - month - day of week)
data_cco$strata <- paste(data_cco$region, 
                         year(data_cco$date), 
                         formatC(month(data_cco$date), width = 2, flag = "0"),
                         wday(data_cco$date, week_start = 1),
                         sep = ":")

# We exclude stratums without cases as they do not contribute to the 
# likelihood in the case-crossover design
keep <- sapply(split(data_cco$mort, data_cco$strata), sum)
keep <- data_cco$strata %in% names(keep[keep != 0])
data_cco <- data_cco[keep,]
cb_cco <- cb_cco[keep,]
rm(keep)

save(data_cco, file = "output/data_casecrossover.RData")
save(cb_cco, file = "output/crossbasis_casecrossover.RData")

# PREPARE DATA FOR THE TIME-SERIES DESIGN

# Define day of week and year
data_ts$day_of_week <- wday(data_ts$date, week_start = 1)
data_ts$year <- year(data_ts$date)

# Create the trend matrix

trend <- onebasis(data_ts$date, fun = "ns", df = dlnm_var$df_trend)

# Create the seasonality matrix

# Ensure that the data is ordered by region and year to maintain alignment 
# between seasonality and data
if(is.unsorted(order(data_ts$region, data_ts$year))) {
  stop("data needs to be ordered by region and year for the next loop")}

seas <- lapply(formatC(1:dlnm_var$n_reg, width = 2, flag = "0"), 
 function(i_reg) {
   lapply(unique(year(data_ts$date)), function(i_year) {
     date <- subset(data_ts, region == i_reg & year(date) == i_year)$date
     seas <- onebasis(yday(date), fun = "ns", df = dlnm_var$df_seas)
     return(seas)
   })
 })
seas <- do.call(rbind, do.call(rbind, seas)) 

save(data_ts, file = "output/data_timeseries.RData")
save(cb_ts, file = "output/crossbasis_timeseries.RData")
save(trend, file = "output/trend_timeseries.RData")
save(seas, file = "output/seasonality_timeseries.RData")