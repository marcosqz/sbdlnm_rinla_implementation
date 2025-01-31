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
# We used R version 4.4.2 and INLA version 24.12.11
################################################################################

################################################################################
# In this R project we implement Bayesian and spatial Bayesian distributed 
# lag non-Linear models (B-DLNM and SB-DLNM) for the case study of short-term 
# associations between temperature and mortality in the city of Barcelona.
################################################################################

################################################################################
# CODE 2: RUN B-DLNMs AND SB-DLNMs
# Model specification and execution: independent B-DLNMs and SB-DLNMs in R and
# R-INLA with simulated mortality data.
################################################################################

# Load libraries
library(INLA)

# 10 threads were used to generate the results
# Use "inla.setOption("num.threads", 10)" only if your computer has at least 
# 10 threads

# inla.setOption("num.threads", 10)

# Load data for case-crossover design
load("output/data_casecrossover.RData")
load("output/crossbasis_casecrossover.RData")

# Load data for time-series design
load("output/data_timeseries.RData")
load("output/crossbasis_timeseries.RData")
load("output/trend_timeseries.RData")
load("output/seasonality_timeseries.RData")

# Load spatial data structure and DLNM configuration
load("output/list_neighbours.RData")
load("output/dlnm_configuration.RData")

# For INLA implementation, add variables identifying different categories

# (1) Case-crossover
colnames(cb_cco) <- paste0("cb", 1:ncol(cb_cco))
data_cco <- cbind(data_cco, cb_cco) # Add cross-basis to the dataset
for (i in 1:dlnm_var$n_coef) {
  col_name <- paste0("id_cb", i)
  data_cco[[col_name]] <- as.numeric(data_cco$region) # One estimate per region for cb coefficients
}; rm(i, col_name)

# (2) Time-series
colnames(cb_ts) <- paste0("cb", 1:ncol(cb_ts)) # Add cross-basis to the dataset
data_ts <- cbind(data_ts, cb_ts)
for (i in 1:dlnm_var$n_coef) {
  col_name <- paste0("id_cb", i)
  data_ts[[col_name]] <- as.numeric(data_ts$region) # One estimate per region for cb coefficients
}; rm(i, col_name)

data_ts$trend <- trend[,1] # Add trend to the dataset
data_ts$id_trend <- as.numeric(data_ts$region)  # One estimate per region for trend coefficients

colnames(seas) <- paste0("seas", 1:dlnm_var$df_seas) 
data_ts <- cbind(data_ts, seas) # Add seasonality to the dataset
for (i in 1:dlnm_var$df_seas) {
  data_ts[[paste0("id_year", i)]] <- lubridate::year(data_ts$date) # One estimate per year and region for seasonality coefficients
  data_ts[[paste0("id_region", i)]] <- as.numeric(data_ts$region)
}; rm(i)

data_ts$intercept_dow <- paste0("dow_", data_ts$day_of_week, "_", data_ts$region) # Add day of week to the dataset

#...............................................................................
### MODEL 1 (Independent B-DLNM – case-crossover design) ####
#...............................................................................

### INLA - MODEL 1 Independent B-DLNM – case-crossover design

inla_formula <- mort ~ -1 + 
  f(strata, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb1, cb1, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb2, cb2, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb3, cb3, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb4, cb4, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb5, cb5, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb6, cb6, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb7, cb7, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb8, cb8, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb9, cb9, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb10, cb10, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb11, cb11, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb12, cb12, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE)))

# In order to sample from the posterior distribution model need to be fitted 
# with config = TRUE
inla_model <- inla(inla_formula,
                   data = data_cco, family = "poisson",
                   control.compute = list(config = TRUE), 
                   control.inla = list(strategy = "laplace",
                                       int.strategy = "grid"))

# Extract the ensemble of sample coefficients of the crossbasis
inla_res <- inla.posterior.sample(1000, inla_model,
                                  selection = list("id_cb1" = 1:dlnm_var$n_reg,
                                                   "id_cb2" = 1:dlnm_var$n_reg,
                                                   "id_cb3" = 1:dlnm_var$n_reg,
                                                   "id_cb4" = 1:dlnm_var$n_reg,
                                                   "id_cb5" = 1:dlnm_var$n_reg,
                                                   "id_cb6" = 1:dlnm_var$n_reg,
                                                   "id_cb7" = 1:dlnm_var$n_reg,
                                                   "id_cb8" = 1:dlnm_var$n_reg,
                                                   "id_cb9" = 1:dlnm_var$n_reg,
                                                   "id_cb10" = 1:dlnm_var$n_reg,
                                                   "id_cb11" = 1:dlnm_var$n_reg,
                                                   "id_cb12" = 1:dlnm_var$n_reg))

cb_res <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  beta_reg <- sapply(inla_res, function(x) {
    sapply(1:dlnm_var$n_coef, function(i) {
      x$latent[paste0("id_cb", i, ":", i_reg),]
    })
  })
  t(beta_reg)
})

# save(cb_res, file = "output/predicted_inla_model1_independent_casecrossover.RData")

rm(inla_formula, inla_model, inla_res, cb_res)

#...............................................................................
### MODEL 2 (Independent B-DLNM – time-series design) ####
#...............................................................................

### INLA - MODEL 2 Independent B-DLNM – time-series design

inla_formula <- mort ~ -1 + 
  f(intercept_dow, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) +
  f(id_cb1, cb1, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb2, cb2, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb3, cb3, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb4, cb4, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb5, cb5, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb6, cb6, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb7, cb7, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb8, cb8, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb9, cb9, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb10, cb10, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb11, cb11, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb12, cb12, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_trend, trend, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) +
  seas1:factor(id_year1):factor(id_region1) + 
  seas2:factor(id_year2):factor(id_region2) +
  seas3:factor(id_year3):factor(id_region3) + 
  seas4:factor(id_year4):factor(id_region4)

# In order to sample from the posterior distribution model need to be fitted 
# with config = TRUE
inla_model <- inla(inla_formula,
                   data = data_ts, family = "poisson",
                   control.compute = list(config = TRUE),
                   control.inla = list(strategy = "laplace",
                                       int.strategy = "grid"))

# Extract the ensemble of sample coefficients of the crossbasis
inla_res <- inla.posterior.sample(1000, inla_model,
                                  selection = list("id_cb1" = 1:dlnm_var$n_reg,
                                                   "id_cb2" = 1:dlnm_var$n_reg,
                                                   "id_cb3" = 1:dlnm_var$n_reg,
                                                   "id_cb4" = 1:dlnm_var$n_reg,
                                                   "id_cb5" = 1:dlnm_var$n_reg,
                                                   "id_cb6" = 1:dlnm_var$n_reg,
                                                   "id_cb7" = 1:dlnm_var$n_reg,
                                                   "id_cb8" = 1:dlnm_var$n_reg,
                                                   "id_cb9" = 1:dlnm_var$n_reg,
                                                   "id_cb10" = 1:dlnm_var$n_reg,
                                                   "id_cb11" = 1:dlnm_var$n_reg,
                                                   "id_cb12" = 1:dlnm_var$n_reg))

cb_res <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  beta_reg <- sapply(inla_res, function(x) {
    sapply(1:dlnm_var$n_coef, function(i) {
      x$latent[paste0("id_cb", i, ":", i_reg),]
    })
  })
  t(beta_reg)
})

# save(cb_res, file = "output/predicted_inla_model2_independent_timeseries.RData")

rm(inla_formula, inla_model, inla_res, cb_res)

#...............................................................................
### MODEL 3 (SB-DLNM – case-crossover design) ####
#...............................................................................

### INLA - MODEL 3 SB-DLNM – case-crossover design

inla_formula <- mort ~ -1 + 
  cb1 + cb2 + cb3 + cb4 + cb5 + cb6 + 
  cb7 + cb8 + cb9 + cb10 + cb11 + cb12 +
  f(strata, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
  f(id_cb1, cb1, model = "bym2", graph = list_neig) + 
  f(id_cb2, cb2, model = "bym2", graph = list_neig) +
  f(id_cb3, cb3, model = "bym2", graph = list_neig) + 
  f(id_cb4, cb4, model = "bym2", graph = list_neig) +
  f(id_cb5, cb5, model = "bym2", graph = list_neig) + 
  f(id_cb6, cb6, model = "bym2", graph = list_neig) +
  f(id_cb7, cb7, model = "bym2", graph = list_neig) +
  f(id_cb8, cb8, model = "bym2", graph = list_neig) +
  f(id_cb9, cb9, model = "bym2", graph = list_neig) +
  f(id_cb10, cb10, model = "bym2", graph = list_neig) + 
  f(id_cb11, cb11, model = "bym2", graph = list_neig) +
  f(id_cb12, cb12, model = "bym2", graph = list_neig)

# In order to sample from the posterior distribution model need to be fitted 
# with config = TRUE
inla_model <- inla(inla_formula,
                   data = data_cco, family = "poisson",
                   control.compute = list(config = TRUE),
                   control.inla = list(strategy = "laplace",
                                       int.strategy = "grid"))

# Extract the ensemble of sample coefficients of the crossbasis
inla_res <- inla.posterior.sample(1000, inla_model,
                                  selection = list(cb1 = 1,
                                                   cb2 = 1,
                                                   cb3 = 1,
                                                   cb4 = 1,
                                                   cb5 = 1,
                                                   cb6 = 1,
                                                   cb7 = 1,
                                                   cb8 = 1,
                                                   cb9 = 1,
                                                   cb10 = 1,
                                                   cb11 = 1,
                                                   cb12 = 1,
                                                   "id_cb1" = 1:dlnm_var$n_reg,
                                                   "id_cb2" = 1:dlnm_var$n_reg,
                                                   "id_cb3" = 1:dlnm_var$n_reg,
                                                   "id_cb4" = 1:dlnm_var$n_reg,
                                                   "id_cb5" = 1:dlnm_var$n_reg,
                                                   "id_cb6" = 1:dlnm_var$n_reg,
                                                   "id_cb7" = 1:dlnm_var$n_reg,
                                                   "id_cb8" = 1:dlnm_var$n_reg,
                                                   "id_cb9" = 1:dlnm_var$n_reg,
                                                   "id_cb10" = 1:dlnm_var$n_reg,
                                                   "id_cb11" = 1:dlnm_var$n_reg,
                                                   "id_cb12" = 1:dlnm_var$n_reg))

cb_res <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  beta_reg <- sapply(inla_res, function(x) {
    sapply(1:dlnm_var$n_coef, function(i) {
      x$latent[paste0("cb", i, ":1"),] + 
        x$latent[paste0("id_cb", i, ":", i_reg),]
    })
  })
  t(beta_reg)
})

# save(cb_res, file = "output/predicted_inla_model3_spatial_casecrossover.RData")

rm(sdunif, inla_formula, inla_model, inla_res, cb_res)

#...............................................................................
### MODEL 4 (SB-DLNM – time-series design) ####
#...............................................................................

### INLA - MODEL 4 SB-DLNM – time-series design

inla_formula <- mort ~ -1 + 
  cb1 + cb2 + cb3 + cb4 + cb5 + cb6 + 
  cb7 + cb8 + cb9 + cb10 + cb11 + cb12 +
  f(intercept_dow, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) +
  f(id_cb1, cb1, model = "bym2", graph = list_neig) + 
  f(id_cb2, cb2, model = "bym2", graph = list_neig) + 
  f(id_cb3, cb3, model = "bym2", graph = list_neig) + 
  f(id_cb4, cb4, model = "bym2", graph = list_neig) + 
  f(id_cb5, cb5, model = "bym2", graph = list_neig) + 
  f(id_cb6, cb6, model = "bym2", graph = list_neig) + 
  f(id_cb7, cb7, model = "bym2", graph = list_neig) + 
  f(id_cb8, cb8, model = "bym2", graph = list_neig) + 
  f(id_cb9, cb9, model = "bym2", graph = list_neig) + 
  f(id_cb10, cb10, model = "bym2", graph = list_neig) + 
  f(id_cb11, cb11, model = "bym2", graph = list_neig) + 
  f(id_cb12, cb12, model = "bym2", graph = list_neig) + 
  f(id_trend, trend, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) +
  seas1:factor(id_year1):factor(id_region1) + 
  seas2:factor(id_year2):factor(id_region2) +
  seas3:factor(id_year3):factor(id_region3) + 
  seas4:factor(id_year4):factor(id_region4)

# In order to sample from the posterior distribution model need to be fitted 
# with config = TRUE
inla_model <- inla(inla_formula,
                   data = data_ts, family = "poisson",
                   control.compute = list(config = TRUE),
                   control.inla = list(strategy = "laplace",
                                       int.strategy = "grid"))

# Extract the ensemble of sample coefficients of the crossbasis
inla_res <- inla.posterior.sample(1000, inla_model,
                                  selection = list(cb1 = 1,
                                                   cb2 = 1,
                                                   cb3 = 1,
                                                   cb4 = 1,
                                                   cb5 = 1,
                                                   cb6 = 1,
                                                   cb7 = 1,
                                                   cb8 = 1,
                                                   cb9 = 1,
                                                   cb10 = 1,
                                                   cb11 = 1,
                                                   cb12 = 1,
                                                   "id_cb1" = 1:dlnm_var$n_reg,
                                                   "id_cb2" = 1:dlnm_var$n_reg,
                                                   "id_cb3" = 1:dlnm_var$n_reg,
                                                   "id_cb4" = 1:dlnm_var$n_reg,
                                                   "id_cb5" = 1:dlnm_var$n_reg,
                                                   "id_cb6" = 1:dlnm_var$n_reg,
                                                   "id_cb7" = 1:dlnm_var$n_reg,
                                                   "id_cb8" = 1:dlnm_var$n_reg,
                                                   "id_cb9" = 1:dlnm_var$n_reg,
                                                   "id_cb10" = 1:dlnm_var$n_reg,
                                                   "id_cb11" = 1:dlnm_var$n_reg,
                                                   "id_cb12" = 1:dlnm_var$n_reg))

cb_res <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  beta_reg <- sapply(inla_res, function(x) {
    sapply(1:dlnm_var$n_coef, function(i) {
      x$latent[paste0("cb", i, ":1"),] + 
        x$latent[paste0("id_cb", i, ":", i_reg),]
    })
  })
  t(beta_reg)
})

# save(cb_res, file = "output/predicted_inla_model4_spatial_timeseries.RData")

rm(sdunif, inla_formula, inla_model, inla_res, cb_res)