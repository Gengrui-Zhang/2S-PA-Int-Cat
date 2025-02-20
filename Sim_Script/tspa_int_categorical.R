library(semTools)
library(lavaan)
library(SimDesign)
library(mnormt)
library(dplyr)
library(tidyverse)
library(semlrtp)
library(MplusAutomation)
library(here)

# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# ========================================= Simulation Conditions ========================================= #

DESIGNFACTOR <- createDesign(
  N = c(10000, 25000, 50000),
  num_cat = c(3, 5), # response categories
  rel = c(0.7, 0.8, 0.9),
  gamma_xm = c(0, 0.3), # Two levels of the interaction effect
  skewness = c("symm", "skew")
)

FIXED <- list(gamma_x = 0.3,
              gamma_m = 0.3,
              cor_xm = 0.3,
              num_x_ind = 3,
              num_m_ind = 12,
              y_int = 1.2,
              # Measurement parameters
              lambda_y = c(.5, .7, .9) # loadings for y indicators
              )

# Temporary directory
FIXED$temp_session <- "/Users/jimmy_z/R Projects/2S-PA-Int-Cat/temp_folder"

# LMS syntax
FIXED$lms_syntax <- "
TITLE: Latent Moderated Structural (LMS) Model

DATA:
  FILE = sim_dat_repid.dat;  
  FORMAT IS FREE;              

VARIABLE:
  NAMES = x1 x2 x3 
  	    m1 m2 m3 m4 m5 m6 
  	    m7 m8 m9 m10 m11 m12 
  	    y1 y2 y3;  
  USEVARIABLES = x1 x2 x3 
  		    m1 m2 m3 m4 m5 m6 
  		    m7 m8 m9 m10 m11 m12 
                 y1 y2 y3;  
  MISSING = ALL (-999);  
  CATEGORICAL = x1 x2 x3 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12; 

ANALYSIS:
  TYPE = RANDOM;    
  ESTIMATOR = MLR;  
  ALGORITHM = INTEGRATION;  

MODEL:
  X BY x1* x2* x3*;   
  M BY m1* m2* m3*;   
  Y BY y1* y2* y3*;   
  XM | X XWITH M;  

  Y ON X* M* XM*;   

  X@1;  ! Fix variance of X
  M@1;  ! Fix variance of M
  X WITH M*;
  Y*;

OUTPUT:
  sampstat standardized tech1;  
  
SAVEDATA:
  RESULTS = lms_results_repid.dat;
"

# ========================================= Data Generation ========================================= #
# Check function
# condition = DESIGNFACTOR[1,]
# fixed_objects = FIXED

# Helper Function
generate_dat <- function(condition, fixed_objects = NULL) {
  N <- condition$N
  num_cats <- condition$num_cat
  skewness <- condition$skewness
  gamma_xm <- condition$gamma_xm
  gamma_x <- fixed_objects$gamma_x
  gamma_m <- fixed_objects$gamma_m
  cor_xm <- fixed_objects$cor_xm
  num_eta1_ind <- fixed_objects$num_x_ind
  num_eta2_ind <- fixed_objects$num_m_ind
  y_int <- fixed_objects$y_int
  y_dist <- 1 - (gamma_x^2 + gamma_m^2 + gamma_xm^2 + 2*gamma_x*gamma_m*cor_xm)
  lambda_y <- fixed_objects$lambda_y
  # Threshold
  if (condition$skewness == "symm") {
    thres_x <- matrix(rep(0, 3), nrow = 1)  # binary
    col_vector <- c(-1.5, -0.5, 0.5, 1.5)
    thres_m <- matrix(rep(col_vector, times = 12), 
                     nrow = 4, ncol = 12, byrow = FALSE)
  } else {
    thres_x <- matrix(rep(0.9, 3), nrow = 1)  # binary
    col_vector <- c(0.05, 0.75, 1.55, 2.55)
    thres_m <- matrix(rep(col_vector, times = 12), 
                      nrow = 4, ncol = 12, byrow = FALSE)
  }
  # Reliability
  rel <- condition$rel
  
  # Generate latent predictors and interaction term
  cov_xm_ey <- matrix(c(1, cor_xm, 0,
                        cor_xm, 1, 0,
                        0, 0, y_dist), nrow = 3)
  eta <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = cov_xm_ey,
                       empirical = TRUE) # x and m are standardized normal variables
  eta <- cbind(eta, eta[, 1] * eta[, 2]) # Add latent interaction term
  etay <- y_int + eta[, -3] %*% c(gamma_x, gamma_m, gamma_xm) + eta[, 3] # Generate latent outcome variable
  
  # Check if the latent variables are simulated correctly 
  # lm(etay ~ eta[, -3])
  
  # Generate outcome variables: y1, y2, y3
  Y <- t(lambda_y %*% t(etay) + matrix(rnorm(N * length(lambda_y), 
                                           sd = rep(sqrt(1 - lambda_y^2), each = N)),
                                     nrow = length(lambda_y), ncol = N))
  
  # Generate measurement indicators (latent)
  lambda_x = 1.7*seq(0.6, 0.8, length.out = num_eta1_ind)  # on normal ogive metric
  lambda_m = 1.7*seq(0.3, 0.85, length.out = num_eta2_ind)  # on normal ogive metric
  err_var_x <- sum(lambda_x)^2*(1 - rel)/rel*seq(0.8, 0.2, length.out = num_eta1_ind)
  err_var_x <- err_var_x/sum(err_var_x)
  err_var_m <- sum(lambda_m)^2*(1 - rel)/rel*seq(0.8, 0.2, length.out = num_eta2_ind)
  err_var_m <- err_var_m/sum(err_var_m)
  x_lat <- eta[, 1] %*% t(lambda_x) + matrix(rnorm(N * length(lambda_x), 
                                                   sd = rep(sqrt(err_var_x), each = N)),
                                             nrow = N, ncol = length(lambda_x))
  m_lat <- eta[, 2] %*% t(lambda_m) + matrix(rnorm(N * length(lambda_m), 
                                                   sd = rep(sqrt(err_var_m), each = N)),
                                             nrow = N, ncol = length(lambda_m))
  
  # Check if simulated continuous items are correctly parameterized
  # test_df <- cbind(x_lat, m_lat, Y)
  # colnames(test_df) <- c("x1", "x2", "x3",
  #                        "m1", "m2", "m3", "m4", "m5", "m6",
  #                        "m7", "m8", "m9", "m10", "m11", "m12",
  #                        "y1", "y2", "y3"
  #                        )
  # mod <- "X =~ x1 + x2 + x3
  #         M =~ m1 + m2 + m3 + m4 + m5 + m6 +
  #             m7 + m8 + m9 + m10 + m11 + m12
  #         Y =~ y1 + y2 + y3
  #         Y ~ X + M"
  # summary(sem(mod, test_df, std.lv = TRUE),
  #         standardized = TRUE)
  
  # Generate measurement indicators (observed)
  #' Obtain categorical items
  X <- vapply(
    seq_along(lambda_x),
    FUN = \(i) {
      findInterval(x_lat[, i], thres_x[, i])
    },
    FUN.VALUE = numeric(N))
  M <- vapply(
    seq_along(lambda_m),
    FUN = \(i) {
      findInterval(m_lat[, i], thres_m[, i])
    },
    FUN.VALUE = numeric(N))
  
  # Rename indicators
  dat <- data.frame(X, M, Y)
  colnames(dat) <- c("x1", "x2", "x3",
                     "m1", "m2", "m3", "m4", "m5", "m6",
                     "m7", "m8", "m9", "m10", "m11", "m12",
                     "y1", "y2", "y3")
  
  # Return data
  return(dat)
}

# Test df
# test_dat <- generate_dat(DESIGNFACTOR[1,], FIXED)
# hist(test_dat$x1)
# hist(test_dat$m1)

# ========================================= Data Analysis ========================================= #

analyze_2spa <- function (condition, dat, fixed_objects = NULL) {
  
  # Prepare factor scores
  fs_y <- get_fs(dat[c(paste0("y", 1:3))], std.lv = TRUE) # continuous indicators
  irt_x <- mirt(dat[c(paste0("x", 1:3))], 
                itemtype = "2PL",
                verbose = FALSE)  # IRT (2-PL) for binary items
  fs_x <- fscores(irt_x, full.scores.SE = TRUE) # EAP scores
  irt_m <- mirt(dat[c(paste0("m", 1:12))], 
                itemtype = "graded",
                verbose = FALSE)  # GRM for items with multiple categories
  fs_m <- fscores(irt_m, full.scores.SE = TRUE) # EAP scores
  
  # Prepare fs_dat
  fs_dat <- data.frame(fs_y[c(1, 3, 4)], fs_x, fs_m) |>
    setNames(c(
      "fs_y", "rel_fs_y", "ev_fs_y",
      "fs_x", "se_fs_x", "fs_m", "se_fs_m"
    )) |>
    # Compute reliability; only needed for 2S-PA
    within(expr = {
      rel_fs_x <- 1 - se_fs_x^2 # reliability
      rel_fs_m <- 1 - se_fs_m^2
      ev_fs_x <- se_fs_x^2 * (1 - se_fs_x^2) # error variance (reliability with incorporated uncertainty)
      ev_fs_m <- se_fs_m^2 * (1 - se_fs_m^2)
    }) |>
    # Add interaction
    within(expr = {
      fs_xm <- fs_x * fs_m
      ld_fs_xm <- rel_fs_x * rel_fs_m
      ev_fs_xm <- rel_fs_x^2 * ev_fs_m + rel_fs_m^2 * ev_fs_x + ev_fs_m * ev_fs_x # formula
    })
  
  # OpenMX model
  tspa_umx <- umxLav2RAM(
    "
      fs_y ~ b1 * fs_x + b2 * fs_m + b3 * fs_xm
      fs_y + fs_x + fs_m + fs_xm ~ 1
      fs_x ~~ vx * fs_x + fs_m + fs_xm
      fs_m ~~ vm * fs_m + fs_xm
    ",
    printTab = FALSE)
    # Loading
    matL <- mxMatrix(
      type = "Diag", nrow = 4, ncol = 4,
      free = FALSE,
      labels = c("data.rel_fs_y", "data.rel_fs_x", "data.rel_fs_m", "data.ld_fs_xm"),
      name = "L"
    )
    # Error
    matE <- mxMatrix(
      type = "Diag", nrow = 4, ncol = 4,
      free = FALSE,
      labels = c("data.ev_fs_y", "data.ev_fs_x", "data.ev_fs_m", "data.ev_fs_xm"),
      name = "E"
    )
    
  # Run the OpenMX model
  tspa_mx <-
      tspa_mx_model(
        mxModel(tspa_umx,
                mxAlgebra(b1 * sqrt(vx), name = "stdx_b1"),
                mxAlgebra(b2 * sqrt(vm), name = "stdx_b2"),
                mxAlgebra(b3 * sqrt(vx) * sqrt(vm), name = "stdx_b3"),
                mxCI(c("stdx_b1", "stdx_b2", "stdx_b3"))),
        data = fs_dat,
        mat_ld = matL, 
        mat_ev = matE)
  tspa_mx_fit <- mxRun(tspa_mx, intervals = TRUE)
}

analyze_lms <- function (condition, dat, fixed_objects = NULL) {
  temp_session <- fixed_objects$temp_session
  dat_name <- file.path(temp_session,
                        sprintf("sim_dat_%s.dat", condition$REPLICATION))
  write.table(dat,
              file = dat_name,
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  mplus_files <- sprintf(
    c("sim_dat_%s.dat", "lms_%s.inp", "lms_%s.out"),
    condition$REPLICATION
  )
  writeLines(gsub("repid", replacement = condition$REPLICATION,
                  x = fixed_objects$lms_syntax),
             con = file.path(temp_session, mplus_files[[2]]))
  MplusAutomation::runModels(
    temp_session,
    filefilter = mplus_files[[2]],
    Mplus_command = "/Applications/Mplus/mplus")
  
  # Extract results
  res_ust <- readModels(file.path(temp_session, mplus_files[[3]]), what="parameters")$parameters$unstandardized
  res_std <- readModels(file.path(temp_session, mplus_files[[3]]), what="parameters")$parameters$std.standardized
  
  est_ust <- res_ust[res_ust$paramHeader == "Y.ON", "est"]
  se_ust <- res_ust[res_ust$paramHeader == "Y.ON", "se"]
  est_std <- res_std[res_std$paramHeader == "Y.ON", "est"]
  se_std <- res_std[res_std$paramHeader == "Y.ON", "se"]
  p <- res_std[res_std$paramHeader == "Y.ON", "pval"]
  
  # Create the output vector
  out <- c(est_ust, se_ust, est_std, se_std, p)  
  names(out) <- c("est_x_ust", "est_m_ust", "est_xm_ust",
                  "se_x_ust", "se_m_ust", "se_xm_ust", 
                  "est_x_std", "est_m_std", "est_xm_std",
                  "se_x_std", "se_m_std", "se_xm_std", 
                  "pval_x", "pval_m", "pval_xm")
  
  # Remove temp files
  file.remove(file.path(temp_session, mplus_files), recursive = TRUE)
  
  # Return
  return(out)
}

# ========================================= Results Summary ========================================= #
# Helper function: robust bias
robust_bias <- function(est, se, par, trim = 0, type = NULL) {
  output <- numeric(ncol(est))
  for (i in seq_len(ncol(est))) {
    if (type == "raw") {
      output[i] <- mean((est[,i] - par), na.rm = TRUE)
    } else if (type == "standardized") {
      output[i] <- (mean(est[,i], na.rm = TRUE) - par)/sd(est[,i], na.rm = TRUE)
    } else if (type == "trim") {
      output[i] <- mean(est[,i], trim = trim, na.rm = TRUE) - par
    } else if (type == "median") {
      output[i] <- (median(est[,i], na.rm = TRUE) - par) / mad(est[,i], na.rm = TRUE)
    } else {
      output[i] <- (mean(est[,i], trim = trim, na.rm = TRUE) - par) / sd(est[,i], na.rm = TRUE)
    }
  }
  names(output) <- colnames(est)
  return(output)
}

# Helper function: relative SE bias
rse_bias <- function(est, se, trim = 0, type = "raw") {
  if (type == "raw") {
    se_mean <- apply(se, 2, mean, na.rm = T)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  } else if (type == "median") {
    se_median <- apply(se, 2, median, na.rm = TRUE)
    mad_sd <- apply(est, 2, function(x) mad(x, na.rm = TRUE))
    rse_bias <- se_median / mad_sd - 1
  } else if (type == "trim") {
    se_mean <- apply(se, 2, mean, trim = trim, na.rm = TRUE)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  }
  return(rse_bias)
}

# Helper function: detecting outliers for SE
outlier_se <- function(se) {
  results <- c()
  for(column in colnames(se)) {
    # Calculate Q1, Q3, and IQR
    Q1 <- quantile(se[,column], 0.25, na.rm = TRUE)
    Q3 <- quantile(se[,column], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    # Determine outliers
    lower_bound <- (Q1 - 1.5 * IQR)
    upper_bound <- (Q3 + 1.5 * IQR)
    outliers <- se[,column][se[,column] < lower_bound | se[,column] > upper_bound]
    # Calculate the percentage of outliers
    percentage <- length(outliers) / sum(!is.na(se[,column])) * 100
    results[column] <- percentage
  }
  return(results)
}

# Helper function for calculating coverage rate, Type I error rate, and power
ci_stats <- function(est, se, par, stats_type, lrt_lo = NULL, lrt_up = NULL) {

  # Calculate the confidence intervals (usd)
  lo_95 <- est - qnorm(.975) * se
  up_95 <- est + qnorm(.975) * se
  ci_est <- vector("list", length = ncol(est))
  names(ci_est) <- colnames(est)
  
  # Construct confidence intervals for each method
  for (i in seq_len(ncol(est))) {
    ci_est[[i]] <- cbind(lo_95[,i], up_95[,i])
  }
  
  # Extract LRT CIs
  if (!is.null(lrt_lo) && !is.null(lrt_up)) {
    ci_lrt <- vector("list", length = ncol(est))
    names(ci_lrt) <- colnames(est)
    
    for (i in seq_len(ncol(est))) {
      ci_lrt[[i]] <- cbind(lrt_lo[,i], lrt_up[,i])
    }
  }
  
  # Determine which statistic to calculate
  if (stats_type == "Coverage") {
    return(sapply(ci_est, function(ci) mean(ci[,1] <= par & ci[,2] >= par)))
  } else if (stats_type == "TypeI") {
    return(sapply(ci_est, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Lrt_TypeI") {
    return(sapply(ci_lrt, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Power") {
    return(sapply(ci_est, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
  } else if (stats_type == "Lrt_Power") {
    return(sapply(ci_lrt, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
  } else {
    stop("Invalid stats_type specified. Please choose from 'Coverage', 'TypeI', or 'Power'.")
  }
}

# Helper function for warning sum
warning_sum <- function(count) {
  apply(count, 2, sum, na.rm = TRUE)
}

# Evaluation Function
evaluate_res <- function (condition, results, fixed_objects = NULL) {

  # Population parameter
  pop_par <- condition$gamma_xm

  # Parameter estimates
  est_std <- results[, grep(".est$", colnames(results))]
  est_usd <- results[, grep(".est_usd$", colnames(results))]
  se_std <- results[, grep(".se_std$", colnames(results))]
  se_usd <- results[, grep(".se_usd$", colnames(results))]
  lrtp_lower <- results[, grep(".lrtp_lower$", colnames(results))]
  lrtp_upper <- results[, grep(".lrtp_upper$", colnames(results))]
  warnings <- results[, grep(".warnings_count$", colnames(results))]

  c(raw_bias = robust_bias(est_std,
                           se_std,
                           pop_par,
                           type = "raw"),
    std_bias = robust_bias(est_std,
                           se_std,
                           pop_par,
                           type = "standardized"),
    trim_bias = robust_bias(est_std,
                            se_std,
                            pop_par,
                            trim = 0.2,
                            type = "trim"), # 20% trimmed mean
    stdMed_bias = robust_bias(est_std,
                              se_std,
                              pop_par,
                              type = "median"),
    raw_rse_bias = rse_bias(est_std,
                            se_std,
                            type = "raw"),
    stdMed_rse_bias = rse_bias(est_std,
                               se_std,
                               type = "median"),
    trim_rse_bias = rse_bias(est_std,
                             se_std,
                             trim = 0.2,
                             type = "trim"),
    outlier_se = outlier_se(se_std),
    coverage_usd = ci_stats(est_usd, 
                        se_usd,
                        pop_par, 
                        "Coverage"),
    coverage_std = ci_stats(est_std, 
                            se_std,
                            pop_par, 
                            "Coverage"),
    type1_usd = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "TypeI"),
    type1_std = ci_stats(est_std, 
                         se_std,
                         pop_par, 
                         "TypeI"),
    type1_lrt = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Lrt_TypeI",
                         lrt_lo = lrtp_lower,
                         lrt_up = lrtp_upper),
    power_usd = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Power"),
    power_std = ci_stats(est_std, 
                         se_std,
                         pop_par, 
                         "Power"),
    power_lrt = ci_stats(est_usd, 
                         se_usd,
                         pop_par, 
                         "Lrt_Power",
                         lrt_lo = lrtp_lower,
                         lrt_up = lrtp_upper),
    rmse = RMSE(na.omit(est_std),
                parameter = pop_par),
    warning_total = warning_sum(warnings)
  )
}

# ========================================= Run Experiment ========================================= #
res <- runSimulation(design = DESIGNFACTOR,
              replications = 2000,
              generate = generate_dat,
              analyse = list(mmr = analyze_mmr,
                             upi = analyze_upi,
                             rapi = analyze_rapi,
                             tspa = analyze_tspa),
              summarise = evaluate_res,
              fixed_objects = FIXED_PARAMETER,
              seed = rep(61543, nrow(DESIGNFACTOR)),
              packages = "lavaan", 
              filename = "continuous_boundMLR_09252024",
              parallel = TRUE,
              ncores = 30,
              save = TRUE,
              save_results = TRUE)
