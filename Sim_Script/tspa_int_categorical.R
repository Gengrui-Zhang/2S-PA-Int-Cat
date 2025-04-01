library(semTools)
library(lavaan)
library(SimDesign)
library(mnormt)
library(dplyr)
library(tidyverse)
library(semlrtp)
library(mirt)
library(umx)
library(OpenMx)
library(R2spa)
library(MplusAutomation)
library(here)

# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# ========================================= Simulation Conditions ========================================= #

DESIGNFACTOR <- createDesign(
  N = c(100, 250, 500),
  # num_cat = c(3, 5),
  rel = 0.8,
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
FIXED$temp_session <- tempdir()

# UPI syntax
FIXED$upi_syntax <- "
  # Measurement Model
    Y =~ y1 + y2 + y3
    X =~ x1 + x2 + x3
    M =~ m1 + m2 + m3 + m4 + m5 + m6 + 
         m7 + m8 + m9 + m10 + m11 + m12 
  # Structural Model
    Y ~ b1*X + b2*M + b3*X:M
  # Define Standardized Coefficients
    X ~~ v1*X
    M ~~ v2*M
  # Standardized coefficients
    stdx_b1 := b1*sqrt(v1)
    stdx_b2 := b2*sqrt(v2)
    stdx_b3 := b3*sqrt(v1)*sqrt(v2)
"

# LMS syntax
FIXED$lms_syntax <- "
TITLE: Latent Moderated Structural (LMS) Model

DATA:
  FILE = lms_simdat_repid.dat;  
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
  M BY m1* m2* m3* m4* m5* m6* m7* m8* m9* m10* m11* m12*;   
  Y BY y1* y2* y3*;   
  XM | X XWITH M;  

  Y ON X* M* XM*;   

  X@1;  ! Fix variance of X
  M@1;  ! Fix variance of M
  X WITH M*;
  Y*;

OUTPUT:
  sampstat standardized tech1 CINTERVAL;  
  
SAVEDATA:
  RESULTS = lms_results_repid.dat;
"

# 2S-PA Syntax
FIXED$tspa_syntax <- "
TITLE: Two-Stage Path Analysis

DATA:
  FILE = 2spa_simdat_repid.dat;  
  FORMAT IS FREE;              

VARIABLE:
  NAMES = fs_y rel_fs_y ev_fs_y
          fs_x se_fs_x fs_m se_fs_m
          ev_fs_m ev_fs_x rel_fs_m rel_fs_x 
          ev_fs_xm ld_fs_xm fs_xm;
  USEVARIABLES = fs_y fs_x fs_m fs_xm;
  CONSTRAINT = rel_fs_y rel_fs_m rel_fs_x ld_fs_xm
               ev_fs_y ev_fs_m ev_fs_x ev_fs_xm;
        
ANALYSIS:
        
MODEL:  !outcome
          Y BY fs_y*1 (ld_fsy);
          fs_y* (ev_fsy);
          
          !predictor
          X BY fs_x*1 (ld_fsx);
          X (var_X);
          fs_x* (ev_fsx);

          !moderator
          M BY fs_m*1 (ld_fsm);
          M (var_M);
          fs_m* (ev_fsm);

          !product indicator
          XM BY fs_xm*1 (ld_fsxm);
          fs_xm* (ev_fsxm);

          !structural model
          Y ON X(b1); 
          Y ON M(b2);
          Y ON XM(b3);
          
  MODEL CONSTRAINT:
          ld_fsy = rel_fs_y;
          ev_fsy = ev_fs_y;  
          ld_fsx = rel_fs_x;
          ev_fsx = ev_fs_x;
          ld_fsm = rel_fs_m;
          ev_fsm = ev_fs_m;
          ld_fsxm = ld_fs_xm;
          ev_fsxm = ev_fs_xm;
          
          New(stdx_b1);
          stdx_b1 = b1*sqrt(var_X);
          New(stdx_b2);
          stdx_b2 = b2*sqrt(var_M);
          New(stdx_b3);
          stdx_b3 = b3*sqrt(var_X)*sqrt(var_M);

OUTPUT:
  sampstat tech1 CINTERVAL;  
  
SAVEDATA:
  RESULTS = 2spa_results_repid.dat;
"

# LMS with factor scores
FIXED$lmsfs_syntax <- "
TITLE: LMS using factor scores

DATA:
  FILE = lmsfs_simdat_repid.dat;
  FORMAT IS FREE;              

VARIABLE:
  NAMES = fs_y rel_fs_y ev_fs_y
          fs_x se_fs_x fs_m se_fs_m
          ev_fs_m ev_fs_x rel_fs_m rel_fs_x 
          ev_fs_xm ld_fs_xm fs_xm; 
  USEVARIABLES = fs_y fs_x fs_m;  
  CONSTRAINT = rel_fs_y rel_fs_m rel_fs_x 
               ev_fs_y ev_fs_m ev_fs_x;

ANALYSIS:
        TYPE = RANDOM;
        ALGORITHM = INTEGRATION;
MODEL:  !outcome
          Y BY fs_y*1 (ld_fsy);
          fs_y* (ev_fsy);
          
          !predictor
          X BY fs_x*1 (ld_fsx);
          X (var_X);
          fs_x* (ev_fsx);

          !moderator
          M BY fs_m*1 (ld_fsm);
          M (var_M);
          fs_m* (ev_fsm);

          !interaction model
          XM | X XWITH M;
          
          !structural model
          Y ON X(b1); 
          Y ON M(b2);
          Y ON XM(b3);
          
  MODEL CONSTRAINT:
          ld_fsy = rel_fs_y;
          ev_fsy = ev_fs_y; 
          ld_fsx = rel_fs_x;
          ev_fsx = ev_fs_x;
          ld_fsm = rel_fs_m;
          ev_fsm = ev_fs_m;
          
          New(stdx_b1);
          stdx_b1 = b1*sqrt(var_X);
          New(stdx_b2);
          stdx_b2 = b2*sqrt(var_M);
          New(stdx_b3);
          stdx_b3 = b3*sqrt(var_X)*sqrt(var_M);

OUTPUT:
  sampstat tech1 CINTERVAL;  
  
SAVEDATA:
  RESULTS = lmsfs_results_repid.dat;
"

# ========================================= Data Generation ========================================= #
# Check function
# condition = DESIGNFACTOR[1,]
# fixed_objects = FIXED

# Helper Function
generate_dat <- function(condition, fixed_objects = NULL) {
  N <- condition$N
  # num_cats <- condition$num_cat
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
    col_vector <- c(-1.5, -0.5, 0.5, 1.5) # multiple categories 
    thres_m <- matrix(rep(col_vector, times = 12), 
                     nrow = 4, ncol = 12, byrow = FALSE)
  } else {
    thres_x <- matrix(rep(0.9, 3), nrow = 1)  # binary
    col_vector <- c(0.05, 0.75, 1.55, 2.55) # multiple categories 
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
  # Set error variance of continuous indicators to 1 
  # err_var_x <- sum(lambda_x)^2*(1 - rel)/rel*seq(0.8, 0.2, length.out = num_eta1_ind)
  # err_var_x <- err_var_x/sum(err_var_x)
  # err_var_m <- sum(lambda_m)^2*(1 - rel)/rel*seq(0.8, 0.2, length.out = num_eta2_ind)
  # err_var_m <- err_var_m/sum(err_var_m)
  x_lat <- eta[, 1] %*% t(lambda_x) + matrix(rnorm(N * length(lambda_x), 
                                                   sd = 1), # sd = 1
                                             nrow = N, ncol = length(lambda_x))
  m_lat <- eta[, 2] %*% t(lambda_m) + matrix(rnorm(N * length(lambda_m), 
                                                   sd = 1),
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

generate_fsdat <- function (dat) {
  # Prepare factor scores
  fs_y <- get_fs(dat[c(paste0("y", 1:3))], std.lv = TRUE) # continuous indicators
  
  irt_x <- mirt(dat[c(paste0("x", 1:3))], # IRT (2-PL) for binary items
                itemtype = "2PL",
                verbose = FALSE,
                technical = list(NCYCLES = 1e4)) # Increase iteration number
  if (!irt_x@OptimInfo$converged) stop('No convergence for MIRT') # Check for convergence
  fs_x <- fscores(irt_x, full.scores.SE = TRUE) # EAP scores
  
  irt_m <- mirt(dat[c(paste0("m", 1:12))], # GRM for items with multiple categories
                itemtype = "graded",
                verbose = FALSE,
                technical = list(NCYCLES = 1e4)) # # Increase iteration number
  if (!irt_m@OptimInfo$converged) stop('No convergence for MIRT') # Check for convergence
  fs_m <- fscores(irt_m, full.scores.SE = TRUE) # EAP scores
  
  # Prepare fs_dat
  fs_dat <- data.frame(fs_y[c(1, 3, 4)], fs_x, fs_m) |>
    setNames(c(
      "fs_y", "rel_fs_y", "ev_fs_y",
      "fs_x", "se_fs_x", "fs_m", "se_fs_m"
    )) |>
    mutate(
      fs_x = fs_x - mean(fs_x, na.rm = TRUE),  # Mean-centered fs_x
      fs_m = fs_m - mean(fs_m, na.rm = TRUE),  # Mean-centered fs_m
    ) |>
    within(expr = {
      rel_fs_x <- 1 - se_fs_x^2 # reliability
      rel_fs_m <- 1 - se_fs_m^2
      ev_fs_x <- se_fs_x^2 * (1 - se_fs_x^2) # error variance (reliability with incorporated uncertainty)
      ev_fs_m <- se_fs_m^2 * (1 - se_fs_m^2)
    }) |>
    # Add interaction
    within(expr = {
      # Double Mean-centered interaction
      fs_xm <- fs_x * fs_m - mean(fs_x * fs_m)
      ld_fs_xm <- rel_fs_x * rel_fs_m
      ev_fs_xm <- rel_fs_x^2 * ev_fs_m + rel_fs_m^2 * ev_fs_x + ev_fs_m * ev_fs_x # formula
    })
  
  return(fs_dat)
}

# ========================================= Data Analysis ========================================= #

# analyze_2spa <- function (condition, dat, fixed_objects = NULL) {
#   
#   # Prepare fs_dat
#   fs_dat <- generate_fsdat(dat)
#   
#   # OpenMX model
#   tspa_umx <- umxLav2RAM(
#     "
#       fs_y ~ b1 * fs_x + b2 * fs_m + b3 * fs_xm
#       fs_y + fs_x + fs_m + fs_xm ~ 1
#       fs_x ~~ vx * fs_x + fs_m + fs_xm
#       fs_m ~~ vm * fs_m + fs_xm
#     ",
#     printTab = FALSE)
#     # Loading
#     matL <- mxMatrix(
#       type = "Diag", nrow = 4, ncol = 4,
#       free = FALSE,
#       labels = c("data.rel_fs_y", "data.rel_fs_x", "data.rel_fs_m", "data.ld_fs_xm"),
#       name = "L"
#     )
#     # Error
#     matE <- mxMatrix(
#       type = "Diag", nrow = 4, ncol = 4,
#       free = FALSE,
#       labels = c("data.ev_fs_y", "data.ev_fs_x", "data.ev_fs_m", "data.ev_fs_xm"),
#       name = "E"
#     )
#     
#   # Run the OpenMX model
#   tspa_mx <-
#       tspa_mx_model(
#         mxModel(tspa_umx,
#                 mxAlgebra(b1 * sqrt(vx), name = "stdx_b1"),
#                 mxAlgebra(b2 * sqrt(vm), name = "stdx_b2"),
#                 mxAlgebra(b3 * sqrt(vx) * sqrt(vm), name = "stdx_b3"),
#                 mxCI(c("b1", "b2", "b3", "stdx_b1", "stdx_b2", "stdx_b3"))),
#         data = fs_dat,
#         mat_ld = matL, 
#         mat_ev = matE)
#   tspa_mx_fit <- mxRun(tspa_mx, intervals = TRUE)
# 
#   # Extract Parameters
#   est_ust <- tspa_mx_fit$output$estimate[c("b1", "b2", "b3")]
#   se_ust <- tspa_mx_fit$output$standardErrors[c("b1", "b2", "b3"), ]
#   est_std <- c(tspa_mx_fit$output$algebras$m1.stdx_b1,
#                tspa_mx_fit$output$algebras$m1.stdx_b2,
#                tspa_mx_fit$output$algebras$m1.stdx_b3)
#   se_std <- c(mxSE(m1.stdx_b1, model = tspa_mx_fit),
#               mxSE(m1.stdx_b2, model = tspa_mx_fit),
#               mxSE(m1.stdx_b3, model = tspa_mx_fit))
#   
#   # Create the output vector
#   out <- c(est_ust, se_ust, est_std, se_std)  
#   names(out) <- c("x_est_ust", "m_est_ust", "xm_est_ust",
#                   "x_se_ust", "m_se_ust", "xm_se_ust", 
#                   "x_est_std", "m_est_std", "xm_est_std",
#                   "x_se_std", "m_se_std", "xm_se_std")
#   
#   # Return
#   return(out)
# }

analyze_upi <- function (condition, dat, fixed_objects = NULL) {
  
  fit_upi <- upi(model = fixed_objects$upi_syntax, 
                 data = dat, 
                 mode = "all",
                 ordered = c("x1", "x2", "x3", "m1", "m2", "m3",
                             "m4", "m5", "m6", "m7", "m8", "m9",
                             "m10", "m11", "m12")) 
  param_upi <- parameterEstimates(fit_upi)
  
  # Extract parameters
  est_ust <- param_upi %>% 
    filter(label %in% c("b1", "b2", "b3")) %>% 
    pull(est)
  se_ust <- param_upi %>% 
    filter(label %in% c("b1", "b2", "b3")) %>% 
    pull(se)
  est_std <- param_upi %>% 
    filter(label %in% c("stdx_b1", "stdx_b2", "stdx_b3")) %>% 
    pull(est)
  se_std <- param_upi %>% 
    filter(label %in% c("stdx_b1", "stdx_b2", "stdx_b3")) %>% 
    pull(se)
  
  # Create the output vector
  out <- c(est_ust, se_ust, est_std, se_std)  
  names(out) <- c("x_est_ust", "m_est_ust", "xm_est_ust",
                  "x_se_ust", "m_se_ust", "xm_se_ust", 
                  "x_est_std", "m_est_std", "xm_est_std",
                  "x_se_std", "m_se_std", "xm_se_std")
  
  # Return
  return(out)
}

analyze_2spamplus <- function (condition, dat, fixed_objects = NULL) {
  
  fs_dat <- generate_fsdat(dat)
  
  # Run Mplus model
  temp_session <- fixed_objects$temp_session
  dat_name <- file.path(temp_session,
                        sprintf("2spa_simdat_%s.dat", condition$REPLICATION))
  write.table(fs_dat,
              file = dat_name,
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  mplus_files <- sprintf(
    c("2spa_simdat_%s.dat", "2spa_%s.inp", "2spa_%s.out", "2spa_results_%s.dat"),
    condition$REPLICATION
  )
  writeLines(gsub("repid", replacement = condition$REPLICATION,
                  x = fixed_objects$tspa_syntax),
             con = file.path(temp_session, mplus_files[[2]]))
  MplusAutomation::runModels(
    temp_session,
    filefilter = mplus_files[[2]],
    Mplus_command = "/opt/mplus/8.7/mplus")
  
  # Extract results
  # Read lines
  res_2spa <- as.numeric(unlist(strsplit(trimws(readLines(file.path(temp_session, mplus_files[[4]]))), "\\s+")))
  # readLines(file.path(temp_session, mplus_files[[3]]))
  
  # Extract Parameters
  est_ust <- res_2spa[13:15]
  se_ust <- res_2spa[38:40]
  est_std <- res_2spa[23:25]
  se_std <- res_2spa[48:50]
  
  # Create the output vector
  out <- c(est_ust, se_ust, est_std, se_std)  
  names(out) <- c("x_est_ust", "m_est_ust", "xm_est_ust",
                  "x_se_ust", "m_se_ust", "xm_se_ust", 
                  "x_est_std", "m_est_std", "xm_est_std",
                  "x_se_std", "m_se_std", "xm_se_std")
  
  # Remove temp files
  file.remove(file.path(temp_session, mplus_files), recursive = TRUE)
  
  # Return
  return(out)
}

analyze_lms <- function (condition, dat, fixed_objects = NULL) {

  temp_session <- fixed_objects$temp_session
  dat_name <- file.path(temp_session,
                        sprintf("lms_simdat_%s.dat", condition$REPLICATION))
  write.table(dat,
              file = dat_name,
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  mplus_files <- sprintf(
    c("lms_simdat_%s.dat", "lms_%s.inp", "lms_%s.out", "lms_results_%s.dat"),
    condition$REPLICATION
  )
  writeLines(gsub("repid", replacement = condition$REPLICATION,
                  x = fixed_objects$lms_syntax),
             con = file.path(temp_session, mplus_files[[2]]))
  MplusAutomation::runModels(
    temp_session,
    filefilter = mplus_files[[2]],
    Mplus_command = "/opt/mplus/8.7/mplus")
  
  # Extract results
  # Read lines
  res_lms <- as.numeric(unlist(strsplit(trimws(readLines(file.path(temp_session, mplus_files[[4]]))), "\\s+")))
  # readLines(file.path(temp_session, mplus_files[[3]]))
  
  est_ust <- res_lms[25:27]
  se_ust <- res_lms[105:107]
  est_std <- res_lms[185:187]
  se_std <- res_lms[265:267]
  
  # Create the output vector
  out <- c(est_ust, se_ust, est_std, se_std)  
  names(out) <- c("x_est_ust", "m_est_ust", "xm_est_ust",
                  "x_se_ust", "m_se_ust", "xm_se_ust", 
                  "x_est_std", "m_est_std", "xm_est_std",
                  "x_se_std", "m_se_std", "xm_se_std")
  
  # Remove temp files
  file.remove(file.path(temp_session, mplus_files), recursive = TRUE)
  
  # Return
  return(out)
}

analyze_lmsfs <- function (condition, dat, fixed_objects = NULL) {
  
  fs_dat <- generate_fsdat(dat)
  
  temp_session <- fixed_objects$temp_session
  dat_name <- file.path(temp_session,
                        sprintf("lmsfs_simdat_%s.dat", condition$REPLICATION))
  write.table(fs_dat,
              file = dat_name,
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  mplus_files <- sprintf(
    c("lmsfs_simdat_%s.dat", "lmsfs_%s.inp", "lmsfs_%s.out", "lmsfs_results_%s.dat"),
    condition$REPLICATION
  )
  writeLines(gsub("repid", replacement = condition$REPLICATION,
                  x = fixed_objects$lmsfs_syntax),
             con = file.path(temp_session, mplus_files[[2]]))
  MplusAutomation::runModels(
    temp_session,
    filefilter = mplus_files[[2]],
    Mplus_command = "/opt/mplus/8.7/mplus")
  
  # Extract results
  # Read lines
  res_lmsfs <- as.numeric(unlist(strsplit(trimws(readLines(file.path(temp_session, mplus_files[[4]]))), "\\s+")))
  # readLines(file.path(temp_session, mplus_files[[3]]))
  
  est_ust <- res_lmsfs[10:12]
  se_ust <- res_lmsfs[29:31]
  est_std <- res_lmsfs[17:19]
  se_std <- res_lmsfs[36:38]
  
  # Create the output vector
  out <- c(est_ust, se_ust, est_std, se_std)  
  names(out) <- c("x_est_ust", "m_est_ust", "xm_est_ust",
                  "x_se_ust", "m_se_ust", "xm_se_ust", 
                  "x_est_std", "m_est_std", "xm_est_std",
                  "x_se_std", "m_se_std", "xm_se_std")
  
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
ci_stats <- function(est, se, par, stats_type) {
  
  # Calculate the confidence intervals (usd)
  ll <- est - qnorm(.975) * se
  ul <- est + qnorm(.975) * se
  ci <- vector("list", length = ncol(est))
  names(ci) <- colnames(est)
  
  # Construct confidence intervals for each method
  for (i in seq_len(ncol(est))) {
    ci[[i]] <- cbind(ll[,i], ul[,i])
  }

  # Determine which statistic to calculate
  if (stats_type == "Coverage") {
    return(sapply(ci, function(ci) mean(ci[,1] <= par & ci[,2] >= par)))
  } else if (stats_type == "TypeI") {
    return(sapply(ci, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Power") {
    return(sapply(ci, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
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
  est_std <- results[, grep("xm_est_std", colnames(results)), drop = FALSE]
  est_ust <- results[, grep("xm_est_ust", colnames(results)), drop = FALSE]
  se_std <- results[, grep("xm_se_std", colnames(results)), drop = FALSE]
  se_ust <- results[, grep("xm_se_ust", colnames(results)), drop = FALSE]

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
    coverage = ci_stats(est_std, se_std, pop_par,
                        "Coverage"),
    type1_std = ci_stats(est_std, se_std, pop_par,
                         "TypeI"),
    power_std = ci_stats(est_std, se_std, pop_par,
                          "Power"),
    rmse = RMSE(na.omit(est_std),
                parameter = pop_par)
    # warning_total = warning_sum(warnings)
  )
}

# ========================================= Run Experiment ========================================= #
# Trial
# trial <- runSimulation(design = DESIGNFACTOR[1,],
#               replications = 3,
#               generate = generate_dat,
#               analyse = list(upi = analyze_upi,
#                              tspa = analyze_2spamplus,
#                              lms = analyze_lms,
#                              lmsfs = analyze_lmsfs),
#               summarise = evaluate_res,
#               fixed_objects = FIXED,
#               seed = rep(61543, nrow(DESIGNFACTOR[1,])),
#               control = list(include_reps = TRUE),
#               packages = "lavaan",
#               parallel = TRUE,
#               ncores = 30)

runSimulation(design = DESIGNFACTOR,
              replications = 2000,
              generate = generate_dat,
              analyse = list(upi = analyze_upi,
                             tspa = analyze_2spamplus,
                             lms = analyze_lms,
                             lmsfs = analyze_lmsfs),
              summarise = evaluate_res,
              fixed_objects = FIXED,
              seed = rep(61543, nrow(DESIGNFACTOR)),
              control = list(include_reps = TRUE),
              packages = "lavaan", 
              filename = "categorical_03252025",
              parallel = TRUE,
              ncores = 30,
              save = TRUE,
              save_results = TRUE)
