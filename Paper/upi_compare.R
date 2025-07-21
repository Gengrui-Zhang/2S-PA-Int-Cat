library(modsem)

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
    beta1 := b1*sqrt(v1)
    beta2 := b2*sqrt(v2)
    beta3 := b3*sqrt(v1)*sqrt(v2)
"

condition <- DESIGNFACTOR[6,]
fixed_objects <- FIXED
dat <- generate_dat(condition, fixed_objects)

fit_upi <- upi(model = "
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
    beta1 := b1*sqrt(v1)
    beta2 := b2*sqrt(v2)
    beta3 := b3*sqrt(v1)*sqrt(v2)
", 
               data = dat, 
               mode = "all",
               ordered = c("x1", "x2", "x3", "m1", "m2", "m3",
                           "m4", "m5", "m6", "m7", "m8", "m9",
                           "m10", "m11", "m12")) 
summary(fit_upi)

fit_modsem <- modsem(
"# Measurement Model
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
  beta1 := b1*sqrt(v1)
  beta2 := b2*sqrt(v2)
  beta3 := b3*sqrt(v1)*sqrt(v2)
  ",
data = dat, 
method = "dblcent") 
summary(fit_modsem)

# ------------------------------------------------------------------------------------

condition <- DESIGNFACTOR[12,]
fixed_objects <- FIXED_PARAMETER
dat <- generate_dat(condition, fixed_objects)

fit_upi <- upi(model = fixed_objects$model, 
               data = dat, 
               mode = "all") 
param_upi <- parameterEstimates(fit_upi)
summary(fit_upi)

fit_modsem <- modsem(fixed_objects$model, data = dat, method = "dblcent") 
summary(fit_modsem)

# PIs
intNames_x <- c("x1", "x2", "x3")
intNames_m <- c("m1", "m2", "m3")
mod_df <- indProd(dat %>% select(-Y), 
                  var1 = intNames_x,
                  var2 = intNames_m,
                  match = FALSE, 
                  meanC = T, 
                  residualC = F, 
                  doubleMC = T,
                  namesProd = as.vector(outer(intNames_x, intNames_m, paste0)))
upi_df <- cbind(mod_df, dat %>% select(Y))
model <- "# Measurement Model
X =~ x1 + x2 + x3
M =~ m1 + m2 + m3
Int =~ x1m1 + x2m1 + x3m1 + x1m2 + x2m2 + x3m2 + x1m3 + x2m3 + x3m3
# Structural Model
Y ~ b1*X + b2*M + b3*Int
# Covariance
x1m1 ~~ th1*x1m2 + th1*x1m3 
x1m2 ~~ th1*x1m3 
x2m1 ~~ th2*x2m2 + th2*x2m3 
x2m2 ~~ th2*x2m3 
x3m1 ~~ th3*x3m2 + th3*x3m3 
x3m2 ~~ th3*x3m3

x1m1 ~~ th4*x2m1 + th4*x3m1 
x2m1 ~~ th4*x3m1 
x1m2 ~~ th5*x2m2 + th5*x3m2 
x2m2 ~~ th5*x3m2 
x1m3 ~~ th6*x2m3 + th6*x3m3 
x2m3 ~~ th6*x3m3
# Define Standardized Coefficients
X ~~ v1*X
M ~~ v2*M
beta1 := b1*sqrt(v1)
beta2 := b2*sqrt(v2)
beta3 := b3*sqrt(v1)*sqrt(v2)"
fit_upi_raw <- sem(model, upi_df)
parameterEstimates(fit_upi_raw)
summary(fit_upi_raw)

model_modsem <- "# Measurement Model
            X =~ x1 + x2 + x3
            M =~ m1 + m2 + m3
            Int =~ x1m1 + x2m1 + x3m1 + x1m2 + x2m2 + x3m2 + x1m3 + x2m3 + x3m3
          # Structural Model
            Y ~ b1*X + b2*M + b3*Int
          # Covariance
            x1m1 ~~ 0*x2m2 + 0*x2m3 + 0*x3m2 + 0*x3m3
            x1m2 ~~ 0*x2m3 + 0*x3m3
            x2m1 ~~ 0*x1m2 + 0*x1m3 + 0*x3m2 + 0*x3m3 
            x2m2 ~~ 0*x1m3 + 0*x3m3
            x3m1 ~~ 0*x1m2 + 0*x1m3 + 0*x2m2 + 0*x2m3
            x3m2 ~~ 0*x1m3 + 0*x2m3
            
            x1m1 ~~ x1m2 + x1m3 + x2m1 + x3m1
            x1m2 ~~ x1m3 + x2m2 + x3m2
            x1m3 ~~ x2m3 + x3m3
            x2m1 ~~ x2m2 + x2m3 + x3m1
            x2m2 ~~ x2m3 + x3m2
            x2m3 ~~ x3m3
            x3m1 ~~ x3m2 + x3m3
            x3m2 ~~ x3m3

          # Define Standardized Coefficients
            X ~~ v1*X
            M ~~ v2*M
            beta1 := b1*sqrt(v1)
            beta2 := b2*sqrt(v2)
            beta3 := b3*sqrt(v1)*sqrt(v2)"
fit_upi_modsem <- sem(model_modsem, upi_df)
parameterEstimates(fit_upi_modsem)
summary(fit_upi_modsem)

# ------------------------------------------------------------

mod_upi <- "
          X =~ 1*x1 + x2 + x3
          M =~ 1*m1 + m2 + m3
          Int =~ x1m1 + x1m2 + x1m3 + x1m4 + x1m5 + x1m6 + 
                 x1m7 + x1m8 + x1m9 + x1m10 + x1m11 + x1m12 + 
                 x2m1 + x2m2 + x2m3 + x2m4 + x2m5 + x2m6 + 
                 x2m7 + x2m8 + x2m9 + x2m10 + x2m11 + x2m12 + 
                 x3m1 + x3m2 + x3m3 + x3m4 + x3m5 + x3m6 + 
                 x3m7 + x3m8 + x3m9 + x3m10 + x3m11 + x3m12
          Y =~ 1*y1 + y2 + y3

          Y ~ b1*X + b2*M + b3*Int
              
          # Variance
          X ~~ v1*X
          M ~~ v2*M
          Int ~~ v3*Int
          
          # Covariance
          X ~~ cov_12*M
          X ~~ cov_13*Int
          M ~~ cov_23*Int
          
          # Disturbance
    			Y ~~ dist_y*Y
          var_y :=  dist_y + (b1^2 * v1 + b2^2 * v2 + b3^2 * v3 + 
                              2 * b1 * b2 * cov_12 + 
                              2 * b1 * b3 * cov_13 + 
                              2 * b2 * b3 * cov_23)
          
          beta1 := b1*sqrt(v1)/sqrt(var_y)
          beta2 := b2*sqrt(v2)/sqrt(var_y)
          beta3 := b3*sqrt(v1)*sqrt(v2)/sqrt(var_y)
          "
fit_upi <- sem(mod_upi, 
               pi_dat, 
               estimator = "MLR")

summary(fit_upi)
