#_______________Load libraries ___________________

library(asreml)
library(tidyverse)

#____________Upload and tidy data ________________


df <- read.csv("data/raw/input_meng.csv")
head(df)

# Check structure of the data to ensure it was read in correctly
df1 <- df %>%
  # Select only needed columns
  select(Entry, Plot, Row, Col, Rep, geno_id_Meng, spread_22, spread_23) %>%
  # Convert spread columns to numeric and handle missing values in one step
  mutate(
    across(starts_with("spread_"), ~as.numeric(na_if(., "."))),
    # Convert multiple columns to factors at once
    across(c(Entry, Plot, Col, Row, geno_id_Meng, Rep), as.factor)
  ) %>%
  # Pivot before arranging for better performance with long format
  pivot_longer(
    cols = starts_with("spread_"),
    names_to = "year",
    names_prefix = "spread_",
    values_to = "spread"
  ) %>%
  # Convert year to factor and arrange data
  mutate(year = as.factor(year)) %>%
  arrange(Col, Row, year)

df1_22 <- df1 %>% filter(year == "22") # subset data by "year" to get 2022 data.
df1_23 <- df1 %>% filter(year == "23") # subset data by "year" to get 2023 data.

####################################
# Single Stage "Combined" Analysis #
####################################
#
# NOTE: Spatial correction performed in
# G structure due to experiement inbalances
# between years.
#
# NOTE: Model Log-likelihood not converged,
# algorithmic stability is not reached,
# this is due to dataset limitations not model shortcomings.

blue_1phase <- asreml(fixed = spread ~ geno_id_Meng + year,
                      random = ~ at(year):ar1(Row):ar1v(Col),
                      residual = ~ idv(units),
                      data = df1,
                      na.action = na.method(y = "omit", x = "omit"),
                      workspace = "8gb")

blue_1phase_pred <- as_tibble(predict.asreml(blue_1phase, classify = "geno_id_Meng")$pvals) %>%
  mutate(
    predicted.value = case_when(
      predicted.value < 0 ~ 0,
      predicted.value > 10 ~ 10,
      TRUE ~ predicted.value
    )
  )

######################
# Two-phase analysis #
######################

# NOTE: spatial correction performed in R structure (residual)

#_______Single blue for '22___________
#
# Note: This model reaches convergence,
# and algorithmic stability is achieved.

blue22 <- asreml(fixed = spread ~ geno_id_Meng,
                 residual = ~ idv(Col):id(Row),
                 data = df1_22,
                 na.action = na.method(y = "include", x = "include"))

blue22_pred <- as_tibble(predict.asreml(blue22, classify = "geno_id_Meng")$pvals) %>%
  mutate(year = as.factor(22))

#_______Single blue for '23__________
#
# Note: This model reaches convergence,
# and algorithmic stability is achieved.

blue23 <- asreml(fixed = spread ~ geno_id_Meng,
                 residual = ~ idv(Col):id(Row),
                 data = df1_23,
                 na.action = na.method(y = "include", x = "include"))

blue23_pred <- as_tibble(predict.asreml(blue23, classify = "geno_id_Meng")$pvals) %>%
  mutate(year = as.factor(23))

#_________Combine BLUEs for '22 and '23___________
#
# Extract reliability (precision) and
# carry over var/covar matrix as weights
# from first stage analysis to second stage analysis.
#
# NOTE: The weights are calculated as the inverse of the variance

df_weight <- blue23_pred %>%
  bind_rows(blue22_pred) %>%
  mutate(wt = 1 / std.error^2)

#________Weighted fixed effect model_____________
#
# NOTE: 2phase achieves perfect convergence
# and succesfully achieves algorithmic stability
# (6 iterations, No warnings).

blue_2phase <- asreml(fixed = predicted.value ~ geno_id_Meng + year,
                      weights = wt,
                      family = asr_gaussian(dispersion = 1),
                      residual = ~ idv(units),
                      data = df_weight,
                      workspace = "8gb",
                      na.action = na.method(y = "omit", x = "omit"))

# Print results
blue_2phase_pred <- as_tibble(predict.asreml(blue_2phase, classify = "geno_id_Meng")$pvals) %>%
  mutate(
    predicted.value = case_when(
      predicted.value < 0 ~ 0,
      predicted.value > 10 ~ 10,
      TRUE ~ predicted.value
    )
  )

blue_2phase_pred_year <- as_tibble(predict.asreml(blue_2phase, classify = "geno_id_Meng:year")$pvals) %>%
  mutate(
    predicted.value = case_when(
      predicted.value < 0 ~ 0,
      predicted.value > 10 ~ 10,
      TRUE ~ predicted.value
    )
  )

####################
# Model comparison #
####################

# Information criteria for 1 phase "combined" model
blue_1phase_var <- summary(blue_1phase)$varcomp
blue_1phase_bic <- summary(blue_1phase)$bic
blue_1phase_aic <- summary(blue_1phase)$aic
blue_1phase_loglik <- summary(blue_1phase)$loglik

# Information criteria for 2 phase model
blue_2phase_var <- summary(blue_2phase)$varcomp
blue_2phase_bic <- summary(blue_2phase)$bic
blue_2phase_aic <- summary(blue_2phase)$aic
blue_2phase_loglik <- summary(blue_2phase)$loglik

# Create table to compare bic, aic, loglik, and varcomp for both models
model_comparison <- tibble(
  Model = c("Two Phase", "One Phase"),
  BIC = c(blue_2phase_bic, blue_1phase_bic),
  AIC = c(blue_2phase_aic, blue_1phase_aic),
  LogLik = c(blue_2phase_loglik, blue_1phase_loglik),
  Parameters = c(length(blue_2phase$vparameters), length(blue_1phase$vparameters))
)

# _______________Print comparison table ______________________
#
# NOTE: The model with the lowest AIC/BIC is preferred - HOWEVER:
# LRT cannot be used to compare models with different nesting
# structures, e.g. (nested vs unested / 2 phase vs 1 phase).
# This phenomenon is largely due to the fact there is a
# significant difference in the degrees of freedom between the
# two models (df = 6).
#
# For comparison, predictive ability is preferred, but would be 
# very overkill for this situation.

print(model_comparison)

####################################################################
#                                                                  #
#          STASTICIAL RECCOMENDATION: 2 phase, weighted model      #
#  Reason: algorithmic stability achieved and convergence reached. #
#                                                                  #
####################################################################

# _______________Save results____________________
write.csv(blue_2phase_pred, "data/processed/blue_2phase_results.csv", row.names = FALSE)
write.csv(blue_2phase_pred_year, "data/processed/blue_2phase_by_year_results.csv", row.names = FALSE)