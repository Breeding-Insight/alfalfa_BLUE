##setwd("~/BIG/alfalfa")

library(tidyverse)
library(asremlPlus)
library(asreml)
library(tidyr)

# Read in data
df <- read.csv("input_2.csv")
head(df)
# Check structure of the data to ensure it was read in correctly
df1 <- df %>%
  select(Entry, Plot, Row, Col, rep, geno_id, spread_22, spread_23) %>%
  mutate(
    spread_22 = as.numeric(na_if(spread_22, ".")),
    spread_23 = as.numeric(na_if(spread_23, ".")),
    across(c(Entry, Plot, Col, Row, ,geno_id, rep), as.factor)
  ) %>%
  arrange(Col, Row) %>%
  pivot_longer(
    cols = starts_with("spread_"),
    names_to = "year",
    names_prefix = "spread_",
    values_to = "spread"
  )

df1$year<-as.factor(df1$year)
# Save to CSV
write.csv(df1, "input_3.csv", row.names = FALSE)
head(df1)

####################################
#          Single EXP analysis   #
####################################
df2 <- df1 %>%
  filter(year == "23") # subset data by "year" to get 2023 data.

#Simple model with spatial correction with Row and Col effects
blue0 <- asreml(fixed = spread ~ geno_id , # entry is fixed = blue
                 random = ~ Row + Col, # Use a simple Row and Col random effects to capture spatial variation
                 residual = ~ idv(units),  # Basic nugget structure for residual without R- structure
                 data = df2,
                 na.action = na.method(y = "include", x = "include"), # include all data to maintain structure for various level of analyses  
                 workspace = "8gb")

blue0<-update(blue0)
var0 <- summary(blue0)$varcomp # Lets analyze variance components to have an idea how to move forward:
var0

# Second model with auto-regressive correlation structure for spatial correction
blue1 <- asreml(fixed = spread ~ geno_id , 
                random = ~ ar1(Row):ar1v(Col), #auto-regressive order-1 structure for row and col, with separate col var 
                residual = ~ idv(units),  
                data = df2,
                na.action = na.method(y = "include", x = "include"),
                workspace = "8gb")


var1<-summary(blue1)$varcomp # Check variance components for the updated model
var1

blue_pred_2023 <- as_tibble(predict.asreml(blue1, classify = "geno_id")$pvals) %>%
  select(geno_id, predicted.value,status, std.error) %>%
  mutate(wt = 1 / (std.error^2),year="2023")

#blue_pred_2023 <- as_tibble(predict.asreml(blue1, classify = "geno_id")$pvals)
write.csv(blue_pred_2023, "blue_pred_2023.csv", row.names = FALSE)

#####Evaluate models with AIC, BIC and Loglikelihood values and tests
summary(blue0)$aic
summary(blue1)$aic 
summary(blue0)$bic
summary(blue1)$bic 
logLik_model1<-summary(blue0)$loglik
logLik_model2<-summary(blue1)$loglik
lrt_1_vs_2 <- 2 * (logLik_model2 - logLik_model1)
p_value_1_vs_2 <- pchisq(lrt_1_vs_2, df = 1, lower.tail = FALSE)
cat("\nLikelihood Ratio Test Results:\n")
cat("Model 1 vs Model 2: Chi-square =", lrt_1_vs_2, " P-value =", p_value_1_vs_2, "\n")

########For both 2022 and 2023 spread data, model "blue1" seems to be a better fit
########Get BLUEs for 22 and 23##################################################
df3 <- df1 %>%
  filter(year == "22")

blue1_1 <- asreml(fixed = spread ~ geno_id , 
                random = ~ ar1(Row):ar1v(Col), #auto-regressive order-1 structure for row and col, with separate col var 
                residual = ~ idv(units),  
                data = df3,
                na.action = na.method(y = "include", x = "include"),
                workspace = "8gb")

blue1_1<-update(blue1_1)
var1_1<-summary(blue1_1)$varcomp # Check variance components for the updated model
var1_1

blue_pred_2022 <- as_tibble(predict.asreml(blue1_1, classify = "geno_id")$pvals) %>%
  select(geno_id, predicted.value,status, std.error) %>%
  mutate(wt = 1 / (std.error^2), year="2022")

#blue_pred_2022 <- as_tibble(predict.asreml(blue1_1, classify = "geno_id")$pvals)
write.csv(blue_pred_2022, "blue_pred_2022.csv", row.names = FALSE) 

##########Combined analysis data for 2022 and 2023##############################
##########Single stage analysis- by combining year- very simple and resolvable##

head(df1)
blue_combined <- asreml(fixed = spread ~ geno_id + year, # genotype as fixed corrected with year effect, which is also fixed
                  random = ~ at(year):ar1(Row):ar1v(Col), #auto-regressive order-1 structure for row and col, with separate col variance 
                  residual = ~ idv(units),  
                  data = df1,
                  na.action = na.method(y = "include", x = "include"),
                  workspace = "8gb")

blue_combined<-update(blue_combined)
summary(blue_combined)$varcomp

blue_combined_pred <- as_tibble(predict.asreml(blue_combined, classify = "geno_id")$pvals)
write.csv(blue_combined_pred, "blue_combined_pred.csv", row.names = FALSE)


##################################### Combined analysis with weights ######
############# Carry the weight from first stage analysis in terms of 1/(standard error^2)#### 
############# use BLUE from stage one as input data with #######################
##############units variance constraint to 1####################################

weighted_blue_dataset <- bind_rows(blue_pred_2022, blue_pred_2023) %>%
  relocate(year, .after = geno_id)

write.csv(weighted_blue_dataset, "weighted_blue_input_dataset.csv", row.names = FALSE)

str(weighted_blue_dataset)
head(weighted_blue_dataset)
weighted_blue_dataset$year <- as.factor(weighted_blue_dataset$year)

###################Fit the second stage model


weighted_blue_model <- asreml(fixed = predicted.value ~ geno_id + year, # year effect corrected
                     weights = wt, #weight for each observation carried forward from first stage analysis
                     family = asr_gaussian(dispersion = 1), # Gaussian family parameter is required for weights
                     residual = ~ units, #Constrained to 1
                     data = weighted_blue_dataset,
                     workspace = "8gb",
                     na.action = na.method(y = "omit", x = "omit"))

summary(weighted_blue_model)$varcomp

pred_weighted_blue <- as_tibble(predict.asreml(weighted_blue_model, classify = "geno_id")$pvals) # 
write.csv(pred_weighted_blue, "pred_weighted_blue_by_geontype.csv", row.names = FALSE)

pred_weighted_blue2 <- as_tibble(predict.asreml(weighted_blue_model, classify = "geno_id:year")$pvals)
write.csv(pred_weighted_blue2, "pred_weighted_blue_by_geontype_year.csv", row.names = FALSE)



