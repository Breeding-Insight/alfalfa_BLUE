library(tidyverse)
library(asremlPlus)
library(asreml)

# Read in data
df <- read.csv("data/raw/input.csv") %>%
  select(Entry, Plot, Col, Row, Rep, spread_cm) %>%
  mutate(Rep = if_else(Rep == "na", "5", Rep)) %>%
  mutate(spread_cm = na_if(spread_cm, ".")) %>%
  mutate(spread_cm = as.numeric(spread_cm)) %>%
  mutate(across(c(Entry, Plot, Col, Row, Rep), as.factor)) %>%
    arrange(Col, Row)

####################################
#          Create BLUE             #
####################################

# Examine different model structures to find the best model for BLUE
# Lets test run a blue0 for var comp, this is a simple model with no spatial structure
blue0 <- asreml(fixed = spread_cm ~ Entry , # entry is fixed = blue
                 random = ~ Rep + Row + Col,
                 residual = ~ idv(units),  # Basic nugget structure for residual
                 data = df,
                 na.action = na.method(y = "omit", x = "omit"),
                 workspace = "8gb")

var0 <- summary(blue0)$varcomp # Lets analyze variance components to have an idea how to move forward:
# All have low z ratios, this is a small experiment! Since z raios are lower than 2 we should probably avoid to many random effects
# Rep has 5 cats (n=5, the minimum number for me to consider it as random)
# Rep also has a low z ratio (<2), but we need to keep it in the model to retain replication for our predictions, lets make it fixed
# Row Col have small z ratios, suggesting they may not be needed
# However - Units! is high, this tells me we have unmodeled spatial variance
# As breeders we know even tho they are not siginificant we should keep them in the model, high units! supports this
# Rather than retain as reandom effect, lets model our residual structure around row / col to account for this

# Stock model - simple BLUE, notice genotype (Entry) is fixed, for demo purposes lets start here then implment our findings above
blue1 <- asreml(fixed = spread_cm ~ Entry + Rep,
                residual = ~ idv(units),  # Basic nugget structure for residual
                data = df,
                na.action = na.method(y = "omit", x = "omit"),
                workspace = "8gb")

test1 <- infoCriteria(blue1, IClikelihood = "full")  # Check for convergence, test AIC (best predictor) and BIC (truest model, often simplest)

# First order autoregressive structure for row col - this captures all the variation you were wanting with a BLUP (and then some, its ***fancy***)
blue2 <- asreml(fixed = spread_cm ~ Entry + Rep,
                residual = ~ ar1v(Col):ar1(Row), # Test 2D AR structure, this may be overkill
                data = df,
                na.action = na.method(y = "include", x = "include") # include all data to maintain data structure to match data structure, you can use these in the gwas! Similiar to marker imputation
)

var2 <- summary(blue2)$varcomp # Notice!! our row col cors are low, the AR structure is overkill and introduced unneeded complexity, but it was a good test

test2 <- infoCriteria(blue2, IClikelihood = "full")

# independent residual structure for row col - row and col effects are discrete. We are "dumbing down" our AR structure
blue3 <- asreml(fixed = spread_cm ~ Entry + Rep,
                residual = ~ idv(Col):id(Row), # notice, we have more col heterogeniety than row, so lets use idv (heterogeneius var) for col and id (homogeneous var) for row
                data = df,
                na.action = na.method(y = "include", x = "include") # include all data to maintain data structure to match data structure, you can use these in the gwas! Similiar to marker imputation
)

var3 <- summary(blue3)$varcomp # Bingo, BIC loves it.
test3 <- infoCriteria(blue3, IClikelihood = "full")

# Lets test some elevated models for fun - nugget effect in random
blue4 <- asreml(fixed = spread_cm ~ Entry + Rep,
                random = ~ units,
                residual = ~ idv(Col):id(Row), 
                data = df,
                na.action = na.method(y = "include", x = "include"),
                ai.sing = TRUE, # We have a singularity, units and residual strucuture are not independent, lets forget this one
                workspace = "8gb")

test3 <- infoCriteria(blue4, IClikelihood = "full") # Check convergence and residual variance

# Elevated autoregressive structure with added random effect for Col and Row, residual is parsed accordingly with our AR structure, this is the most complex model for test purposes
blue5 <- asreml(fixed = spread_cm ~ Entry + Rep,
                random = ~ Col + Row + idv(units),
                residual = ~ ar1v(Col):ar1(Row),
                data = df,
                na.action = na.method(y = "include", x = "include"),
                workspace = "8gb")
 
var5 <- summary(blue5)$varcomp # Low z ratios in these random effects, probably best to not include them and avoid redundancy with spatial resiudals
test5 <- infoCriteria(blue5, IClikelihood = "full") # Check convergence and residual variance

# Evaluate all models
test <- rbind(test1, test2, test3, test4, test5) # Options 2 (lowest BIC) and 5 (Lowest AIC) are best

# Create BLUE with best model
# Lets use blue3 w/lowest BIC, idepentent residuals, no random effect redundancy, simple and accurate, best for this design and objective
# If this were a breeding situation, I would pursue AIC for accuracy, but since this is GWAS phase one, lets put a higher precident on BIC
# Our analysis of variance componenets also supports this decision
# Below is how we predict models using our "winner" model
blue_pred <- as_tibble(predict.asreml(blue3, classify = "Entry:Rep")$pvals) # classify entry within rep to retain replication as Meng wanted

write.csv(blue_pred, "/Users/aja294/Documents/alfalfa/data/processed/blues_23.csv", row.names = FALSE) # Save predictions for later use

####################################
#     REPEAT FOR WEIGHTED BLUE     #
####################################

# Caculate reliability or prediction to carry over as var-covar mat into next step - Best stastical practice
df_weight <- df %>%
  left_join(blue_pred, by = c("Entry", "Rep"), multiple = "all") %>%
  rename(spread = spread_cm) %>%
  mutate(Site = as.factor("year23")) %>%
  mutate(wt = 1 / std.error^2) # inverse of SE = reliability per ob to carry into next model

# Perform balanced BLUE using BLUE option 2, examine family and weights arguements for new var/covar matrices
blue3 <- asreml(fixed = spread ~ Entry + Rep, # Retain rep, we still need to account for systematic block effects if we are to average over it
                weights = wt, # We carry over var/cov matrix from previous model, this is the reliability of each observation
                family = asr_gaussian(dispersion = 1), # Gaussian family parameter is required for weights, this allows weight to soley determine observation reliability
                residual = ~ units, # Lets not double count our row col residuals, we already accounted for this in the previous model, instead lets use a basic homogeneous var with unit nugget
                data = df_weight,
                workspace = "8gb",
                na.action = na.method(y = "omit", x = "omit"))

blue_balanced <-  as_tibble(predict.asreml(blue2, classify = "Entry:Rep")$pvals) # classify entry within rep to retain replication as Meng wanted

write.csv(blue_balanced, "/Users/aja294/Documents/alfalfa/data/processed/balanced_blue_23.csv", row.names = FALSE) # Save predictions for later use

####################################
#         REPEAT FOR 5-13          #
####################################

# Repeat for spread ratings on 5-13, which seems to be the most reliable of the two ratings taken that year
# Less missing data in 5-13, fewer "0" values, and more consistent spread ratings
df2 <- read.csv("data/raw/input.csv") %>%
  select(Entry, Plot, Col, Row, Rep, spread_5.13) %>%
  mutate(spread_5.13 = na_if(spread_5.13, "."), Rep = if_else(Rep == "na", "5", Rep)) %>%
  mutate(spread_5.13 = as.numeric(spread_5.13)) %>%
  mutate(across(c(Entry, Plot, Col, Row, Rep), as.factor)) %>%
  arrange(Col, Row)

blue2 <- asreml(fixed = spread_5.13 ~ Entry + Rep,
                residual = ~ idv(Col):id(Row),
                data = df2,
                na.action = na.method(y = "include", x = "include") # include all data to maintain data structure to match data structure, you can use these in the gwas! Similiar to marker imputation
)

test2 <- infoCriteria(blue2, IClikelihood = "full")

blue_5.13 <- as_tibble(predict.asreml(blue2, classify = "Entry:Rep")$pvals) # classify entry within rep to retain replication as Meng wanted

write.csv(blue_5.13, "/Users/aja294/Documents/alfalfa/data/processed/blues_22.csv", row.names = FALSE) # Save predictions for later use

df2_weight <- df2 %>%
  left_join(blue_5.13, by = c("Entry", "Rep"), multiple = "all") %>%
  mutate(wt = 1 / std.error^2) # inverse of SE = reliability per ob to carry into next model

blue2 <- asreml(fixed = spread_5.13 ~ Entry + Rep,
                weights = wt, # We carry over var/cov matrix from previous model, this is the reliability of each observation
                family = asr_gaussian(dispersion = 1), # Gaussian family parameter is required for weights
                residual = ~ units,
                data = df2_weight,
                workspace = "8gb",
                na.action = na.method(y = "include", x = "include"))

blue_balanced <-  as_tibble(predict.asreml(blue2, classify = "Entry:Rep")$pvals) # classify entry within rep to retain replication as Meng wanted

write.csv(blue_balanced, "/Users/aja294/Documents/alfalfa/data/processed/balanced_blue_22.csv", row.names = FALSE) # Save predictions for later use


####################################
#       Create Multi-Env          #
####################################

# THIS IS BEST STATISTICAL PRACTICE

df3_multi <- df2_weight %>%
  rename(spread = spread_5.13) %>%
  mutate(spread = spread * 10) %>%
  mutate(Site = as.factor("year22")) %>%
  bind_rows(df_weight) %>%
  select(Entry, Site, Plot, Col, Row, Rep, spread, wt)

# We already have weights from our last models, so start at phase 2 (balancing) here

blue_multi <- asreml(fixed = spread ~ Entry + Site + Entry:Site, # Interaction included but genotype remains fixed
                     random = ~ at(Site):Rep, # Nested variance structure for rep within loc, Rep must be random now
                     weights = wt, # We carry over var/cov matrix from previous model, this is the reliability of each observation
                     family = asr_gaussian(dispersion = 1), # Gaussian family parameter is required for weights
                     residual = ~ units,
                     data = df3_multi,
                     workspace = "8gb",
                     na.action = na.method(y = "omit", x = "omit"))

pred_blue_multi <- as_tibble(predict.asreml(blue_multi, classify = "Entry:Rep:Site")$pvals) # I reccomend using these BLUES for GWAS, they are the most accurate

var_multi <- summary(blue_final)$varcomp
test_multi <- infoCriteria(blue_final, IClikelihood = "full")
# AIC and BIC will be higher since this is ***much more complex***,
# and should not be compared with relativity to the previous versions,
# but we know from initial studies and sound statistical theory
# that this structure is the best approach.
write.csv(pred_blue_multi, "/Users/aja294/Documents/alfalfa/data/processed/multi_env_balanced_BLUE.csv", row.names = FALSE)
