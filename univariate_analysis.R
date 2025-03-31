# About ========================================================================

# Title: Linear modelling for Maria Island Data (KSM308)
# Date: April 2025
# Author: Freddie Heather

# Imports data
# Visualises data
# Fits generalised linear models
# Generates figures

# Important: 
# 1. Make sure you are working inside an R project
# 2. Make sure the datafiles "multivariate_analysis.R" and "predictor_table.csv" are within your project working directory

# Packages =====================================================================

# install.packages("tidyverse") # run this line of code if you haven't installed the tidyverse on your computer
library(tidyverse)
# install.packages("performance", dependencies = TRUE) # for model diagnostics
library(performance)
# install.packages("MASS", dependencies = TRUE) # for the glm.nb function


# Import =======================================================================

# Prior to importing - have a look at your data in excel.

resp_vars <- read_csv("univariate_response_vars.csv") 
pred_vars <- read_csv("predictor_table.csv")

# Data checking ================================================================

resp_vars # has NA instead of zeros
pred_vars # has some blank columns

# Check to make sure the columns are of the right class (dbl or chr)

# Data wrangling ===============================================================

# mutate_at: replace SOME cols with NA values with zero
resp_vars_clean <- 
  resp_vars %>% 
  mutate_at(.vars = vars(`Acanthaluteres vittiger`:`Pictilabrus laticlavius`), 
            .funs = replace_na, 
            replace = 0) 

# mutate_all: replace ALL the NA values with zero
resp_vars_clean <- 
  resp_vars %>% 
  mutate_all(.funs = replace_na, replace = 0) 

pred_vars_clean <- 
  pred_vars %>% 
  discard(~all(is.na(.))) %>% # gets rid of those blank cols
  rename(Sample = sample) # to match the col-name in the other table

# joining the two datasets by the ID variable common to both ('Sample')
# this step is so that we have the data in the same dataframe so that we can 
# visualise the data and also perform stats on them
clean_data <- 
  resp_vars_clean %>% 
  left_join(pred_vars_clean, by = "Sample")

# Data checking ===============================================================

glimpse(clean_data)

# to convert a column from one data type (character to another):
# clean_data <- clean_data %>% mutate(x = as.numeric(x)) # to dbl
# clean_data <- clean_data %>% mutate(x = as.character(x)) # to chr

# EXAMPLE 1 (Numerical vs numerical) ===========================================

# STEP 1: visualise data -------------------------------------------------------
clean_data %>% 
  ggplot() + 
  aes(
    x = Depth, 
    y = Richness
  ) +
  geom_point()  + 
  stat_smooth(se = TRUE) + # smoothed trend line? standard error?
  labs(
    x = "Depth (m)",
    y = "Richness (# species)"
  ) +
  theme_bw() # themes from https://ggplot2.tidyverse.org/reference/ggtheme.html

# STEP 2: interpret ------------------------------------------------------------
# maybe a slight negative decline in richness with depth

# STEP 3: picking a suitable model ---------------------------------------------

# Framework:
# glm_obj <- 
#   glm(formula = {RESPONSE_VAR} ~ {PREDICTOR_VAR1} + {PREDICTOR_VAR2}, 
#       data = clean_data, 
#       family = "{CHOICE_OF_DISTRIBUTION}")

rich_depth_mod <-
  glm(formula = Richness ~ Depth,
      data = clean_data,
      family = "poisson") # dealing with count data: try possion

# Choosing the distribution family
"poisson" # count data (e.g. number of species), default: link = "log"
"quasipoisson" # similar to poisson but a bit more overdispersed
"gaussian" # normally distributed data, default: link = "identity" (no transform)
"binomial" # yes/no, true/false, 0/1 type data, default: link = "logit" (logistic)

rich_depth_mod %>% summary()
rich_depth_mod %>% performance::check_model() # looks good

# no significant affect of depth on site richness

# Interpretation of the model summary: -----------------------------------------

# t or Z value = Estimate divided by SE
# this corresponds to "how many SE's the coeff is from zero"
# Higher Z value = low error for estimate value = low pval = high significance
# More information at: https://www.statology.org/interpret-glm-output-in-r/

# what do the parameters actually mean
# for the poisson, the back transform function = exp()
exp(2.64983)
exp(2.64983 + 1*-0.06428) # depth = 1m
exp(2.64983 + 2*-0.06428) # depth = 2m

# Good info on interpreting poisson model output
# https://stats.stackexchange.com/questions/11096/how-to-interpret-coefficients-in-a-poisson-regression


# STEP 4: plotting stats on figure ---------------------------------------------

# creating a function to backtransform the data
rich_depth_mod_INVERSE <- rich_depth_mod$family$linkinv # back-transformation (= exp())

clean_data %>% 
  # adding in the model fit:
  mutate(pred_trans = predict(rich_depth_mod, newdata = clean_data, type = "link", se.fit=TRUE)$fit, 
         pred_se_trans = predict(rich_depth_mod, newdata = clean_data, type = "link", se.fit=TRUE)$se.fit, 
         pred = rich_depth_mod_INVERSE(pred_trans), 
         lo = rich_depth_mod_INVERSE(pred_trans - 1.96*pred_se_trans), 
         hi = rich_depth_mod_INVERSE(pred_trans + 1.96*pred_se_trans)) |> 
  ggplot() + 
  aes(
    x = Depth, 
    y = Richness
  ) +
  geom_point()  + 
  # stat_smooth(se = TRUE) +  # we don't need this anymore
  geom_line(aes(y = pred), col = "red") +
  geom_line(aes(y = lo), col = "red", lty = 2) +
  geom_line(aes(y = hi), col = "red", lty = 2) +
  labs(
    x = "Depth (m)",
    y = "Richness (# species)"
  ) +
  theme_bw()

# EXAMPLE 2 (Numerical vs categorical) =========================================

# STEP 1: visualise data -------------------------------------------------------

# Categorical vs numerical plot
clean_data %>% 
  summarise(richness_mean = mean(Richness, na.rm = TRUE), 
            .by = MPA) %>% # .by is used to define the grouping variable 
  ggplot() + 
  aes(
    x = MPA, 
    y = richness_mean
  ) +
  geom_col() +
  labs(
    x = "MPA status",
    y = "Richness (# species)"
  ) +
  theme_bw()

# STEP 2: interpret ------------------------------------------------------------
# I don't expect a significant difference here. Hard to tell without idea of spread

# STEP 3: picking a suitable model ---------------------------------------------

# Framework:
# glm_obj <- 
#   glm(formula = {RESPONSE_VAR} ~ {PREDICTOR_VAR1} + {PREDICTOR_VAR2}, 
#       data = clean_data, 
#       family = "{CHOICE_OF_DISTRIBUTION}")

rich_mpa_mod <-
  glm(formula = Richness ~ MPA,
      data = clean_data,
      family = "poisson") # dealing with count data: try possion

rich_mpa_mod %>% summary()
rich_mpa_mod %>% performance::check_model() # looks good

# significant affect of MPA on site richness

# STEP 4: plotting stats on figure ---------------------------------------------

# creating a function to backtransform the data
rich_mpa_mod_INVERSE <- rich_mpa_mod$family$linkinv # back-transformation (= exp())

clean_data %>% 
  summarise(richness_mean = mean(Richness, na.rm = TRUE), 
            .by = MPA) %>%
  mutate(pred_trans = predict(rich_mpa_mod, newdata = ., type = "link", se.fit=TRUE)$fit, 
         pred_se_trans = predict(rich_mpa_mod, newdata = ., type = "link", se.fit=TRUE)$se.fit, 
         pred = rich_mpa_mod_INVERSE(pred_trans), 
         lo = rich_mpa_mod_INVERSE(pred_trans - 1.96*pred_se_trans), 
         hi = rich_mpa_mod_INVERSE(pred_trans + 1.96*pred_se_trans)) |> 
  ggplot() + 
  aes(
    x = MPA, 
    y = richness_mean
  ) +
  geom_col() +
  # adding the appropirate errorbars
  geom_errorbar(aes(y = pred, 
                    ymin = lo, 
                    ymax = hi),
                width = 0.1) +
  labs(
    x = "MPA status",
    y = "Richness (# species)"
  ) +
  theme_bw()


# EXAMPLE 3 (Numerical vs categorical - less straightforward) ==================

# STEP 1: visualise data -------------------------------------------------------

# Categorical vs numerical plot
clean_data %>% 
  summarise(abundance_mean = mean(N, na.rm = TRUE), 
            .by = MPA) %>% # .by is used to define the grouping variable 
  ggplot() + 
  aes(
    x = MPA, 
    y = abundance_mean
  ) +
  geom_col() +
  labs(
    x = "MPA status",
    y = "Richness (# species)"
  ) +
  theme_bw()

# STEP 2: interpret ------------------------------------------------------------
# Reserve has lower abundance. Maybe a signficant difference?? 

# STEP 3: picking a suitable model ---------------------------------------------

# Framework:
# glm_obj <- 
#   glm(formula = {RESPONSE_VAR} ~ {PREDICTOR_VAR1} + {PREDICTOR_VAR2}, 
#       data = clean_data, 
#       family = "{CHOICE_OF_DISTRIBUTION}")

abun_mpa_mod <-
  glm(formula = N ~ MPA,
      data = clean_data,
      family = "poisson") # dealing with count data: try possion

abun_mpa_mod %>% summary()
abun_mpa_mod %>% performance::check_model() # not as good.

# lets try another model.
# alternatives to poisson include "quasipoission" and "negative binomial"
# these are just slight more flexible distributions (they can handle overdispersion better)
abun_mpa_mod_qp <-
  glm(formula = N ~ MPA,
      data = clean_data,
      family = "quasipoisson")

abun_mpa_mod_qp %>% summary()
abun_mpa_mod_qp %>% performance::check_model() # not that sure about it. maybe negbinom


# what about negative binomial?
# we need a slighly different glm function from the MASS package
abun_mpa_mod_nb <-
  MASS::glm.nb(formula = N ~ MPA,
      data = clean_data)

abun_mpa_mod_nb %>% summary()
abun_mpa_mod_nb %>% performance::check_model() # a bit better. let's go with this

# no significant affect of MPA on site abundance

# STEP 4: plotting stats on figure ---------------------------------------------

# creating a function to backtransform the data
abun_mpa_mod_nb_INVERSE <- abun_mpa_mod_nb$family$linkinv 

clean_data %>% 
  summarise(abundance_mean = mean(N, na.rm = TRUE), 
            .by = MPA) %>%
  mutate(pred_trans = predict(abun_mpa_mod_nb, newdata = ., type = "link", se.fit=TRUE)$fit, 
         pred_se_trans = predict(abun_mpa_mod_nb, newdata = ., type = "link", se.fit=TRUE)$se.fit, 
         pred = abun_mpa_mod_nb_INVERSE(pred_trans), 
         lo = abun_mpa_mod_nb_INVERSE(pred_trans - 1.96*pred_se_trans), 
         hi = abun_mpa_mod_nb_INVERSE(pred_trans + 1.96*pred_se_trans)) |> 
  ggplot() + 
  aes(
    x = MPA, 
    y = abundance_mean
  ) +
  geom_col() +
  # adding the appropriate errorbars
  geom_errorbar(aes(y = pred, 
                    ymin = lo, 
                    ymax = hi),
                width = 0.1) +
  labs(
    x = "MPA status",
    y = "Abundance"
  ) +
  theme_bw()
