# About ========================================================================

# Title: Linear modelling for Maria Island Data (KSM308)
# Date: March 2024
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

# Import =======================================================================

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
clean_data <- 
  resp_vars_clean %>% 
  left_join(pred_vars_clean, by = "Sample")

# Data visualisation ===========================================================

# Numerical vs numerical plot
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

# Categorical vs numerical plot
clean_data %>% 
  summarise(richness_mean = mean(Richness, na.rm = TRUE), 
            richness_sd = sd(Richness, na.rm = TRUE), 
            n = n(),
            .by = MPA) %>% 
  mutate(richness_se = richness_sd/sqrt(n)) %>% 
  ggplot() + 
  aes(
    x = MPA, 
    y = richness_mean
  ) +
  geom_col() +
  geom_errorbar(aes(ymin = richness_mean - richness_se, 
                    ymax = richness_mean + richness_se), 
                width = 0.2) +
  labs(
    x = "Depth (m)",
    y = "Richness (# species)"
  ) +
  theme_bw() # themes from https://ggplot2.tidyverse.org/reference/ggtheme.html

# Reserve vs abundance
clean_data %>% 
  summarise(abundance_mean = mean(N, na.rm = TRUE), 
            abundance_sd = sd(N, na.rm = TRUE), 
            n = n(),
            .by = MPA) %>% 
  mutate(abundance_se = abundance_sd/sqrt(n),
         ci = abundance_se*qt((1-alpha)/2 + .5, n-1)) %>% 
  ggplot() + 
  aes(
    x = MPA, 
    y = abundance_mean
  ) +
  geom_col() +
  geom_errorbar(aes(ymin = abundance_mean - abundance_se, 
                    ymax = abundance_mean + abundance_se), 
                width = 0.2) +
  labs(
    x = "Depth (m)",
    y = "Abundance (# individuals)"
  ) +
  theme_bw() # themes from https://ggplot2.tidyverse.org/reference/ggtheme.html

# Reserve vs large fish abundance
clean_data %>% 
  summarise(abundance_mean = mean(`N>25`, na.rm = TRUE), 
            abundance_sd = sd(`N>25`, na.rm = TRUE), 
            n = n(),
            .by = MPA) %>% 
  mutate(abundance_se = abundance_sd/sqrt(n)) %>% 
  ggplot() + 
  aes(
    x = MPA, 
    y = abundance_mean
  ) +
  geom_col() +
  geom_errorbar(aes(ymin = abundance_mean - abundance_se, 
                    ymax = abundance_mean + abundance_se), 
                width = 0.2) +
  labs(
    x = "Depth (m)",
    y = "Large Fish Abundance (# individuals)"
  ) +
  theme_bw() # themes from https://ggplot2.tidyverse.org/reference/ggtheme.html

# Linear modelling =============================================================

glm1 <- 
  glm(formula = Richness ~ MPA, 
      data = clean_data, 
      family = "poisson") # poisson as predictor variable is count data

summary(glm1)
# t or Z value = Estimate divided by SE
# Higher Z value = low error for estimate value = low pval = high significance
# More information at: https://www.statology.org/interpret-glm-output-in-r/


# Say we wanted another variable in the model
glm2 <- 
  glm(formula = Richness ~ MPA + Depth, 
      data = clean_data, 
      family = "poisson") 

summary(glm2)

# Model for abundance
glm3 <- 
  glm(formula = N ~ MPA, 
      data = clean_data, 
      family = "poisson") # poisson as predictor variable is count data

summary(glm3)

# Good info on interpreting poisson model output
# https://stats.stackexchange.com/questions/11096/how-to-interpret-coefficients-in-a-poisson-regression

exp(5.61)
exp(5.61 + 1*-0.25056) # reserve = TRUE
exp(5.61 + 0*-0.25056) # reserve = FALSE


# Choosing the family
"poisson" # count data (e.g. number of species)
"gaussian" # normally distributed data
"binomial" # yes/no, true/false, 0/1 type data


# Numerical predictor
glm4 <- 
  glm(formula = N ~ Depth, 
      data = clean_data, 
      family = "poisson") # poisson as predictor variable is count data

summary(glm4)

# Figures ======================================================================

# Calculation of confidence itervals
# If normally distribution (guassian) its easy
# Otherwise a bit more complex
# https://stackoverflow.com/questions/40985366/prediction-of-poisson-regression

# Normal: CI = mean + 1.96*(sd/sqrt(n))
# Poission: CI = inverse_link_function(mean + 1.96*(sd/sqrt(n)))


# Categorical
glm3_inverse_func <- glm3$family$linkinv
mod_out <- 
  tibble(MPA = c("Reserve", "Fished")) %>% 
  mutate(pred_trans = predict(glm3, newdata = list(MPA = MPA), type = "link", se.fit=TRUE)$fit, 
         pred_se_trans = predict(glm3, newdata = list(MPA = MPA), type = "link", se.fit=TRUE)$se.fit, 
         pred = glm3_inverse_func(pred_trans), 
         lo = glm3_inverse_func(pred_trans - 1.96*pred_se_trans), 
         hi = glm3_inverse_func(pred_trans + 1.96*pred_se_trans))

clean_data %>% 
  summarise(abundance_mean = mean(N, na.rm = TRUE), 
            abundance_sd = sd(N, na.rm = TRUE), 
            n = n(),
            .by = MPA) %>%
  ggplot() + 
  aes(
    x = MPA, 
    y = abundance_mean
  ) +
  geom_col() +
  geom_errorbar(aes(y = pred, 
                    ymin = lo, 
                    ymax = hi), 
                data = mod_out, 
                width = 0.1) +
  labs(
    x = "Depth (m)",
    y = "Large Fish Abundance (# individuals)"
  ) +
  theme_bw()

# Good article on Error bars
# https://www.nature.com/articles/nmeth.2659


# Two numerical values
clean_data %>% 
  mutate(pred_n = exp(predict(glm4))) %>% 
  ggplot() + 
  aes(
    x = Depth, 
    y = N
  ) +
  geom_point()  + 
  geom_line(aes(y = pred_n), col = "red") +
  stat_smooth(se = TRUE) + # smoothed trend line? standard error?
  labs(
    x = "Visibility (m)",
    y = "Abundance (# individuals)"
  ) +
  theme_bw() # themes from https://ggplot2.tidyverse.org/reference/ggtheme.html

glm4_inverse_func <- glm4$family$linkinv

clean_data %>% 
  mutate(pred_trans = predict(glm4, newdata = clean_data, type = "link", se.fit=TRUE)$fit, 
         pred_se_trans = predict(glm4, newdata = clean_data, type = "link", se.fit=TRUE)$se.fit, 
         pred = glm4_inverse_func(pred_trans), 
         lo = glm4_inverse_func(pred_trans - 1.96*pred_se_trans), 
         hi = glm4_inverse_func(pred_trans + 1.96*pred_se_trans)) %>% 
  ggplot() + 
  aes(
    x = Depth, 
    y = N
  ) +
  geom_point()  + 
  geom_line(aes(y = pred), col = "red") +
  geom_line(aes(y = lo), col = "red", lty = 2) +
  geom_line(aes(y = hi), col = "red", lty = 2) +
  # stat_smooth(se = TRUE) + # smoothed trend line? standard error?
  labs(
    x = "Depth (m)",
    y = "Abundance (# individuals)"
  ) +
  theme_bw() 

