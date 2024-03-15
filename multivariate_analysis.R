# About ========================================================================

# Title: Multivariate analysis for Maria Island Data (KSM308)
# Date: March 2024
# Author: Freddie Heather

# Imports data
# Visualises data  via ordination (nMDS)
# Relate environmental data to abundance counts (adonis)
# Generates figures

# Ecological distances: https://www.youtube.com/watch?v=xyufizOpc5I
# Great tutorial on nMDS: https://riffomonas.org/code_club/2022-02-14-nmds
# Another good tutorial on nMDS in R: https://www.rpubs.com/RGrieger/545184

# Important: 
# 1. Make sure you are working inside an R project
# 2. Make sure the datafiles "multivariate_analysis.R" and "predictor_table.csv" are within your project working directory

# Packages =====================================================================

# install.packages("tidyverse") # run this line of code if you haven't installed the tidyverse on your computer
# install.packages("vegan")
# install.packges("ggrepel")
library(tidyverse)
library(vegan)
library(ggrepel)

# Import =======================================================================

resp_vars <- read_csv("multivariate_response_vars.csv") 
pred_vars <- read_csv("predictor_table.csv")

# Data checking ================================================================

resp_vars # has NA instead of zeros
pred_vars # has some blank columns

# Check to make sure the columns are of the right class (dbl or chr)

# Data wrangling ===============================================================

# mutate_at: replace SOME cols with NA values with zero
resp_vars_clean <- 
  resp_vars %>% 
  mutate_at(.vars = vars(-sample), # everything except 'sample' col
            .funs = replace_na, 
            replace = 0) 

# mutate_all: replace ALL the NA values with zero
resp_vars_clean <- 
  resp_vars %>% 
  mutate_all(.funs = replace_na, replace = 0) 

pred_vars_clean <- 
  pred_vars %>% 
  discard(~all(is.na(.))) %>%  # gets rid of those blank cols
  filter(!is.na(sample))

# Data visualisation ===========================================================

abundance_matrix <- 
  resp_vars_clean %>% 
  # mutate_if(is.numeric, log1p) %>% # log(x+1) transformation
  column_to_rownames("sample") %>% 
  as.matrix() 

# Calculating distance
# Bray curtis dissimilarity: https://www.statisticshowto.com/bray-curtis-dissimilarity/
distance_matrix <- 
  abundance_matrix %>% 
  vegdist(method = "bray")

set.seed(1)
nmds <- metaMDS(distance_matrix)
# stress = 0.1519975 

nmds_vals <- 
  scores(nMDS) %>% 
  as_tibble(rownames = "sample") %>% 
  left_join(pred_vars_clean, by = join_by(sample))

elps_center <- 
  nmds_vals %>% 
  summarise(NMDS1 = mean(NMDS1), 
            NMDS2 = mean(NMDS2), 
            .by = MPA)

# play with colouring points to identify drivers
nmds_vals %>% 
  ggplot() +
  aes(x = NMDS1, 
      y = NMDS2, 
      col = MPA, # change colour of points
  ) +
  geom_point(
    # aes(pch = as.character(Depth)) # change shape of point
  ) +
  stat_ellipse() +
  geom_point(data = elps_center, size = 5)

# Detecting statistical difference =============================================

# Using the adonis function

# Great tutorial on adonis: https://www.youtube.com/watch?v=1ETBgbXl-BM
# Note the use of adonis2 here (adonis used in the above video is depreciated).

# Response var NEEDS to be a distance matrix
class(distance_matrix) # class == "dist"
adonis2(formula = distance_matrix ~ nmds_vals$MPA) # no sig effect of MPA

# if we have two predictor vars
adonis2(formula = distance_matrix ~ nmds_vals$MPA + nmds_vals$vis) 


# Adding environemtnal vars or significant species =============================

pred_vars_clean %>% 
  column_to_rownames("sample")

env <- envfit(nmds, pred_vars_clean %>% column_to_rownames("sample"), na.rm = TRUE)
env <- envfit(nmds, resp_vars_clean %>% column_to_rownames("sample"), na.rm = TRUE)
plot(nmds)
plot(env)

vectors <- 
scores(env, display = "vectors") %>% 
  as_tibble(rownames = "species") %>% 
  mutate(pval = env$vectors$pvals) %>% 
  filter(pval < 0.05)


nmds_vals %>% 
  ggplot() +
  aes(x = NMDS1, 
      y = NMDS2
  ) +
  geom_point(aes(col = MPA)) +
  geom_segment(
    aes(
    x = 0, 
    y = 0, 
    xend = NMDS1,
    yend = NMDS2
  ),
  data = sig_vectors, 
  arrow = arrow(length = unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = sig_vectors,
                           aes(label = species), 
                           cex = 5, 
                           direction = "both")
