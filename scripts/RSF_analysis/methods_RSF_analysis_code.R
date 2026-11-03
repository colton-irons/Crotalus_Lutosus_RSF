# load packages ----------------------
library(pacman)
p_load(dplyr, tidyterra, lubridate, terra, sf, cowplot,
        raster, tidyverse, mgcv, viridisLite, gratia, ggcorrplot)



# load spatial and analysis data --------------

#don't include this dataset in GitHub
snake_usedpts <- st_read("/Users/coltonirons/Desktop/GRBA rattlesnake proj/thin data/snakeusedpts.utm.sf.gpkg")

#maybe include in GitHub?
snake_usedptsDATE <- st_read("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/thin data/thinned_snakeptsNO2009.shp")

### reading in all the data points combining background and used locations 
combined_data25693pts <- read_csv("~/Desktop/GRBA rattlesnake Proj/testing_trainingData/25707PApts/combined_data.csv")

raster_stack <- rast("/Users/coltonirons/Desktop/GitHub/Crotalus_Lutosus_RSF/scripts/RSF_analysis/bake_mastercovarFINALFINAL.tif")

study_area <- vect("~/Desktop/GRBA rattlesnake Proj/Rattlesnake scripts /bakeSA_740m.shp")


# check data --------------

table(snake_usedpts$YEAR)
table(snake_usedpts$YEAR,snake_usedpts$PIT4)

# Standardize four-digit PIT identifiers -------------------------------------

format_pit4 <- function(x) {
  stringr::str_pad(
    stringr::str_trim(as.character(x)),
    width = 4,
    side = "left",
    pad = "0"
  )
}


combined_data25693pts$PIT4 <-
  format_pit4(combined_data25693pts$PIT4)

# Plot spatial inputs --------------

plot(raster_stack[["elevation"]])
plot(study_area,add=T)
plot(snake_usedpts$geom,pch=20,col="gray",add=T)


# Assess correlation among broader predictors before deciding on final model -----------------

candidate_predictors <- combined_data25693pts |>
  dplyr::select(
    elevation,
    roughness,
    SHR,
    LTR,
    TRE,
    Dist2Water,
    Dist2Hiber,
    LSTseason.r
  ) |>
  tidyr::drop_na()

# Pearson correlation matrix
cor_matrix <- stats::cor(
  candidate_predictors,
  method = "pearson"
)

# Variance inflation factors
vif_results <- data.frame(
  predictor = names(candidate_predictors),
  VIF = diag(solve(cor_matrix)),
  row.names = NULL
)

vif_results

# Plot correlation matrix
correlation_plot <- ggcorrplot::ggcorrplot(
  cor_matrix,
  lab = TRUE,
  lab_size = 3,
  colors = c("blue", "white", "red"),
  hc.order = TRUE,
  type = "upper"
)

correlation_plot

# Examine correlation between LST and elevation
stats::cor(
  candidate_predictors$LSTseason.r,
  candidate_predictors$elevation
)

# prepare data for GAM ---------------


## prepare weights for downweighted poisson regression approach ------------

n25000_background <- 25000   # set number of background points
total_area_ha <- study_area$area_ha
bg_weight <- total_area_ha / n25000_background
1/bg_weight    # approx 15 background points per hectare (seems like plenty!) (each random point represents approximately 625 m2- on the same order of magnitude as 30*30m pixel...)

combined_data25693pts$Weights <- 1e-6
combined_data25693pts$Weights[combined_data25693pts$pres == 0] <- bg_weight

bg_points <- subset(combined_data25693pts,pres==0)


combined_data25693pts$PIT4.f <- as.factor(combined_data25693pts$PIT4)
combined_data25693pts$YEAR.f <- as.factor(combined_data25693pts$YEAR)

table(combined_data25693pts$pres, combined_data25693pts$PIT4)



# Center and scale model predictors ---------------

standardize <- function(x) {
  as.numeric(scale(x))
}

combined_data25693pts <- combined_data25693pts |>
  dplyr::mutate(
    elevation_sc  = standardize(elevation),
    LSTseason_sc  = standardize(LSTseason.r),
    roughness_sc  = standardize(roughness),
    Dist2Water_sc = standardize(Dist2Water),
    SHR_sc        = standardize(SHR),
    TRE_sc        = standardize(TRE)
  )

# Save means and standard deviations for back-transformation
mean_elev <- mean(combined_data25693pts$elevation, na.rm = TRUE)
sd_elev   <- sd(combined_data25693pts$elevation, na.rm = TRUE)

mean_LST <- mean(combined_data25693pts$LSTseason.r, na.rm = TRUE)
sd_LST   <- sd(combined_data25693pts$LSTseason.r, na.rm = TRUE)

mean_rough <- mean(combined_data25693pts$roughness, na.rm = TRUE)
sd_rough   <- sd(combined_data25693pts$roughness, na.rm = TRUE)

mean_Dist2Water <- mean(combined_data25693pts$Dist2Water, na.rm = TRUE)
sd_Dist2Water <- sd(combined_data25693pts$Dist2Water, na.rm = TRUE)

mean_SHR <- mean(combined_data25693pts$SHR, na.rm = TRUE)
sd_SHR   <- sd(combined_data25693pts$SHR, na.rm = TRUE)

mean_TRE <- mean(combined_data25693pts$TRE, na.rm = TRUE)
sd_TRE   <- sd(combined_data25693pts$TRE, na.rm = TRUE)

#check 
hist(combined_data25693pts$elevation)
hist(combined_data25693pts$elevation_sc)

hist(combined_data25693pts$LSTseason.r)
hist(combined_data25693pts$LSTseason_sc)

hist(combined_data25693pts$roughness)
hist(combined_data25693pts$roughness_sc)

hist(combined_data25693pts$Dist2Water)
hist(combined_data25693pts$Dist2Water_sc)

hist(combined_data25693pts$SHR)
hist(combined_data25693pts$SHR_sc)

hist(combined_data25693pts$TRE)
hist(combined_data25693pts$TRE_sc)

# run final GAM model ----------------------
set.seed(939)
gam_final <- gam(
  pres / Weights  ~ 
    ## Global smooths 
    s(elevation,bs="cr") +
    s(LSTseason.r,bs="cr") +
    s(roughness,bs="cr") +
    s(Dist2Water,bs="cr") +
    s(SHR,bs="cr") +
    s(TRE,bs="cr") +
    s(YEAR.f, bs = 're')+ # random intercept for year
    s(Easting, Northing, k = 200, bs = "gp"),  # Spatial random effect  
  data = combined_data25693pts,
  method="REML",
  family = poisson(link = "log"),
  weights = Weights
)

summary(gam_final) ## Dev.Expl =  10.8%
gam.check(gam_final)
draw(gam_final)  

#saveRDS(gam_final,"GAMfinalmod.RDS")

##### Supplemental Materials  -------------- 

# visualize partial effects of spatial random effect and fixed effect on year

sms <- smooths(gam_final)
#draw(gam_final, select = sms[1], rug = FALSE) + ggtitle("") visualize elevation

draw(
  gam_final, 
  select = sms[8], 
  rug = FALSE, 
  ylab = "Northing", 
  xlab = "Easting"
) + ggtitle("")

draw(
  gam_final, 
  select = sms[7], 
) + ggtitle("")


# run final GAM model scaled and centered  ----------------------
set.seed(935)
gam_finalmod_sc <- gam(
  pres / Weights  ~ 
    ## Global smooths 
    s(elevation_sc,bs="cr") +
    s(LSTseason_sc,bs="cr") +
    s(roughness_sc,bs="cr") +
    s(Dist2Water_sc,bs="cr") +
    s(SHR_sc,bs="cr") +
    s(TRE_sc,bs="cr") +
    s(YEAR.f, bs = 're')+ # random intercept for year
    s(Easting, Northing, k = 200, bs = "gp"),  # Spatial random effect  
  data = combined_data25693pts,
  method="REML",
  family = poisson(link = "log"),
  weights = Weights
)

summary(gam_finalmod_sc) ## Dev.Expl =  10.8%
gam.check(gam_finalmod_sc)
draw(gam_finalmod_sc)  

#saveRDS(gam_finalmod_sc,"GAMfinalmod_sc.RDS")

# Function to calculate and plot log-RSS ---------------

make_log_rss_plot <- function(
    model,
    model_data,
    predictor_sc,
    predictor_raw,
    predictor_mean,
    predictor_sd,
    rug_data,
    x_label,
    n = 200
) {
  
  # Sequence covering the observed standardized predictor range
  predictor_sequence <- seq(
    min(model_data[[predictor_sc]], na.rm = TRUE),
    max(model_data[[predictor_sc]], na.rm = TRUE),
    length.out = n
  )
  
  # Hold all other continuous predictors at their means (standardized zero)
  newdata <- data.frame(
    elevation_sc = 0,
    LSTseason_sc = 0,
    roughness_sc = 0,
    Dist2Water_sc = 0,
    SHR_sc = 0,
    TRE_sc = 0,
    YEAR.f = factor(
      levels(model_data$YEAR.f)[1],
      levels = levels(model_data$YEAR.f)
    ),
    Easting = mean(model_data$Easting, na.rm = TRUE),
    Northing = mean(model_data$Northing, na.rm = TRUE)
  )
  
  # Repeat the reference row across the predictor sequence
  newdata <- newdata[rep(1, n), ]
  newdata[[predictor_sc]] <- predictor_sequence
  
  # Reference condition: predictor at its mean
  refdata <- newdata
  refdata[[predictor_sc]] <- 0
  
  # Linear-predictor matrices
  X_predictor <- stats::predict(
    model,
    newdata = newdata,
    type = "lpmatrix"
  )
  
  X_reference <- stats::predict(
    model,
    newdata = refdata,
    type = "lpmatrix"
  )
  
  X_difference <- X_predictor - X_reference
  
  # Calculate log-RSS and uncertainty
  coefficients <- stats::coef(model)
  covariance_matrix <- stats::vcov(model)
  
  results <- data.frame(
    predictor_sc = predictor_sequence,
    predictor_raw = (
      predictor_sequence * predictor_sd + predictor_mean
    ),
    log_RSS = drop(X_difference %*% coefficients)
  )
  
  results$log_RSS_se <- sqrt(
    rowSums(
      (X_difference %*% covariance_matrix) *
        X_difference
    )
  )
  
  results$log_RSS_lb <- results$log_RSS -
    1.96 * results$log_RSS_se
  
  results$log_RSS_ub <- results$log_RSS +
    1.96 * results$log_RSS_se
  
  # Construct figure
  plot <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = predictor_raw, y = log_RSS)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = log_RSS_lb,
        ymax = log_RSS_ub
      ),
      fill = "gray",
      alpha = 0.6
    ) +
    ggplot2::geom_line(
      linewidth = 1.5,
      color = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    ggplot2::geom_rug(
      data = rug_data,
      ggplot2::aes(x = .data[[predictor_raw]]),
      inherit.aes = FALSE,
      sides = "b"
    ) +
    ggplot2::labs(
      x = x_label,
      y = "log-RSS"
    ) +
    ggplot2::theme_minimal()
  
  list(
    data = results,
    plot = plot
  )
}

# Elevation ---------------

elevation_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "elevation_sc",
  predictor_raw = "elevation",
  predictor_mean = mean_elev,
  predictor_sd = sd_elev,
  rug_data = snake_usedpts,
  x_label = "Elevation (m)"
)

plot_a <- elevation_rss$plot
plot_a

# Land-surface temperature ---------------

lst_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "LSTseason_sc",
  predictor_raw = "LSTseason.r",
  predictor_mean = mean_LST,
  predictor_sd = sd_LST,
  rug_data = snake_usedpts,
  x_label = "Land Surface Temperature (°C)"
)

plot_b <- lst_rss$plot
plot_b

# roughness ---------------

rough_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "roughness_sc",
  predictor_raw = "roughness",
  predictor_mean = mean_rough,
  predictor_sd = sd_rough,
  rug_data = snake_usedpts,
  x_label = "Roughness (m)"
)

plot_c <- rough_rss$plot
plot_c

# Distance to water ---------------

D2water_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "Dist2Water_sc",
  predictor_raw = "Dist2Water",
  predictor_mean = mean_Dist2Water,
  predictor_sd = sd_Dist2Water,
  rug_data = snake_usedpts,
  x_label = "Distance to Water (m)"
)

plot_d <- D2water_rss$plot
plot_d

# shrub ---------------

shrub_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "SHR_sc",
  predictor_raw = "SHR",
  predictor_mean = mean_SHR,
  predictor_sd = sd_SHR,
  rug_data = snake_usedpts,
  x_label = "Shrub Cover (%)"
)

plot_e <- shrub_rss$plot
plot_e

# tree ---------------

tree_rss <- make_log_rss_plot(
  model = gam_finalmod_sc,
  model_data = combined_data25693pts,
  predictor_sc = "TRE_sc",
  predictor_raw = "TRE",
  predictor_mean = mean_TRE,
  predictor_sd = sd_TRE,
  rug_data = snake_usedpts,
  x_label = "Tree Cover (%)"
)

plot_f <- tree_rss$plot
plot_f


library(patchwork)

# create multi-panel figure with response curves ---------------

multi_plot <- wrap_plots(
  plot_a, plot_b,
  plot_c, plot_d,
  plot_e, plot_f,
  ncol = 2
) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold")
  )

multi_plot

#ggsave(
  #"RSF_logRSS_curves.png",
  #multi_plot,
  #width = 10,
  #height = 12,
  #units = "in",
  #dpi = 600
#)



# Snake-level grouped cross-validation ---------------

set.seed(321)

n_folds <- 10

# Model specification
rsf_formula <- pres / Weights ~
  s(elevation, bs = "cr") +
  s(LSTseason.r, bs = "cr") +
  s(roughness, bs = "cr") +
  s(Dist2Water, bs = "cr") +
  s(SHR, bs = "cr") +
  s(TRE, bs = "cr") +
  s(YEAR.f, bs = "re") +
  s(Easting, Northing, k = 200, bs = "gp")

# Fit the full model
final_mod_fast <- mgcv::bam(
  formula = rsf_formula,
  data = combined_data25693pts,
  family = poisson(link = "log"),
  weights = Weights,
  method = "fREML",
  discrete = TRUE,
  nthreads = 1
)

summary(final_mod_fast)
gratia::draw(final_mod_fast)

# Assign entire snakes to folds
snake_ids <- unique(stats::na.omit(snake_usedpts$PIT4))
snake_ids <- sample(snake_ids)

snake_folds <- split(
  snake_ids,
  rep(seq_len(n_folds), length.out = length(snake_ids))
)

# Confirm the number of snakes held out in each fold
lengths(snake_folds)

# Store full-model and cross-validated predictions
combined_data25693pts$predicted_full <- stats::predict(
  final_mod_fast,
  type = "response"
)

combined_data25693pts$predicted_CV <- NA_real_

# Snake-level cross-validation
for (fold in seq_along(snake_folds)) {
  
  test_snakes <- snake_folds[[fold]]
  
  train_index <- !combined_data25693pts$PIT4 %in% test_snakes
  test_index  <- combined_data25693pts$PIT4 %in% test_snakes
  
  train_data <- combined_data25693pts[train_index, ]
  test_data  <- combined_data25693pts[test_index, ]
  
  fold_model <- update(
    final_mod_fast,
    data = train_data
  )
  
  combined_data25693pts$predicted_CV[test_index] <-
    stats::predict(
      fold_model,
      newdata = test_data,
      type = "response"
    )
}

# Confirm that every observation received a CV prediction
sum(is.na(combined_data25693pts$predicted_CV))

# Evaluate discrimination ---------------

full_roc <- pROC::roc(
  response = combined_data25693pts$pres,
  predictor = combined_data25693pts$predicted_full,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

cv_roc <- pROC::roc(
  response = combined_data25693pts$pres,
  predictor = combined_data25693pts$predicted_CV,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

full_auc <- as.numeric(pROC::auc(full_roc))
cv_auc   <- as.numeric(pROC::auc(cv_roc))

full_auc #0.9111
cv_auc #0.8235

plot(cv_roc)

# Continuous Boyce index ---------------

calculate_boyce <- function(predictions, presence) {
  
  valid <- !is.na(predictions) & !is.na(presence)
  
  ecospat::ecospat.boyce(
    fit = predictions[valid],
    obs = predictions[valid & presence == 1]
  )
}

boyce_results <- list(
  full = calculate_boyce(
    predictions = combined_data25693pts$predicted_full,
    presence = combined_data25693pts$pres
  ),
  cross_validated = calculate_boyce(
    predictions = combined_data25693pts$predicted_CV,
    presence = combined_data25693pts$pres
  )
)

# Extract Boyce-index values
boyce_values <- vapply(
  boyce_results,
  FUN = function(x) x$cor,
  FUN.VALUE = numeric(1)
)

boyce_values

plot(
  boyce_results$cross_validated$HS,
  boyce_results$cross_validated$F.ratio,
  type = "l",
  xlab = "Predicted relative use",
  ylab = "Predicted-to-expected ratio"
)

# Create study-area prediction map ---------------

# Convert raster predictors to a data frame while retaining all cells
prediction_data <- terra::as.data.frame(
  raster_stack,
  xy = TRUE,
  na.rm = FALSE
) |>
  dplyr::mutate(
    Easting = x,
    Northing = y,
    YEAR.f = factor(
      "2024",
      levels = levels(combined_data25693pts$YEAR.f)
    )
  )

# Predict on the link scale, excluding the year random effect
prediction_data$predicted_log <- stats::predict(
  final_mod_fast,
  newdata = prediction_data,
  type = "link",
  exclude = "s(YEAR.f)"
)

# Place predictions into a raster with the original geometry
prediction_raster <- raster_stack[[1]]

terra::values(prediction_raster) <-
  prediction_data$predicted_log

names(prediction_raster) <- "predicted_log"

# Center predictions on the study-area mean
mean_lp <- terra::global(
  prediction_raster,
  fun = "mean",
  na.rm = TRUE
)[1, 1]

prediction_raster_centered <- prediction_raster - mean_lp

names(prediction_raster_centered) <-
  "centered_log_resource_use"

mean_lp

terra::global(
  prediction_raster_centered,
  fun = "mean",
  na.rm = TRUE
)

# Prepare data for plotting ---------------

prediction_df <- terra::as.data.frame(
  prediction_raster_centered,
  xy = TRUE,
  na.rm = TRUE
)

names(prediction_df)[3] <- "centered_lp"

# Ensure used locations match the raster CRS
used_sf <- sf::st_transform(
  snake_usedpts,
  crs = terra::crs(prediction_raster_centered)
)

prediction_range <- range(
  prediction_df$centered_lp,
  na.rm = TRUE
)

terra::plot(
  prediction_raster_centered,
  main = "Predicted Relative Resource Use",
  col = viridisLite::viridis(100)
)

terra::points(
  terra::vect(used_sf),
  pch = 20,
  col = "red",
  cex = 0.5
)


#writeRaster(prediction_raster_centered, "Final_RSF_Predictive_Map.tif")



