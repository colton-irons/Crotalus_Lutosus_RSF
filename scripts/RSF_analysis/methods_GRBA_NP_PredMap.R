# load packages ----------------------
library(pacman)
p_load(dplyr, readxl, tidyterra, lubridate, terra, sf,
       tidyverse, mgcv, viridisLite, gratia)


## read in GRBA national park polygon as vector
grba_nps.v <- vect("/Users/coltonirons/Desktop/GRBA rattlesnake proj/grbanps.shp")


#read in DEM
master_DEM<- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/GRBA_Rasters/DEM_GRBANP.tif")
grba_aspect.rast <- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/GRBA_Rasters/aspectgrba_DEM.tif")

## read in ruggedness layer provided by Tom Dilts
rugged_grbanp.rast <- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/vrml_mosaic.tif")


### read in distance to water, roads, and hibernacula 
road_grbanp.rast <- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/GRBA_Rasters/grbanp_road.tif")
water_grbanp.rast <- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/GRBA_Rasters/grbanp_stream.tif")
hiber_grbanp.rast <- rast("/Users/coltonirons/Desktop/GRBA rattlesnake Proj/rasters/GRBA_Rasters/grbanp_den.tif")


# Crop and mask vegetation to the park boundary ----------------------
easternNV.veg.rast <- rast("/Users/coltonirons/Desktop/Veg/easternNV_veg.tif")

grba_nps_veg_crs <- terra::project(
  grba_nps.v,
  terra::crs(easternNV.veg.rast)
)

grbanp_veg <- easternNV.veg.rast |>
  terra::crop(grba_nps_veg_crs) |>
  terra::mask(grba_nps_veg_crs)

terra::plot(grbanp_veg)

# Calculate northwest-facing aspect ----------------------

northwest_aspect <- cos(
  grba_aspect.rast[["aspect"]] * pi / 180 -
    7 * pi / 4
)

names(northwest_aspect) <- "NW_aspect"

terra::plot(
  northwest_aspect,
  main = "Northwest-facing aspect"
)

# Align terrain and northwest aspect to the park vegetation grid -------------

terrain_grbanp <- terra::project(
  master_DEM,
  grbanp_veg,
  method = "cubicspline"
)

terrain_grbanp <- terra::mask(
  terrain_grbanp,
  grbanp_veg
)

northwest_aspect_grbanp <- terra::project(
  northwest_aspect,
  terrain_grbanp[["elevation"]],
  method = "bilinear"
)


northwest_aspect_grbanp <- terra::mask(
  northwest_aspect_grbanp,
  terrain_grbanp[["elevation"]]
)

# Combine aligned predictors
raster_stack_grbanp <- c(
  terrain_grbanp,
  northwest_aspect_grbanp,
  grbanp_veg
)

terra::plot(raster_stack_grbanp)

elevation_template <- raster_stack_grbanp[["elevation"]]

# Align distance rasters to the park-wide template ---------------------------

# incorporate distance to water, roads, and hibernacula 
bakewater_grbanp <- project(water_grbanp.rast, elevation_template, method = 'cubicspline')
bakewater_grbanp <- mask(bakewater_grbanp, elevation_template)
bakeroad_grbanp <- project(road_grbanp.rast, elevation_template, method = 'cubicspline')
bakeroad_grbanp <- mask(bakeroad_grbanp, elevation_template)
bakehiber_grbanp <- project(hiber_grbanp.rast,elevation_template, method = 'cubicspline')
bakehiber_grbanp <- mask(bakehiber_grbanp, elevation_template)

GRBA_dist2lyr.r <- c(bakewater_grbanp, bakeroad_grbanp, bakehiber_grbanp)


names(GRBA_dist2lyr.r)[1] <- "Dist2Water"
names(GRBA_dist2lyr.r)[2] <- "Dist2Roads"
names(GRBA_dist2lyr.r)[3] <- "Dist2Hiber"
plot(GRBA_dist2lyr.r)

names(GRBA_dist2lyr.r)

terra::compareGeom(
  elevation_template,
  GRBA_dist2lyr.r,
  stopOnError = FALSE
)

# Align ruggedness to the park-wide template ---------------------------------

ruggedness_grbanp <- rugged_grbanp.rast |>
  terra::project(
    elevation_template,
    method = "cubicspline"
  ) |>
  terra::mask(elevation_template)

names(ruggedness_grbanp) <- "ruggedness"

terra::plot(
  ruggedness_grbanp,
  main = "Terrain ruggedness"
)



# Calculate mean active-season land-surface temperature ----------------------

GRBAmastertemp_lst.rast_1 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2010ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_1)
GRBAmastertemp_lst.rast_1$mean <- GRBAmastertemp_lst.rast_1$climateEngine2010ActSeas.LST
GRBAmastertemp_lst.rast_2 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2011ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_2)
GRBAmastertemp_lst.rast_2$mean <- GRBAmastertemp_lst.rast_2$climateEngine2011ActSeas.LST
GRBAmastertemp_lst.rast_3 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2012ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_3)
GRBAmastertemp_lst.rast_3$mean <- GRBAmastertemp_lst.rast_3$climateEngine2012ActSeas.LST
GRBAmastertemp_lst.rast_4 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2013ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_4)
GRBAmastertemp_lst.rast_4$mean <- GRBAmastertemp_lst.rast_4$climateEngine2013ActSeas.LST
GRBAmastertemp_lst.rast_5 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2014ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_5)
GRBAmastertemp_lst.rast_5$mean <- GRBAmastertemp_lst.rast_5$climateEngine2014ActSeas.LST
GRBAmastertemp_lst.rast_6 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2023ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_6)
GRBAmastertemp_lst.rast_6$mean <- GRBAmastertemp_lst.rast_6$climateEngine2023ActSeas.LST
GRBAmastertemp_lst.rast_7 <- rast("/Users/coltonirons/Desktop/TempRock_raster/LST_years/ClimateEngine2024ActSeas.LST.tif")
plot(GRBAmastertemp_lst.rast_7)
GRBAmastertemp_lst.rast_7$mean <- GRBAmastertemp_lst.rast_7$climateEngine2024ActSeas.LST





GRBAprefix_lst <- "GRBAmastertemp_lst.rast_" 
GRBAvarlist_lst <- names(GRBAmastertemp_lst.rast_1$mean)
seq<- 1:7


avgLST <- list()

y=1;v=1
for(v in 1:length(GRBAvarlist_lst)){
  thisv <- GRBAvarlist_lst[v]
  thisl <- list()
  for(y in 1:length(seq)){
    thisy = seq[y]
    thisl[[y]] <- eval(parse(text= paste0(GRBAprefix_lst,thisy,"$",thisv) ))
  }
  thisl <- rast(thisl)
  avgLST[[thisv]] <- app(thisl,"mean")
}

GRBAavgLST <- rast(avgLST)

GRBAavgLST_utm <- project(
  GRBAavgLST,
  elevation_template
)


st_crs(GRBAavgLST_utm)
class(GRBAavgLST_utm)
plot(GRBAavgLST_utm)
res(GRBA_dist2lyr.r)
res(GRBAavgLST)


names(GRBAavgLST)[1] <- "LSTseason.r"

grbalst <- project(GRBAavgLST$LSTseason.r, terrain_grbanp$elevation, method = 'cubicspline')
grbalst <- mask(grbalst, terrain_grbanp$elevation)
plot(grbalst)



# Prepare rock-cover and LST predictors for park-wide prediction --------------
grbarock.rast <- rast("/Users/coltonirons/Desktop/TempRock_raster/rocklayer/BlockStRock.tif")
plot(grbarock.rast)

names(grbarock.rast)[1] <- "rock.r"

grbarocklyr <- project(grbarock.rast$rock.r, terrain_grbanp$elevation, method = 'bilinear')
grbarocklyr <- mask(grbarocklyr, terrain_grbanp$elevation)
plot(grbarocklyr)
grbarocklyr

GRBArocktemp.r <- c(grbalst, grbarocklyr)
plot(GRBArocktemp.r)




# Combine all raster layers into a single SpatRaster --------------
raster_stack_grbanp_final <- c(raster_stack_grbanp, GRBA_dist2lyr.r, ruggedness_grbanp,GRBArocktemp.r)
plot(raster_stack_grbanp_final)

### plotting tail end of raster stack 
# Get the number of layers
n <- nlyr(raster_stack_grbanp_final)

# Plot the last 6 layers
plot(raster_stack_grbanp_final[[ (n-5):n ]])

### reading in the full points for 25693 points
combined_data25693 <- read_csv("~/Desktop/GRBA rattlesnake Proj/testing_trainingData/25707PApts/combined_data.csv")
sum(is.na(combined_data25693))

## read in GRBA national park polygon as spatial feature
grba_nps.sf <- st_read("/Users/coltonirons/Desktop/GRBA rattlesnake proj/grbanps.shp")

n25000_background <- 25000   # set number of background points
total_area_ha <- grba_nps.sf$area_ha
bg_weight <- total_area_ha / n25000_background
1/bg_weight    # approx 15 background points per hectare (seems like plenty!) (each random point represents approximately 625 m2- on the same order of magnitude as 30*30m pixel...)

combined_data25693$Weights <- 1e-6
combined_data25693$Weights[combined_data25693$pres == 0] <- bg_weight

bg_points <- subset(combined_data25693,pres==0)

table(bg_points$YEAR,bg_points$PIT4)   # bg points are fairly well distributed by year and snake ID

combined_data25693$PIT4.f <- as.factor(combined_data25693$PIT4)
combined_data25693$YEAR.f <- as.factor(combined_data25693$YEAR)

table(combined_data25693$pres, combined_data25693$PIT4)





# Fit no-GP model for park-wide prediction --------------


# Confirm required model variables are properly formatted
stopifnot(
  is.factor(combined_data25693$YEAR.f),
  all(combined_data25693$Weights > 0)
)

park_model_formula <- pres / Weights ~
  s(elevation, bs = "cr") +
  s(LSTseason.r, bs = "cr") +
  s(roughness, bs = "cr") +
  s(Dist2Water, bs = "cr") +
  s(SHR, bs = "cr") +
  s(TRE, bs = "cr") +
  s(YEAR.f, bs = "re")

park_model_no_gp <- mgcv::gam(
  formula = park_model_formula,
  data = combined_data25693,
  family = poisson(link = "log"),
  weights = Weights,
  method = "REML"
)

# Model summary and diagnostics
summary(park_model_no_gp) #  Deviance explained = 4.82%
mgcv::gam.check(park_model_no_gp)

# Plot estimated smooth effects
gratia::draw(park_model_no_gp)

# Predict original study area with the no-GP model --------------


SA_raster_stack <- rast("/Users/coltonirons/Desktop/GitHub/Crotalus_Lutosus_RSF/scripts/RSF_analysis/bake_mastercovarFINALFINAL.tif")


study_area_data <- terra::as.data.frame(
  SA_raster_stack,
  xy = TRUE,
  na.rm = FALSE
)

# Confirm that all model predictors are available
required_predictors <- c(
  "elevation",
  "LSTseason.r",
  "roughness",
  "Dist2Water",
  "SHR",
  "TRE"
)

missing_predictors <- setdiff(
  required_predictors,
  names(study_area_data)
)

if (length(missing_predictors) > 0) {
  stop(
    "Missing study-area predictors: ",
    paste(missing_predictors, collapse = ", ")
  )
}

# Add a valid year factor required by the fitted model
study_area_data$YEAR.f <- factor(
  "2024",
  levels = levels(combined_data25693$YEAR.f)
)

# Predict on the link scale while excluding the year random effect
study_area_data$predicted_log_no_gp <- stats::predict(
  park_model_no_gp,
  newdata = study_area_data,
  type = "link",
  exclude = "s(YEAR.f)"
)

# Calculate study-area reference mean ----------------------------------------

study_mean_lp <- mean(
  study_area_data$predicted_log_no_gp,
  na.rm = TRUE
)

if (!is.finite(study_mean_lp)) {
  stop("The study-area mean prediction is not finite.")
}

# Center study-area predictions
study_area_data$centered_log_use <-
  study_area_data$predicted_log_no_gp - study_mean_lp

# Verify that centered predictions have a mean near zero
centered_mean <- mean(
  study_area_data$centered_log_use,
  na.rm = TRUE
)

print(study_mean_lp)
print(centered_mean)


# Prepare Great Basin National Park predictor data ----------------------------

park_data <- terra::as.data.frame(
  raster_stack_grbanp_final,
  xy = TRUE,
  na.rm = FALSE
)

required_predictors <- c(
  "elevation",
  "LSTseason.r",
  "roughness",
  "Dist2Water",
  "SHR",
  "TRE"
)

# Confirm that all required layers are present
missing_predictors <- setdiff(
  required_predictors,
  names(park_data)
)

if (length(missing_predictors) > 0) {
  stop(
    "Missing park-wide predictors: ",
    paste(missing_predictors, collapse = ", ")
  )
}

# Identify cells with complete predictor data
park_data$complete_predictors <- stats::complete.cases(
  park_data[required_predictors]
)

table(park_data$complete_predictors)

# Clamp park predictors to study-area training ranges ------------------------

model_predictors <- c(
  "elevation",
  "LSTseason.r",
  "roughness",
  "Dist2Water",
  "SHR",
  "TRE"
)

# Calculate training range for each predictor
training_ranges <- lapply(
  model_predictors,
  function(variable) {
    range(
      combined_data25693[[variable]],
      na.rm = TRUE
    )
  }
)

names(training_ranges) <- model_predictors

# Preserve the original park-wide values
park_prediction_data <- park_data

# Store information about clamping
clamping_summary <- data.frame(
  predictor = model_predictors,
  training_min = NA_real_,
  training_max = NA_real_,
  cells_below = NA_integer_,
  cells_above = NA_integer_
)

for (i in seq_along(model_predictors)) {
  
  variable <- model_predictors[i]
  limits <- training_ranges[[variable]]
  original_values <- park_prediction_data[[variable]]
  
  clamping_summary$training_min[i] <- limits[1]
  clamping_summary$training_max[i] <- limits[2]
  
  clamping_summary$cells_below[i] <- sum(
    original_values < limits[1],
    na.rm = TRUE
  )
  
  clamping_summary$cells_above[i] <- sum(
    original_values > limits[2],
    na.rm = TRUE
  )
  
  park_prediction_data[[variable]] <- pmin(
    pmax(original_values, limits[1]),
    limits[2]
  )
}

clamping_summary

# Add placeholder year factor required for prediction ------------------------

reference_year <- "2024"
year_levels <- levels(combined_data25693$YEAR.f)

if (!reference_year %in% year_levels) {
  stop("The reference year is not a level of YEAR.f.")
}

park_prediction_data$YEAR.f <- factor(
  rep(reference_year, nrow(park_prediction_data)),
  levels = year_levels
)

stopifnot(
  !anyNA(park_prediction_data$YEAR.f)
)


# Generate park-wide predictions with the no-GP model ------------------------

complete_cells <- stats::complete.cases(
  park_prediction_data[model_predictors]
)

park_prediction_data$predicted_log_no_gp <- NA_real_

park_prediction_data$predicted_log_no_gp[
  complete_cells
] <- stats::predict(
  park_model_no_gp,
  newdata = park_prediction_data[
    complete_cells,
  ],
  type = "link",
  exclude = "s(YEAR.f)"
)

# Center park-wide predictions using the study-area reference mean ------------

park_prediction_data$centered_log_use <-
  park_prediction_data$predicted_log_no_gp -
  study_mean_lp

# Summarize centered park predictions
park_prediction_summary <- summary(
  park_prediction_data$centered_log_use
)

park_centered_mean <- mean(
  park_prediction_data$centered_log_use,
  na.rm = TRUE
)

print(park_prediction_summary)
print(park_centered_mean)

# Create centered park-wide prediction raster --------------------------------

stopifnot(
  nrow(park_prediction_data) ==
    terra::ncell(raster_stack_grbanp_final)
)

park_prediction_raster <-
  raster_stack_grbanp_final[[1]]

terra::values(park_prediction_raster) <-
  park_prediction_data$centered_log_use

names(park_prediction_raster) <-
  "Relative_Resource_Use"

park_prediction_raster

terra::nlyr(park_prediction_raster)
names(park_prediction_raster)
terra::minmax(park_prediction_raster)

terra::global(
  park_prediction_raster,
  fun = c("mean", "min", "max"),
  na.rm = TRUE
)

terra::plot(
  park_prediction_raster,
  main = "Predicted Relative Resource Use",
  col = viridisLite::viridis(100)
)

#writeRaster(park_prediction_raster, "Final_RSF_ParkpredmapCentered.tif")

