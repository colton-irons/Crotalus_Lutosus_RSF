# 2. Load the package
library(devtools)
library(rfUtilities)
library(dplyr)
library(randomForest)

# 1. Install rfUtilities from CRAN
install.packages("~/Desktop/rflab/rfUtilities_2.1-5.tar.gz")


ALL25693pts <- read.csv("/Users/coltonirons/Desktop/GitHub/Crotalus_Lutosus_RSF/Data/usedVSavailablepts_variables.csv")


raster_stackFINAL <- rast("/Users/coltonirons/Desktop/GitHub/Crotalus_Lutosus_RSF/Data/SA_mastercovarRastLyr.tif")

ALL25693pts.test <- ALL25693pts %>% 
  dplyr::select(  -PIT4, -YEAR, -Easting, -Northing)


b=1000
(dim(ALL25693pts.test[ALL25693pts.test$pres == 1,])[1] / dim(ALL25693pts.test)[1]) * 100
ALL25693pts.test$pres <- as.factor(ALL25693pts.test$pres)
rf.model <- rf.modelSel(x=ALL25693pts.test[,2:ncol(ALL25693pts.test)], y=ALL25693pts.test[,"pres"], imp.scale="mir", ntree=b)








