# R CODE FOR COMPILING THE OUTPUT OF THE ANALYSIS IN XXXX.R AND PUBLISHING IT IN A RASTER MAP.
# Author: Gabriel E. García Peña

# READ DATA (output of XXXXX.R)

# HAZARD ON THE SUSTAINABLE PATH IN THRE YEARS: 2025,2050,AND 2100
cell.sum.SSP1.2025 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_2025.csv") # SSP1: Sustainable development 2025
cell.sum.SSP1.2050 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_2050.csv") # SSP1: Sustainable development 2050
cell.sum.SSP1.2100 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_2100.csv") # SSP1: Sustainable development 2100

# HAZARD ON THE SHARED SOCIOECONOMIC PATHWAYS, FOR THRE YEAR 2050
cell.sum.SSP3.2050 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_SSP3_2050.csv") # SSP3: Region rivalry
cell.sum.SSP4.2050 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_SSP4_2050.csv") # SSP4: Inequality
cell.sum.SSP5.2050 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_worst_2050.csv") # SSP5: fossil-fuelled development

# HAZARD ON THE SHARED SOCIOECONOMIC PATHWAYS, FOR THRE YEAR 2100
cell.sum.SSP3.2100 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_SSP3_2100.csv") # SSP3: Region rivalry
cell.sum.SSP4.2100 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_SSP4_2100.csv") # SSP4: Inequality
cell.sum.SSP5.2100 = read.csv("~/LUCIDA/output/cell_analysis/cell_sum_N10_SSP5_2100.csv") # SSP5: fossil-fuelled development

# DATA ON RODENT RICHNES AT EACH GEOGRAPHIC CELL
richness = cell.sum.SSP4.2050$reservoir.richness
names(richness) = cell.sum.SSP4.2050$cell

# ASSIGN VARIABLES FOR THE YEAR 2050
SSP1.2025 <- cell.sum.SSP1.2025$cell.sum.reservoirs # BASELINE

SSP1.2050 <- cell.sum.SSP1.2050$cell.sum.reservoirs
SSP3.2050 <- cell.sum.SSP3.2050$cell.sum.reservoirs
SSP4.2050 <- cell.sum.SSP4.2050$cell.sum.reservoirs
SSP5.2050 <- cell.sum.SSP5.2050$cell.sum.reservoirs

# ASSIGN NAMES cell
names(SSP1.2025) <- cell.sum.SSP1.2025$cell
names(SSP1.2050) <- cell.sum.SSP1.2050$cell
names(SSP3.2050) <- cell.sum.SSP3.2050$cell
names(SSP4.2050) <- cell.sum.SSP4.2050$cell
names(SSP5.2050) <- cell.sum.SSP5.2050$cell

# generar indice:
idx = unique(c(names(SSP1.2025)
               , names(SSP1.2050)
               , names(SSP3.2050)
               , names(SSP4.2050)
               , names(SSP5.2050)
               , names(SSP1.2100)
               , names(SSP3.2100)
               , names(SSP4.2100)
               , names(SSP5.2100)
               ))


# INDEX SHARED SOCIECONOMIC PATHWAYS (SSP) FOR MAP ALGEBRA

SSP1_2025 = SSP1.2025[idx]
names(SSP1_2025) = idx

SSP1_2050 = SSP1.2050[idx]
names(SSP1_2050) = idx

SSP4_2050 = SSP4.2050[idx]
names(SSP4_2050) = idx

SSP5_2050 = SSP5.2050[idx]
names(SSP5_2050) = idx


# HAZARD CHANGE UNDER CONSTRASTING SCENARIOS OF LAND USE
# SUSTAINABLE
SSP1_2050_SSP1_2025 = SSP1.2050[idx]-SSP1.2025[idx]
names(SSP1_2050_SSP1_2025) = idx

# FOSSIL FUELLED SHARED SOCIOECONOMIC PATHWAY
SSP5_2050_SSP1_2025 = SSP5.2050[idx]-SSP1.2025[idx]
names(SSP5_2050_SSP1_2025) = idx

# HAZARD CHANGE IN A SSP4 OF DEEPENING INEQUALITY
# SPP5.2050 - SSP1.2050
SSP4_2050_SSP1_2025 = SSP4.2050[idx]-SSP1.2025[idx]
names(SSP4_2050_SSP1_2025) = idx

# AVOIDED COSTS
# SPP5.2050-SSP1.2050
SSP5_2050_SSP1_2050 = SSP5.2050[idx]-SSP1.2050[idx]
names(SSP5_2050_SSP1_2050) = idx

# SPP4.2050 - SSP1.2050
SSP4_2050_SSP1_2050 = SSP4.2050[idx]-SSP1.2050[idx]
names(SSP4_2050_SSP1_2050) = idx


# REPORT RESULTS IN RASTER
# CREATE RASTERS
# ref: https://stackoverflow.com/questions/19627344/how-to-create-a-raster-from-a-data-frame-in-r

library(raster)
# create spatial points data frame
# CHOOSE THE DATA TO BE PUBLISHED IN A RASTER. NOTE THAT THIS FUNCTION CAN BE PASSED TO A LIST. 
# X = data.frame(SSP4_2100_SSP1_2100) # COSTS
# 
# X = data.frame(SSP3.2050) # HAZARD
# 
# X = data.frame(SSP4_2050_SSP1_2025) # COSTS
# 
# 
# X = data.frame(SSP5_2050)
# 
# X = data.frame(SSP1_2050_SSP1_2025) # COSTS
# X = data.frame(SSP5_2050_SSP1_2025) # COSTS
# X = data.frame(SSP5_2050_SSP1_2050) # COSTS
# 

X = data.frame(richness)

coord = do.call(rbind, strsplit(rownames(X), " "))
lat = as.numeric(coord[,1])
lon = as.numeric(coord[,2])

spg <- data.frame(X, lat, lon) 
coordinates(spg) <- ~ lon + lat

# coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE

# coerce to raster
rasterDF <- raster(spg)
plot(rasterDF)
path_r= "OATH TO RESULTS"

# si los datos son integer
#writeRaster(rasterDF, file=paste(path_r, names(X), ".tif"), format = "GTiff", overwrite=T, datatype="INT1U") # usar este type para integers

# si los datos son costos: -50 a 50
writeRaster(rasterDF, file=paste(path_r, names(X), ".tif"), format = "GTiff", overwrite=T, datatype="FLT4S") # , datatype="INT1U") # usar este type para integers

