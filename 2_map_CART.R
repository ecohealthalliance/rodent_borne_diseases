############################################
# Please make reference of this work as: Gabriel E. García-Peña, André V. Rubio, Hugo Mendoza, Miguel Fernández, Matthew T. Milholland, A. Alonso Aguirre, Gerardo Suzán, Carlos Zambrana-Torrelio. 2021. Land-use change and rodent-borne diseases: Hazards on the shared socioeconomic pathways. Phil. Trans. R. Soc. B.
# Code developer: Gabriel E. García-Peña
#
#
# OBJECTIVES OF THIS R CODE:
# 1) Read historical data copiled into de file: data_gbif_luc.rds # SEE FILE: 1_clean_compile_data.R
# 2) Perform (RPART aka CART) Recursive partitioning / Classification tree analysis (CART), and predict presence in future land-use scenarios.
# 3) Publish resuts for each species as vectors in a map (geojson)
#
#
# Notes: All raster data was obtained from https://luh.umd.edu/data.shtml
# Historical land-use dataset (LUH2 v2h) covers the period 850-2015
############################################
# REFERENCE: https://gis.stackexchange.com/questions/276494/select-pixel-values-in-a-raster-using-r
################################################################

# Load libraries needed
require(raster)
require(rgdal)
require(sp)
require(rpart)

# Read dataset data_gbif_luc.rds with information on sampling (GBIF) and land use. The file was generated with the script 1_clean_compile_data.R
subsam<-readRDS("data_gbif_luc.rds")


# Function (f1) to extract raster values of a variable (varX) in a given year (band), within the IUCN distribution area of each rodent species.
f1<-function(X, x, path, val, band, varX) {
  rasX<-raster(path, band=band, varname = varX)
  p = split(X, X$binomial)
  ix<-p[[x]]
  
  # CROP RASTER
  r2 <- mask(rasX, ix)
  ix <- r2 %in% 1
  r3 <- mask(r2, ix, maskvalue = val, inverse=T)
  r3
}


# CART. Clasification tree analysis or Recursive partitioning (rpart)

# Read species distribution from the multipolygon dataset MAMMALS_TERRESTRIALS_ONLY
X<-readOGR(dsn = "PATH TO THE RASTER FILE", layer = "MAMMALS_TERRESTRIAL_ONLY")

# Select native distributions of Rodentia
rodents<-X[X$order_=="RODENTIA",]
rodents_native<-rodents[rodents$origin==1,]
# Reference: https://nc.iucnredlist.org/redlist/content/attachment_files/Mapping_Standards_Version_1.18_2019.pdf#%5B%7B%22num%22%3A55%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22XYZ%22%7D%2C87%2C436%2C0%5D


# Generate GENERAR LA GRADILLA DE CELDAS id
require(raster)

# Paths to raster values where the prediction will take place.
# For example: RCP8.5 SSP5 (from REMIND-MAGPIE) http://gsweb1vh2.umd.edu/LUH2/LUH2_v2f/MAGPIE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc
path<-"MAGPIE-SSP5-RCP85-2015-2100.nc" # del 2015 al 2100. WORST SCENARIO RCP8.5


# function  f.predict that predicts species presence on the raster, within the distribution area of each species.
f.predict<-function(x) {
    
  all.vars<-c("primf", "primn", "secdf", "secdn"
              , "range", "pastr", "c3per", "c3ann"
              , "c3nfx", "c4per", "c4ann", "urban") # List of variables to include in the analysis
  
  r1<-f1(X = X, x = x , path = path, varX = "primf", val = 0:1, band = 85) # choose a band (year) to predict. Note that the band depends on the raster being used.
  r2<-f1(X = X, x = x , path = path, varX = "primn", val = 0:1, band = 85) 
  r3<-f1(X = X, x = x , path = path, varX = "secdf", val = 0:1, band = 85)
  r4<-f1(X = X, x = x , path = path, varX = "secdn", val = 0:1, band = 85)
  r5<-f1(X = X, x = x , path = path, varX = "range", val = 0:1, band = 85)
  r6<-f1(X = X, x = x , path = path, varX = "pastr", val = 0:1, band = 85)
  r7<-f1(X = X, x = x , path = path, varX = "c3per", val = 0:1, band = 85)
  
  r8<-f1(X = X, x = x , path = path, varX = "c3ann", val = 0:1, band = 85)
  r9<-f1(X = X, x = x , path = path, varX = "c3nfx", val = 0:1, band = 85)
  r10<-f1(X = X, x = x , path = path, varX = "c4per", val = 0:1, band = 85)
  r11<-f1(X = X, x = x , path = path, varX = "c4ann", val = 0:1, band = 85)
  r12<-f1(X = X, x = x , path = path, varX = "urban", val = 0:1, band = 85)
  
  id1 <- rasterToPolygons(r1)
  id2 <- rasterToPolygons(r2)
  id3 <- rasterToPolygons(r3)
  id4 <- rasterToPolygons(r4)
  id5 <- rasterToPolygons(r5)
  id6 <- rasterToPolygons(r6)
  id7 <- rasterToPolygons(r7)
  id8 <- rasterToPolygons(r8)
  id9 <- rasterToPolygons(r9)
  id10 <- rasterToPolygons(r10)
  id11 <- rasterToPolygons(r11)
  id12 <- rasterToPolygons(r12)
  
  # MERGE DATA 
  # GET NAMES OF POLYGONS
  
  f.poli<-function(v) {
    P@polygons[[v]]@ID
  }
  
  P = id1
  v = 1:length(P)
  id1.x<-sapply(v, f.poli)
  
  P = id2
  v = 1:length(P)
  id2.x<-sapply(v, f.poli)
  
  P = id3
  v = 1:length(P)
  id3.x<-sapply(v, f.poli)
  
  P = id4
  v = 1:length(P)
  id4.x<-sapply(v, f.poli)
  
  P = id5
  v = 1:length(P)
  id5.x<-sapply(v, f.poli)
  
  P = id6
  v = 1:length(P)
  id6.x<-sapply(v, f.poli)
  
  P = id7
  v = 1:length(P)
  id7.x<-sapply(v, f.poli)
  
  P = id8
  v = 1:length(P)
  id8.x<-sapply(v, f.poli)
  
  P = id9
  v = 1:length(P)
  id9.x<-sapply(v, f.poli)
  
  P = id10
  v = 1:length(P)
  id10.x<-sapply(v, f.poli)
  
  P = id11
  v = 1:length(P)
  id11.x<-sapply(v, f.poli)
  
  P = id12
  v = 1:length(P)
  id12.x<-sapply(v, f.poli)
  
  id1@data <- data.frame(id1@data
                         , id2@data[match(id1.x, id2.x),]
                         , id3@data[match(id1.x, id3.x),]
                         , id4@data[match(id1.x, id4.x),]
                         , id5@data[match(id1.x, id5.x),]
                         , id6@data[match(id1.x, id6.x),]
                         , id7@data[match(id1.x, id7.x),]
                         , id8@data[match(id1.x, id8.x),]
                         , id9@data[match(id1.x, id9.x),]
                         , id10@data[match(id1.x, id10.x),]
                         , id11@data[match(id1.x, id11.x),]
                         , id12@data[match(id1.x, id12.x),])
  
  names(id1@data)<-all.vars
  
# Collect presence data and a sample of the same size, from other rodent species
  presence.x<-subsam[grep(x, subsam$species), all.vars]   
  other.rodent<-subsam[-grep(x, subsam$species), all.vars]
  
# Random sample
  s<-other.rodent[sample(1:nrow(other.rodent), nrow(presence.x)),] 
  presence.x$occurence = rep(1, nrow(presence.x))
  s$occurence = rep(0, nrow(s))
  
# Recursive partitioning / Classification tree analysis (CART)  
  dat<-rbind(presence.x, s)
  names(dat)<-c(all.vars, "occurence")
  
# Model
  fx<-occurence ~ primn + primf + secdn + secdf + pastr + range + c3ann + c3per + c4ann + c4per + c3nfx  
  tr.raw<-rpart(fx, data=dat, method="class")

###################################################################!!!!!!!!!  
#https://www.statmethods.net/advstats/cart.html
# Prune the tree to avoid over parametrization, based on the cross validation error
  cp=tr.raw$cptable[which.min(tr.raw$cptable[,"xerror"]),"CP"]
  ptree<- prune(tr.raw, cp = tr.raw$cptable[which.min(tr.raw$cptable[,"xerror"]),"CP"])
  
  V<-ptree # Pruned tree for further analyses
  
  
  names(id1)<-all.vars
  guess<-id1@data
  
  y<-predict(V, guess)
  
  # ADD COLUMN OF ABSENCE PRESENCE TO THE POLYGON id1. 
  # CREATE NEW DATAFRAME P
  
  P<-id1
  P@data<-data.frame(P@data, y)
  
  # SAVE P AS geoJson
  
  require(leaflet)
  require(leafletR)
  
  
  toGeoJSON(P, x)
  
}


# Run the function f.predict()
# f.predict requires:
#
# subsam = The data with GBIF locations and historical data on land use
# path = a full path to the raster with the land use scenarios. Note: land use scenarios are 
# X = a polygon of the native distribution of rodent species.
# nms = a list of species to model (named as in X, and the data on land use. For example, 
#
# Example
path<-"PATH TO THE RASTER" # del 2015 al 2100. 
X = rodents_native
nms= c("Aconaemys porteri", "Microsciurus alfari", "Spermophilus dauricus", "Oligoryzomys fornesi")


f.predict(nms[1]) # The prediction is generated and published in a map (geojson), saved in the working directory with the name of the species
