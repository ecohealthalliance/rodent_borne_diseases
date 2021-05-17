# LAS  PREDICCIONES EN EL PRESENTE (2015) SE HICIERON CON ESTOS ARCHIVOS: LUH2 v2h Release (10/14/16): The updated release of the historical land-use forcing dataset (LUH2 v2h) covers the period 850-2015 and corrects all known issues and notices identified with the previous version (LUH2 v1.0h). This dataset replaces the previously released dataset (LUH2 v1.0h). This product is the result of a series of prototypes released previously, uses the established data format, and will connect smoothly to gridded products for the future.
# LAS PREDUCCIONES EN L FUTURO SE HACE CON ESTOS OTROS ARCHIVOS: LUH2 v2f Release (12/21/17)
# LA PAGINA DE REFERENCIA ES: https://luh.umd.edu/data.shtml
#

# MATERIALES
# DATOS DE GBIF CON LAS VARIABLES EN CADA AÑO DE MUESTREO EN GBIF. 

# Make predictions for each species inside its native distribution
# FILTER RASTER BASED ON OBSERVED VALUES
# METHOD REF: https://gis.stackexchange.com/questions/276494/select-pixel-values-in-a-raster-using-r
#

require(raster)
require(rgdal)
require(sp)
require(rpart)

# Read data
subsam<-readRDS("~/LUCIDA/data/data_gbif_luc.rds")

# FUNCION f1 PARA EXTRAER VALORES (val) DEL RASTER DE USO DE SUELO (path, band) DENTRO DEL AREA DE DISTRIBUCIÓN DE LA IUCN (X)
f1<-function(X, x, path, val, band, varX) {
  rasX<-raster(path, band=band, varname = varX)
  # LOAD FILE WITH IUCN DISTRIBUTIONS IN /data/iucn/
  p = split(X, X$binomial)
  ix<-p[[x]]
  
  # CROP RASTER
  r2 <- mask(rasX, ix)
  # plot(r2)
  ix <- r2 %in% 1
  r3 <- mask(r2, ix, maskvalue = val, inverse=T)
  r3
}



# CLASSIFICATION TREE
########## ausencias
# LEER MULTIPOLIGONOS CON LAS DISTRIBUCIONES DE LA IUCN
#X<-readOGR(dsn = "~/LUCIDA/data/iucn/", layer = "roedors_mammals_01") # ADD MULTI POLYGONS
X<-readOGR(dsn = "~/LUCIDA/data/iucn/", layer = "MAMMALS_TERRESTRIAL_ONLY") # ADD MULTI POLYGONS
rodents<-X[X$order_=="RODENTIA",]
rodents_native<-rodents[rodents$origin==1,]

##############
# METADATA reference: 
# https://nc.iucnredlist.org/redlist/content/attachment_files/Mapping_Standards_Version_1.18_2019.pdf#%5B%7B%22num%22%3A55%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22XYZ%22%7D%2C87%2C436%2C0%5D
##############

# ############## GENERAR LA GRADILLA DE CELDAS id
require(raster)

# PATHS A LAS COBERTURAS
#path<-"~/LUCIDA/data/LUH2 v2h Release_10_14_16/states.nc" # del 850 al 2015
#path<-"~/LUCIDA/data/LUH2 v2f_12_21_17/LUH_2015_2100.nc" # del 2015 al 2100. ESCENARIO VERDE RCP2.6

# RCP8.5 SSP5 (from REMIND-MAGPIE) http://gsweb1vh2.umd.edu/LUH2/LUH2_v2f/MAGPIE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc
path<-"~/LUCIDA/data/LUH2_SSP5_RCP85/MAGPIE-SSP5-RCP85-2015-2100.nc" # del 2015 al 2100. WORST SCENARIO RCP8.5



# LA FUNCION QUE HACE TODAS LAS PREDICCIONES...
f.predict<-function(x) {
  
  # OBTENER SERIE DE TIEMPO DE LOS VALORES DE USO DE SUELO DENTRO DEL ARES DE DISTRIBUCION DE LA ESPECIE FOCAL (x)
  
  # NOTE: THE ORDER OF r* MUST FOLLOW all.vars
  
  all.vars<-c("primf", "primn", "secdf", "secdn"
              , "range", "pastr", "c3per", "c3ann"
              , "c3nfx", "c4per", "c4ann", "urban") 
  
  r1<-f1(X = X, x = x , path = path, varX = "primf", val = 0:1, band = 85) # choose a band (850-2015) to predict
  r2<-f1(X = X, x = x , path = path, varX = "primn", val = 0:1, band = 85) # YEAR 2050, 2025, 2075
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
  
  # 
  presence.x<-subsam[grep(x, subsam$species), all.vars]       # agregar all.vars en subsam!
  other.rodent<-subsam[-grep(x, subsam$species), all.vars]
  
  # Select random sample
  s<-other.rodent[sample(1:nrow(other.rodent), nrow(presence.x)),] 
  
  presence.x$occurence = rep(1, nrow(presence.x))
  s$occurence = rep(0, nrow(s))
  
  # GENERATE THE MODEL
  
  dat<-rbind(presence.x, s)
  
  names(dat)<-c(all.vars, "occurence")
  
  # fx<-occurence ~ crop +secd + past + gothr  
  fx<-occurence ~ primn + primf + secdn + secdf + pastr + range + c3ann + c3per + c4ann + c4per + c3nfx  
  #V<-rpart(fx, data=dat, method="class")
  
  tr.raw<-rpart(fx, data=dat, method="class")
  ###################################################################!!!!!!!!!  
  #https://www.statmethods.net/advstats/cart.html
  # USE THRESHOLD TO AVOID OVEPARAMTRIZATION # COMMENTS FROM REVIEWERS
  cp=tr.raw$cptable[which.min(tr.raw$cptable[,"xerror"]),"CP"]
  ptree<- prune(tr.raw, cp = tr.raw$cptable[which.min(tr.raw$cptable[,"xerror"]),"CP"])
  
  #V<-rpart(fx, data=dat, method="class") # cambia al incluir el threshold
  V<-ptree
  
  
  
  
  # names(id1)<-c("gothr","secd", "past", "crop")
  names(id1)<-all.vars
  guess<-id1@data
  
  y<-predict(V, guess)
  # colnames(y)<-paste(x, c(".0", ".1"), sep="") # !cambiar nombres para hacer un archivo multiespecies
  
  # ADD COLUMN OF ABSENCE PRESENCE TO THE POLYGON id1. 
  # CREATE NEW DATAFRAME P
  
  P<-id1
  P@data<-data.frame(P@data, y)
  
  # SAVE P AS geoJson
  
  require(leaflet)
  require(leafletR)
  
  
  toGeoJSON(P, x)
  
}
