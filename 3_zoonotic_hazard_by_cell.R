############################################
# Please make reference of this work as: Gabriel E. García-Peña, André V. Rubio, Hugo Mendoza, Miguel Fernández, Matthew T. Milholland, A. Alonso Aguirre, Gerardo Suzán, Carlos Zambrana-Torrelio. 2021. Land-use change and rodent-borne diseases: Hazards on the shared socioeconomic pathways. Phil. Trans. R. Soc. B.
# Code developer: Gabriel E. García-Peña
#
#
# OBJECTIVES OF THIS R CODE:
# 1) Compile all the predictions from geojson files generated with 2_map_CART.R
# 2) Filter out species with less than 10 observations
# 3) Calculate zoonotic hazard, and richness for each geographic cell
# 
################################################################
# Requirements:
# nms = a list of species for which a geojson exists
# path = access to geojson files
# ref = data frame with values of equilibrium prevalence.


# Required libraries
require(geojsonio)
require(rgdal)


nms1 = c("NAMES OF SPECIES MODELED (geojson")


# Function that reads the files in the path, based on the list nms1
  f3<-function(x){
    d <- readOGR(dsn = paste("PATH TO DIRECTORY WITH geojson", nms1[x], ".geojson", sep="")) # SSP3 2100
    
    cord<-coordinates(d)
    lon<-cord[,1]
    lat<-cord[,2]
    d@data$lat = round(lat*1000)/1000
    d@data$lon = round(lon*1000)/1000
    d@data$cell= paste(d@data$lat, d@data$lon)
    d@data$specie = nms1[x]
    Q = d@data
    Q
  }  

  X<-lapply(1:length(nms1), f3) # Generates a list with cells and presence prediction for each species
  saveRDS(X, "~/LUCIDA/output/cell_analysis/SSP3_2100.rds") # Save the data colected as back up


# Asure that sample size > 10 observations in all species analysed
  nms<-gsub("_", " ", nms1) # list of species in nms1 with "_" between genera and species names


# Read GBIF data
  subsam<-readRDS("data_gbif_luc.rds")


# Function to check sample size
  f.check<-function(x){
    z<-table(match(subsam$species, nms[x])) # Nota: deberia ser con subsam en vez de sam, pero debe ser igual. subsam es una submuestra con las especies de interés.
    names(z)<-nms[x]
    z
  }

# Vector with samples sizes of each species
  sample.size<-sapply(c(1:length(nms)), f.check)

# Indicate species with sample size > 10
  sample_size<-sample.size[sample.size>10]

# Fetch data of species with sample size > 10
  names(X)<-gsub("_", " ", nms)
  X2<-X[names(sample_size)]

  D<-do.call(rbind, X2) # Working data frame 

# Read data on equilibrium prevalence from Hann et al. 2020, and add compile them into D
  ref<-read.csv("~/LUCIDA/maps/no_headers/equilibrium_prevalence_27feb2021.csv", header=T)
  D$eq_prevalence<-ref$equilibrium_prevalence[match(D$specie, ref$spp_name_corrected)]

# Check taxonomic synonyms
# reservoirs$species<-gsub("Spermophilus_beecheyi", "Otospermophilus_beecheyi", reservoirs$species)
# reservoirs$species<-gsub("Spermophilus_variegatus", "Otospermophilus_variegatus", reservoirs$species)
# reservoirs$species<-gsub("Spermophilus_lateralis", "Callospermophilus_lateralis", reservoirs$species)
# reservoirs$species<-gsub("Oryzomys_megacephalus", "Hylaeamys_megacephalus", reservoirs$species)
# reservoirs$species<-gsub("Oryzomys_nitidus", "Euryoryzomys_nitidus", reservoirs$species)

# Multiply presence probability (X1) times equilibrium prevalence
  D$X1.num<-as.numeric(D$X1)
  D$Y = D$X1.num * D$eq_prevalence
  D_nona<-na.exclude(D)  # Exlude species without equilibrium prevalence

# Calculate zoonotic hazard in ach geographic cell
  cell.sum.reservoirs<-aggregate(D_nona$Y, list(D_nona$cell), sum)
  names(cell.sum.reservoirs)<-c("cell", "cell.sum.reservoirs")

# Estimate richness of taxa considered in each geographic cell
  r<-as.data.frame(table(D$cell))
  names(r)<-c("cell", "reservoir.richness")

# Compile data on richness
  cell.sum.reservoirs$reservoir.richness<-r$reservoir.richness[match(cell.sum.reservoirs$cell, r$cell)]

# Add coordinates

  cell.sum.reservoirs$lat<-D$lat[match(cell.sum.reservoirs$cell, D$cell)]
  cell.sum.reservoirs$lon<-D$lon[match(cell.sum.reservoirs$cell, D$cell)]


# Save results of hazard as csv
  write.csv(cell.sum.reservoirs, "cell_sum_N10_SSP_year.csv") # PROSPERO
