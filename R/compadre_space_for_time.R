# statistics cited in the article
#1. space for time for populations
#2. histogram space distance same study in compadre 
#3. histogram space distance same SPECIES in compadre
setwd('C:/Users/ac22qawo/Dropbox/space_for_time')
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(ggthemes)
library(testthat)
library(geosphere)
options(stringsAsFactors = F)


# Space-for-time COMPADRE replication
load('C:/CODE/plant_review_analyses/Data/COMPADRE_v.5.0.0.RData')

# studies in compadre
studies_v <- compadre$metadata$SpeciesAuthor %>% unique

# space and time replication by study in compadre
space_x_time <- function( spp ){
  
  st_df <- compadre$metadata %>% 
    subset( SpeciesAuthor == spp ) %>% 
    subset( MatrixComposite == 'Individual' )
  
  sp_rep_n <- select(st_df, MatrixPopulation, Lat, Lon) %>% 
                unique %>% 
                select( Lat, Lon ) %>% 
                unique %>% 
                nrow
  
  data.frame( yr_rep = st_df$MatrixStartYear %>% unique %>% length,
              sp_rep = sp_rep_n )
  
}

# space for time replication
rep_comp <- lapply(studies_v, space_x_time) %>% 
              bind_rows %>% 
              subset( !(yr_rep == 0 & sp_rep == 0) )

# space for time replication plot
ggplot( rep_comp ) + 
  geom_point( aes( x = sp_rep,
                   y = yr_rep),
              alpha = 0.2 ) +
  theme_minimal() + 
  labs( x = 'Spatial replication (sites)',
        y = 'Temporal replication (years)' ) +
  ggsave( 'space_x_time_rep.png' )


#2. histogram space distance same study in compadre ----------------

# studies in compadre
studies_v <- compadre$metadata$SpeciesAuthor %>% unique

# spp=studies_v[140]

# space and time replication by study in compadre
space_x_time <- function( spp ){
  
  st_df <- compadre$metadata %>% 
    subset( SpeciesAuthor == spp ) %>% 
    subset( MatrixComposite == 'Individual' ) %>% 
    dplyr::select(Lat,Lon) %>% 
    unique %>% 
    rename( Longitude = Lon,
            Latitude  = Lat ) %>% 
    subset( !(is.na(Latitude) | is.na(Longitude)) )
  
 
  if( nrow(st_df) > 1){
    
    # create distance matrix
    all_dist <- distm(dplyr::select(st_df, Longitude, Latitude),
                      dplyr::select(st_df, Longitude, Latitude),
                      fun = distVincentyEllipsoid) %>% 
      as.numeric %>% 
      unique %>% 
      Filter( function(x) x!=0, .)
    
    # return results
    data.frame( study  = spp,
                km     = all_dist / 1000 ) %>% 
      return()
  }
  else{
    return(NULL)
  }
  
}

# produce distances
dist_studies <- sapply(studies_v, space_x_time) %>% bind_rows

# put out graph
tiff('distance_within_study.tiff',
     height=6.3, width=6.3, unit='in', res=600,
     compression='lzw')

hist( unique(dist_studies$km),
      xlab = 'Distance (Km)',
      main = 'Distance of sites within same STUDY' )
abline( v = mean(unique(dist_studies$km)),
        lwd=2, col='red',lty=2)

dev.off()


#3. histogram space distance same SPECIES in compadre ----------

spp_v <- compadre$metadata$SpeciesAccepted %>% unique 
                
# space and time replication by study in compadre
space_x_time <- function( spp ){
  
  st_df <- compadre$metadata %>% 
    subset( SpeciesAccepted == spp ) %>% 
    subset( MatrixComposite == 'Individual' ) %>% 
    dplyr::select(Lat,Lon) %>% 
    unique %>% 
    rename( Longitude = Lon,
            Latitude  = Lat ) %>% 
    subset( !(is.na(Latitude) | is.na(Longitude)) )
  
  
  if( nrow(st_df) > 1){
    
    # create distance matrix
    all_dist <- distm(dplyr::select(st_df, Longitude, Latitude),
                      dplyr::select(st_df, Longitude, Latitude),
                      fun = distVincentyEllipsoid) %>% 
      as.numeric %>% 
      unique %>% 
      Filter( function(x) x!=0, .)
    
    # return results
    data.frame( study  = spp,
                km     = all_dist / 1000 ) %>% 
      return()
  }
  else{
    return(NULL)
  }
  
}

# produce distances per SPECIES
dist_spp <- sapply(spp_v, space_x_time) %>% bind_rows

tiff('distance_within_species.tiff',
     height=6.3, width=6.3, unit='in', res=600,
     compression='lzw')

hist( unique(dist_spp$km),
      xlab = 'Distance (Km)',
      main = 'Distance of sites within same SPECIES' )
abline( v = mean(unique(dist_studies$km)),
        lwd=2, col='red',lty=2)

dev.off()


# # abstract numbers -----------------------------
# 
# tab1 <- read.csv('Data/table_1_update.csv')
# 
# # number of studies
# tab1$Authors %>% 
#   unique %>% 
#   strsplit(';') %>% 
#   unlist %>% 
#   trimws %>% 
#   length
