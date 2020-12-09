library(leaflet)
library(dplyr)
library(measurements)
library(geosphere)
library(ggplot2)
library(ggmap)
library(SPEI)
library(cowplot)

# climate
clim_df   <- data.table::fread('data/chelsa_mon_7913.csv')
 
# coordinates on the sea surface
sea_df <- clim_df %>% 
            subset( variable == 'prec' ) %>% 
            subset( value == 65535 ) %>% 
            dplyr::select( lat, lon ) %>% 
            unique %>% 
            as.data.frame


# clean bad values
prec_df   <- subset(clim_df, variable == 'prec' ) %>% 
               mutate( value = replace(value, value == 65535, NA),
                       value = replace(value, value > 1000, NA) )


# mean temperature
mean_df   <- subset(clim_df, variable == 'tmean' ) %>% 
               mutate( value = (value/10) - 273 ) %>% 
               dplyr::select( lat, lon, variable,
                               year, month, value )
               right_join( sea_df ) #%>% 
               left_join( prec_df )

# calculate SPEI 
all_df   <- mean_df %>% 
              rename( tmean = value ) 
              
coord_i  <- sea_df[1,]


# create placeholder for output    
spei_df <- mean_df %>% 
  subset( lat ==  coord_i$lat & lon == coord_i$lon ) %>% 
  subset( !is.na(value) ) %>% 
  mutate( PET = thornthwaite(value, unique(lat) ) ) %>% 
  mutate( BAL = prec - PET )

# last year/month combination
yr_mon  <- spei_df %>% 
  dplyr::select( year, month ) %>% 
  unique %>% 
  arrange( year, month ) %>% 
  .[nrow(spei_df),] %>% 
  as.numeric

# transform in 
bal_ts <- ts( spei_df$BAL, 
              end = yr_mon, 
              frequency = 12 )

spei_v <- spei( bal_ts, 12)$fitted %>% as.numeric

# final spei file
spei_df %>% 
  mutate( spei = spei_v ) %>% 
  dplyr::select( SpeciesAuthor, MatrixPopulation, 
                 Longitude, Latitude, year, month, spei )

  

# get yearly standard deviations 
mean_sd    <- mean_df %>% 
                # calculate yearly means
                group_by( year, lat, lon ) %>% 
                summarise( t_mean = mean( value ) ) %>% 
                ungroup %>% 
                # calculate SD across years at each coord
                group_by( lat, lon ) %>% 
                summarise( t_sd = sd( t_mean ) ) %>% 
                ungroup

prec_sd <- prec_df %>% 
            # calculate yearly means
            group_by( year, lat, lon ) %>% 
            summarise( p_sum = sum( value ) ) %>% 
            ungroup %>% 
            # calculate SD across years at each coord
            group_by( lat, lon ) %>% 
            summarise( p_sd = sd( p_sum ),
                       p_m  = mean( p_sum ) ) %>% 
            ungroup %>% 
            mutate( p_cv = p_sd / p_m )


# No precipitation data for marine sites 
# get terrestrial coords!
sea_crd <- subset(prec_sd, is.na(p_sd) ) %>% 
                select(lat,lon) %>% 
                unique


# # e.g. Lupine's site
# subset(mean_sd, 
#        lat > 38 & lat < 39 ) %>% 
#   subset( lon > -123 & lon < -122 )
# 
# # map of the whole grid (INCLUDES OCEAN)
# ggplot( anti_join(prec_sd, sea_crd) ) +
#   geom_point( aes( lon, lat,
#                    color = p_cv
#                    # color = p_sd
#                    # color = p_m
#                    ),
#               size = 3,
#               alpha = 0.5 ) +
#   scale_color_viridis_c()

# # 
# ggplot( anti_join(mean_sd, sea_crd) ) +
#   geom_point( aes( lon, lat,
#                    color = p_sd,
#                    color = p_sd),
#   size = 3,
#   alpha = 0.5 ) +
#   scale_color_viridis_c()


# get yearly values
yr_df     <- mean_df %>% 
               dplyr::select(lat,lon,year,value) %>% 
               group_by( lat, lon, year ) %>% 
               summarise( tmean = mean(value, na.rm=T) ) %>% 
               ungroup %>% 
               anti_join( sea_df )

yr_pre_df <- prec_df %>% 
               dplyr::select(lat,lon,year,value) %>% 
               group_by( lat, lon, year ) %>% 
               summarise( pmean = mean(value, na.rm=T) ) %>% 
               ungroup %>% 
               anti_join( sea_df )


coord_df <- yr_df %>% 
              dplyr::select( lat, lon ) %>% 
              unique %>% 
              anti_join( sea_df )

# scale values
scale_temp <- function(ii){
  
  subset(yr_df,
         lat == coord_df$lat[ii] & 
         lon == coord_df$lon[ii] ) %>% 
    mutate( scaled = scale(tmean) ) %>% 
    rename( value  = tmean )
          
}

scale_prec <- function(ii){
  
  subset(yr_pre_df,
         lat == coord_df$lat[ii] & 
         lon == coord_df$lon[ii] ) %>% 
    mutate( scaled = scale(pmean) ) %>% 
    rename( value  = pmean )
  
}


yr_t_sc <- lapply(1:nrow(coord_df), scale_temp) %>% bind_rows
yr_p_sc <- lapply(1:nrow(coord_df), scale_prec) %>% bind_rows

# calculate distances

# good coordinates:
# 967 (lower pacific coast) 
# 3000 (seattle)
# 1585 (central plains)
# 1950 (Boston)
# 39 (Miami)


id <- 39


# get interactions
int_df <- expand.grid( x = id,
                       y = 1:nrow(coord_df) )

# scale values
calc_corr <- function(ii, yr_df){
  
  if( int_df$x[ii] != int_df$y[ii] ){
    x_df <- subset(yr_df,
                   lat == coord_df$lat[int_df$x[ii]] & 
                   lon == coord_df$lon[int_df$x[ii]] ) %>% 
              dplyr::select( -lat, -lon ) %>% 
              rename( scaled1 = scaled,
                      value1  = value )
    y_df <- subset(yr_df,
                   lat  == coord_df$lat[int_df$y[ii]] & 
                   lon == coord_df$lon[int_df$y[ii]] ) %>% 
              dplyr::select( -lat, -lon ) %>% 
              rename( scaled2 = scaled,
                      value2  = value )
    
    cor_df <- inner_join( x_df, y_df )
    
    return( cor(cor_df$scaled1, cor_df$scaled2) )
    
  }else{
    return(NA)
  }
  
}

cor_t <- sapply(1:nrow(int_df), calc_corr, yr_t_sc)
cor_p <- sapply(1:nrow(int_df), calc_corr, yr_p_sc)


# create distance matrix
mat <- distm( dplyr::select(coord_df, lon, lat)[id,],
              dplyr::select(coord_df, lon, lat),
              fun = distVincentyEllipsoid )

# correlation of versus distance!
cor_dist_df <- data.frame( cor    = cor_t ,
                           meters = as.numeric(mat) ) %>% 
                  mutate( Kilometers = meters / 1000 ) %>% 
                  bind_cols( coord_df )
  
# plot correlation by distance!
mod <- lm(cor ~ Kilometers, data=cor_dist_df)


# plot(cor ~ Kilometers, data=cor_dist_df)
# abline(mod)

ggplot(cor_dist_df) +
  geom_point( aes( x = Kilometers,
                   y = cor ) ) +
  theme_minimal() +
  geom_hline( yintercept = 0, lty = 2 ) +
  labs( y = 'Correlation in annual anomaly') +
  ggsave( paste0('results/tempcorr_vs_distance_',id,'.tiff'),
          width = 6.3, height = 6.3, compression = 'lzw' )

# set up the map
coords.data <- read.csv("C:/CODE/download_chelsa_data/data/grid_usa.csv")
map_bounds <- c(left = -125,  bottom = 25,
                right = -65,  top = 50)
coords.map <- get_stamenmap(map_bounds, zoom = 7, maptype = "toner-lite")
coords.map <- ggmap(coords.map, extent="device", legend="none")

# map of delta_T ("multiple")
map_df <- cor_dist_df %>% 
            # remove focal coordinate
            subset( !(Kilometers == 0) ) %>% 
            anti_join( sea_df )

coords.map + 
  geom_tile(data = map_df, 
            aes(x = lon, 
                y = lat, 
                fill = cor),
            alpha = 0.8) + 
  geom_point( data = subset(cor_dist_df, Kilometers == 0),
              aes( x = lon,
                   y = lat),
              color ='red' ) +
  scale_fill_viridis_c() +
  labs( fill = 'Correlation' ) +
  ggsave( paste0('results/cor_t_',id,'central.tiff'),
          width=6.3,height=6.3,compression='lzw')


# # map of sd_T ("sd in annual temperature")
# coords.map + 
#   geom_tile(data = subset(comp_df, !is.nan(multiple) ), 
#             aes(x = lon, 
#                 y = lat, 
#                 fill = sd_t),
#             alpha = 0.8) + 
#   scale_fill_viridis_c() +
#   labs( fill = "sd(T)" ) +
#   ggsave( 'results/climate/sd_annual_t_usa.tiff',
#           width=6.3,height=6.3,compression='lzw')


# -----------------------------------------------------------

ids <- c(967, 3000, 1585, 1950, 39)
# 967 (lower pacific coast) 
# 3000 (seattle)
# 1585 (central plains)
# 1950 (Boston)
# 39 (Miami)


# loop through the 5 sites.
for( ii in 1:5){
  
  id <- ids[ii]
  
  
  # get interactions
  int_df <- expand.grid( x = id,
                         y = 1:nrow(coord_df) )
  
  # scale values
  calc_corr <- function(ii, yr_df){
    
    if( int_df$x[ii] != int_df$y[ii] ){
      x_df <- subset(yr_df,
                     lat == coord_df$lat[int_df$x[ii]] & 
                       lon == coord_df$lon[int_df$x[ii]] ) %>% 
        dplyr::select( -lat, -lon ) %>% 
        rename( scaled1 = scaled,
                value1  = value )
      y_df <- subset(yr_df,
                     lat  == coord_df$lat[int_df$y[ii]] & 
                       lon == coord_df$lon[int_df$y[ii]] ) %>% 
        dplyr::select( -lat, -lon ) %>% 
        rename( scaled2 = scaled,
                value2  = value )
      
      cor_df <- inner_join( x_df, y_df )
      
      return( cor(cor_df$scaled1, cor_df$scaled2) )
      
    }else{
      return(NA)
    }
    
  }
  
  cor_t <- sapply(1:nrow(int_df), calc_corr, yr_t_sc)
  cor_p <- sapply(1:nrow(int_df), calc_corr, yr_p_sc)
  
  
  # create distance matrix
  mat <- distm( dplyr::select(coord_df, lon, lat)[id,],
                dplyr::select(coord_df, lon, lat),
                fun = distVincentyEllipsoid )
  
  # correlation of versus distance!
  cor_dist_df <- data.frame( cor    = cor_t ,
                             meters = as.numeric(mat) ) %>% 
    mutate( Kilometers = meters / 1000 ) %>% 
    bind_cols( coord_df )
  
  
  # correlation of versus distance - temperature and precipitation
  cor_t_dist_df <- data.frame( cor    = cor_t ,
                             meters = as.numeric(mat) ) %>% 
                    mutate( Kilometers = meters / 1000 ) %>% 
                    bind_cols( coord_df )
  
  cor_p_dist_df <- data.frame( cor    = cor_p,
                               meters = as.numeric(mat) ) %>% 
                    mutate( Kilometers = meters / 1000 ) %>% 
                    bind_cols( coord_df )
  
  
  # correlation of versus distance - temperature and precipitation
  map_t_df <- cor_t_dist_df %>% 
                # remove focal coordinate
                subset( !(Kilometers == 0) ) %>% 
                anti_join( sea_df )
  map_p_df <- cor_p_dist_df %>% 
                # remove focal coordinate
                subset( !(Kilometers == 0) ) %>% 
                anti_join( sea_df )
  
  
  # plot(cor ~ Kilometers, data=cor_dist_df)
  # abline(mod)
  
  tt <- ggplot(cor_t_dist_df) +
          geom_point( aes( x = Kilometers,
                           y = cor ),
                      alpha = 0.3 ) +
          theme_minimal() +
          geom_hline( yintercept = 0, lty = 2 ) +
          labs( y = 'Temp. anomalies corr.' )
  
  
  pp <- ggplot(cor_p_dist_df) +
          geom_point( aes( x = Kilometers,
                           y = cor ),
                      alpha = 0.3 ) +
          theme_minimal() +
          geom_hline( yintercept = 0, lty = 2 ) +
          labs( y = 'Prec. anomalies corr.')
  
  # map temperature
  map_t <- coords.map + 
            geom_tile(data = map_t_df, 
                      aes(x = lon, 
                          y = lat, 
                          fill = cor),
                      alpha = 0.8) + 
            geom_point( data = subset(cor_dist_df, Kilometers == 0),
                        aes( x = lon,
                             y = lat),
                        color ='red' ) +
            scale_fill_viridis_c() +
            labs( fill = 'Temp. corr.' )
  
  # map temperature
  map_p <- coords.map + 
    geom_tile(data = map_p_df, 
              aes(x = lon, 
                  y = lat, 
                  fill = cor),
              alpha = 0.8) + 
    geom_point( data = subset(cor_dist_df, Kilometers == 0),
                aes( x = lon,
                     y = lat),
                color ='red' ) +
    scale_fill_viridis_c() +
    labs( fill = 'Prec. corr.' )
  
  h_4 <- plot_grid( tt,  map_t, 
                    pp,  map_p,
                    labels = 'AUTO',
                    label_size = 18,
                    align='h', 
                    nrow = 2, ncol = 2 )
  
  ggsave( paste0('results/temp_prec_',id,'.tiff'),
         h_4,
         dpi=600, width=6.3, height=6.3 )

}
