coo = read.csv('C:/Users/ac22qawo/chelsa/all_coord.csv')

cia = subset(coo, project == 'space_for_time' )
select(cia,lat,lon) %>% dim
select(cia,lat,lon) %>% unique %>% dim

  
clm_df <- read.csv('C:/CODE/space_for_time/data/chelsa_mon_7913.csv')
prj_df <- read.csv('C:/CODE/space_for_time/data/grid_usa_projections.csv')
  


sd_df <- clm_df %>% 
  subset( variable == 'tmean') %>% 
  mutate( value = (value/10)-273 ) %>%
  # subset( lat > 38 & lat < 39 ) %>% 
  # subset( lon > -123 & lon < -122 ) %>% 
  select( year, month, lat, lon, variable, value ) %>% 
  group_by(year, lat, lon ) %>% 
  summarise( mean_t = mean(value) ) %>% 
  ungroup %>% 
  group_by(lat, lon ) %>% 
  summarise( sd_t = sd(mean_t) )

cia_df <- full_join(sd_df, prj_df)

plot(value~sd_t,data=cia_df)

cia_df %>% head

cia_df <- sft_df %>% 
  subset( variable == 'tmean') %>% 
  subset( lat > 36 & lat < 40 ) %>% 
  subset( lon > -100 & lon < -95 ) %>% 
  mutate( value = (value/10)-273 ) %>% 
  select( year, month, lat, lon, variable, value )

# cia_df %>% 
#   group_by(year, lat, lon ) %>% 
#   summarise( mean_t = mean(value) ) %>% 
#   ungroup %>% 
#   mutate( cia = paste0(lat,lon) ) %>% 
#   mutate( cia = as.factor(cia) ) %>% 
#   mutate( group = as.numeric(cia) ) %>% 
#   subset( group == 5 ) %>% 
#   ggplot() +
#   geom_histogram( aes(mean_t) )

all.equal(sft_df, clm_df)

# sft_df %>% 
clm_df %>% 
  subset( variable == 'tmean') %>% 
  mutate( value = (value/10)-273 ) %>%
  # subset( lat > 38 & lat < 39 ) %>% 
  # subset( lon > -123 & lon < -122 ) %>% 
  select( year, month, lat, lon, variable, value ) %>% 
  group_by(year, lat, lon ) %>% 
  summarise( mean_t = mean(value) ) %>% 
  ungroup %>% 
  group_by(lat, lon ) %>% 
  summarise( sd_t = sd(mean_t) ) %>% 
  ungroup %>% 
  # .$sd_t %>% hist
  ggplot() +
  geom_point( aes( lon, lat,
                   color = sd_t
  ),
  size = 3,
  alpha = 0.5 ) +
  scale_color_viridis_c()
  

select( lat, lon, variable, year, month, value) %>% 
  unique %>% 
  dim