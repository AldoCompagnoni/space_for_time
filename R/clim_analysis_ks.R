# Use Adler's data to look at patterns of abundance
library(tidyverse)
library(testthat)
options(stringsAsFactors = F)

# read data 
ks        <- read.csv('data/adler/kansas/allrecords.csv')
ks_q      <- read.csv('data/adler/kansas/quadrat_info.csv')

# prec/temp anomalies
prec_anom <- read.csv('data/adler/kansas/cruts_hays_prec_anom.csv')
temp_anom <- read.csv('data/adler/kansas/cruts_hays_temp_anom.csv')


# Format --------------------------------------------------

# extract just last two digits...
last_two_digits <- function( plotyr ){
  regmatches(plotyr, 
             gregexpr("[[:digit:]]{2}$", 
                      plotyr) ) %>% unlist
}

# formatted kansas dataset
ks_df <- ks %>% 
          mutate( yr      = last_two_digits( plotyear ) ) %>% 
          mutate( nc      = nchar(plotyear)-2 ) %>% 
          mutate( quadrat = substr(plotyear, 1, nc) ) %>% 
          left_join( ks_q ) %>% 
          mutate( yr      = paste0('19',yr) %>% as.numeric )
    

# Analysis ------------------------------------------------

# most abundant species
spp_abund <- ks_df %>% 
  count( species, quadrat, yr) %>% 
  mutate( pres = as.numeric(n > 0) ) %>% 
  group_by( species ) %>% 
  summarise( n    = sum(n),
             pres = sum(pres) ) %>% 
  arrange( desc(n) ) %>% 
  subset( !( species %in% c('Bare ground', 
                            'Short grass',
                            'Unknown') ) )
  
# Species counts
ks_n    <- ks_df %>% 
            count( species, quadrat, yr) %>% 
            complete( species, nesting(quadrat, yr),
                      fill = list(n = 0) )

# species area
ks_a    <- ks_df %>% 
            group_by( species, quadrat, yr) %>% 
            summarise( area = sum(area,na.rm=T) ) %>% 
            ungroup %>% 
            complete( species, nesting(quadrat, yr),
                      fill = list(n = 0) )

# test all plot/year cases are represented
plot_df <- select( ks_df, quadrat, yr) %>% 
            unique %>% 
            arrange( quadrat, yr )

# test all cases are present
select(ks_n, quadrat, yr) %>% 
  unique %>% 
  arrange(quadrat,yr) %>% 
  all.equal( plot_df ) %>% 
  expect_true


# most abundant species, BY GROUP -------------------------------  
ks_n %>%
  subset( species == 'Bouteloua gracilis' ) %>%
  subset( n != 0 ) %>%
  rename( quadrat = plot ) %>% 
  left_join( ks_q ) %>% 
  group_by( group, yr ) %>% 
  summarise( n = sum(n) ) %>% 
  ungroup %>% 
  ggplot() +
  geom_line( aes(x=yr,y=n,
                 color = group,
                 group = group) ) 

# Growth rates by N -------------------------------------

# by plot
ks_nt1     <- ks_n %>% 
                subset( !(n %in% 0) ) %>% 
                rename( nt1 = n ) %>% 
                mutate( yr = as.numeric(yr) )
  
ks_nt0     <- ks_n %>% 
                subset( !(n %in% 0) ) %>% 
                rename( nt0 = n ) %>% 
                mutate( yr = as.numeric(yr) ) %>% 
                mutate( yr = yr + 1 )
  
ks_gr_n    <- left_join( ks_nt0, ks_nt1 ) %>% 
                mutate( gr = log(nt1 / nt0) ) %>% 
                subset( !is.na(gr) ) 


# Growth rates by AREA -------------------------------------

# by plot
ks_at1     <- ks_a %>% 
                subset( !(area %in% 0 | is.na(area) ) ) %>% 
                rename( at1 = area ) %>% 
                mutate( yr = as.numeric(yr) )

ks_at0     <- ks_a %>% 
                subset( !(area %in% 0 | is.na(area)) ) %>% 
                rename( at0 = area ) %>% 
                mutate( yr = as.numeric(yr) ) %>% 
                mutate( yr = yr + 1 )

ks_gr_a    <- left_join( ks_at0, ks_at1 ) %>% 
                mutate( gr = log(at1 / at0) ) %>% 
                subset( !is.na(gr) ) 


