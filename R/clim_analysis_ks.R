# Use Adler's data to look at patterns of abundance
library(tidyverse)
library(testthat)
library(plotly)
library(lme4)
options(stringsAsFactors = F)

# read data 
ks        <- read.csv('data/adler/kansas/allrecords.csv')
ks_q      <- read.csv('data/adler/kansas/quadrat_info.csv') %>% 
              dplyr::select(-shapefiles,-quadX,-quadY)
                                                                          
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


# most abundant species, BY GROUP, ONLY BOUTEOULA GRACILIS -------------------------------  
p_bogr <- ks_n %>%
  subset( species == 'Bouteloua gracilis' ) %>%
  subset( n != 0 ) %>%
  left_join( ks_q ) %>% 
  group_by( group, yr ) %>% 
  summarise( n = sum(n) ) %>% 
  ungroup %>% 
  ggplot() +
  geom_line( aes(x=yr,y=n,
                 color = group,
                 group = group) ) 

# explore interactively (exclude groups with little BOGR)
ggplotly(p_bogr)



# Growth rates by N -------------------------------------

# by plot
ks_nt1     <- ks_n %>% 
                left_join( ks_q ) %>% 
                subset( !(n %in% 0) ) %>% 
                rename( nt1 = n ) %>% 
                mutate( yr = as.numeric(yr) )
  
ks_nt0     <- ks_n %>% 
                left_join( ks_q ) %>% 
                subset( !(n %in% 0) ) %>% 
                rename( nt0 = n ) %>% 
                mutate( yr = as.numeric(yr) ) %>% 
                mutate( yr = yr + 1 )
  
ks_gr_n    <- left_join( ks_nt0, ks_nt1 ) %>% 
                mutate( gr = log(nt1 / nt0) ) %>% 
                subset( !is.na(gr) ) %>% 
                left_join( rename(prec_anom, yr=MatrixEndYear) )

# by group
ks_group_n    <- left_join( ks_nt0, ks_nt1 ) %>% 
  group_by(species, group, yr) %>%
  summarise(nt1 = sum(nt1),
            nt0 = sum(nt0)) %>%
  mutate( gr = log(nt1 / nt0) ) %>% 
  subset( !is.na(gr) ) %>% 
  left_join( rename(prec_anom, yr=MatrixEndYear) )


# Growth rates by AREA -------------------------------------

# by plot
ks_at1     <- ks_a %>% 
                left_join( ks_q ) %>% 
                subset( !(area %in% 0 | is.na(area) ) ) %>% 
                rename( at1 = area ) %>% 
                mutate( yr = as.numeric(yr) )

ks_at0     <- ks_a %>% 
                left_join( ks_q ) %>% 
                subset( !(area %in% 0 | is.na(area)) ) %>% 
                rename( at0 = area ) %>% 
                mutate( yr = as.numeric(yr) ) %>% 
                mutate( yr = yr + 1 )

ks_gr_a    <- left_join( ks_at0, ks_at1 ) %>% 
                mutate( gr = log(at1 / at0) ) %>% 
                subset( !is.na(gr) ) %>% 
                left_join( rename(prec_anom, yr=MatrixEndYear) )

# by group
ks_group_a    <- left_join( ks_at0, ks_at1 ) %>%
  group_by(species, group, yr) %>%
  summarise(at1 = sum(at1, na.rm = T),
            at0 = sum(at0, na.rm = T)) %>%
  mutate( gr = log(at1 / at0) ) %>% 
  subset( !is.na(gr) & is.finite(gr) ) %>% 
  left_join( rename(prec_anom, yr=MatrixEndYear) )


# Only analyze bouteoula gracilis ----------------------------------------

# Ideal functional form 
# gr ~ ppt_t0|group + nt0 versus gr ~ ppt_t0|group
# gr ~ ppt_t0|group + at0 versus gr ~ ppt_t0|group

# also, could we pick 

# correlations are strong
cor_mat <- ks_gr_n %>%
  subset( species == 'Bouteloua gracilis' ) %>%
  group_by( yr, group ) %>% 
  summarise( gr = mean(gr,na.rm=T) ) %>% 
  ungroup %>% 
  # MOST ABUNDANT groups. I've picked these out using ggplotly
  subset( group %in% c('et1', 'et2', 'et3', 'lb2', 'lb3',
                       'sg1', 'sg1a', 'sg2', 'sg3') ) %>% 
  pivot_wider( names_from = group,
               values_from = gr ) %>% 
  dplyr::select( -yr) %>% 
  as.matrix %>% 
  cor(use='na.or.complete') 
  

cor_mat[upper.tri(cor_mat)] %>% as.numeric %>% hist

  # subset( n != 0 ) %>%
  ggplot() +
  geom_line( aes(x=yr,y=gr,
                 color = group,
                 group = group) ) 

# estimate climatic effects  
mod <- lmer( gr ~ at0 + ppt_t0 + (ppt_t0|group) , data= ks_gr_a) 
mod <- lmer( gr ~ ppt_t0 + (0 + ppt_t0|group), data= ks_gr_a)

mod <- lmer( gr ~ nt0 + ppt_t0 + (ppt_t0|group) , data= ks_gr_n) 
mod <- lmer( gr ~ ppt_t0 + (ppt_t0|group) , data= ks_gr_n) 

ks_gr_a_group <- ks_gr_a %>% 
  group_by(yr,group) %>% 
  summarise( gr = mean(gr,na.rm=T) ) %>% 
  ungroup  %>% 
  left_join( rename(prec_anom, yr=MatrixEndYear) )

ks_gr_n_group <- ks_gr_n %>% 
  group_by(yr,group) %>% 
  summarise( gr = mean(gr,na.rm=T) ) %>% 
  ungroup  %>% 
  left_join( rename(prec_anom, yr=MatrixEndYear) )

ggplot(ks_gr_a) +
  geom_point( aes(ppt_t0, gr) )


# estimate climatic effects 10 years at a time -------------------

# 
yr_v <- ks_gr_a_group$yr %>% unique %>% sort

# fit linear model to a subset of years
lm_subset <- function( ii ){
  
  lm(gr ~ ppt_t0,
     data=subset(ks_gr_a, yr %in% yr_v[(c(1:10)+(ii-1))]) ) %>% 
    coef %>% 
    .[2]
   
}

# retrieve betas
betas <- sapply(1:29, lm_subset)



# Variance in growth rates per group  -------------------

# by numbers
gr_n <- ks_gr_n %>%
  filter(species == "Bouteloua gracilis") %>%
  group_by(group) %>%
  summarise(sd_gr = sd(gr),
            mean_abun = mean(nt0))

plot(gr_n$mean_abun, gr_n$sd_gr, xlab = "mean abundance per group (in N)", ylab = "sd")


# by area
gr_a <- ks_gr_a %>%
  filter(species == "Bouteloua gracilis") %>%
  group_by(group) %>%
  summarise(sd_gr = sd(gr),
            mean_abun = mean(at0))

plot(gr_a$mean_abun, gr_a$sd_gr, xlab = "mean abundance per group (in A)", ylab = "sd")
