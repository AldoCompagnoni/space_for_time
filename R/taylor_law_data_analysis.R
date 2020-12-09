# Use Adler's data to look at patterns of abundance
library(tidyverse)
library(testthat)
options(stringsAsFactors = F)

# read data 
ks    <- read.csv('data/adler/kansas/allrecords.csv')
ks_q  <- read.csv('data/adler/kansas/quadrat_info.csv')


# Format --------------------------------------------------

# extract just last two digits...
last_two_digits <- function( plotyr ){
  regmatches(plotyr, 
             gregexpr("[[:digit:]]{2}$", 
                      plotyr) ) %>% unlist
}

# formatted kansas dataset
ks_df <- ks %>% 
          mutate( yr   = last_two_digits( plotyear ) ) %>% 
          mutate( nc   = nchar(plotyear)-2 ) %>% 
          mutate( plot = substr(plotyear, 1, nc) )
          
# INTERMISSION: comparison with BioTIME          
ks_df %>% 
  subset( species == 'Lesquerella ovalifolia') %>% 
  subset( yr == 72 ) %>% 
  subset( plot %in% c("e1q4-1", "e1q4-2") ) %>% 
  count( plot ) %>% 
  .$n



# Analysis ------------------------------------------------

# most abundant species
spp_abund <- ks_df %>% 
  count( species, plot, yr) %>% 
  group_by( species ) %>% 
  summarise( n = sum(n) ) %>% 
  arrange( desc(n) ) %>% 
  subset( !( species %in% c('Bare ground', 
                            'Short grass',
                            'Unknown') ) )
  
# Species counts
ks_n    <- ks_df %>% 
            count( species, plot, yr) %>% 
            complete( species, nesting(plot, yr),
                      fill = list(n = 0) )

# test all plot/year cases are represented
plot_df <- select( ks_df, plot, yr) %>% 
            unique %>% 
            arrange( plot, yr )

# test all cases are present
select(ks_n, plot, yr) %>% 
  unique %>% 
  arrange(plot,yr) %>% 
  all.equal( plot_df ) %>% 
  expect_true

# Kansas replication
ks_n %>%
  subset( n != 0 ) %>%
  ggplot( aes(yr,plot) ) +
  geom_point() +
  ggsave( 'results/taylors_law/rep_ks.tiff',
          width = 6.3, height = 6.3, compression = 'lzw') 
  

# mean by plots -------------------------------------

# means by species/plot
ks_bplot_m <- ks_n %>%
  group_by( species, plot ) %>% 
  summarise( n_mean = mean(n) ) %>% 
  ungroup %>% 
  mutate( n_log_mean = log(n_mean) )

# sd by species/plot  
ks_bplot_sd <- ks_n %>%
  group_by( species, plot ) %>% 
  summarise( n_sd = sd(n) ) %>% 
  ungroup %>% 
  mutate( n_log_sd = log(n_sd) )

# plot the log means
ks_bplot_m %>% 
  subset( species %in% spp_abund$species[1:20] ) %>% 
  ggplot( aes(n_log_mean) ) + 
  geom_histogram() +
  facet_wrap( ~ species ) +
  theme( strip.text = element_text(size=7.5) ) +
  ggsave( 'results/taylors_law/abund_by_spp_plot_ks.tiff',
          width = 6.3, height = 6.3, compression = 'lzw') 
  
# put both things together                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
ks_bplot <- left_join( ks_bplot_m, ks_bplot_sd ) %>% 
              subset( !(n_mean == 0 & n_sd == 0) ) %>% 
              subset( species %in% spp_abund$species[1:20] ) 
              

# Taylor's law using plot-level data
ggplot(ks_bplot) +
  geom_point( aes(n_log_mean, n_log_sd) ) +
  facet_wrap( ~ species ) +
  theme( strip.text = element_text(size=7.5) ) +
  ggsave( 'results/taylors_law/taylorlaw_bplot_ks.tiff',
  width = 6.3, height = 6.3, compression = 'lzw')



# Growth rates -------------------------------------

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
  
ks_gr      <- left_join( ks_nt0, ks_nt1 ) %>% 
                mutate( gr = log(nt1 / nt0) ) %>% 
                subset( !is.na(gr) ) 


# Keitt et al.'s "Growth rate Taylor's law" -----------------------\

# plot-level growth rate summaries
ks_gr_plot <- ks_gr %>% 
                group_by( species, plot ) %>% 
                summarise( gr_sd     = sd(gr),
                           gr_log_sd = sd(gr) %>% log ) %>% 
                ungroup %>% 
                subset( !is.na(gr_log_sd) ) %>% 
                left_join( ks_bplot ) %>% 
                subset( species %in% spp_abund$species[1:20] ) 

# Keitt's "growth rate Taylor's law" , by plot
ggplot(ks_gr_plot) +
  geom_point( aes(n_log_mean, gr_log_sd) ) +
  facet_wrap( ~ species ) +
  theme( strip.text = element_text(size=7.5) ) +
  ggsave( 'results/taylors_law/taylorlaw_gr_ks.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )


anova( lm(gr ~ species, data = ks_gr) ) %>% summary



# growth rate by species
ks_gr_spp <- ks_gr %>% 
  group_by( species ) %>% 
  summarise( gr_sd     = sd(gr),
             gr_log_sd = sd(gr) %>% log ) %>% 
  ungroup %>% 
  subset( !is.na(gr_log_sd) ) %>% 
  left_join( ks_m ) %>% 
  subset( species %in% spp_abund$species[1:20] ) 

# Keitt's et al. 2002 test of Taylor's law
ggplot( ks_gr_spp ) +
  geom_point( aes(n_log_mean, gr_log_sd) ) +
  ggsave( 'results/taylors_law/Keitts_taylor_test.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )

