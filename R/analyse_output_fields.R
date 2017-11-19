library(tidyverse)
library(stringr)
library(lubridate)
library(rnaturalearth)

source('functions.R')

scenario_file = '../output/scenario'
tstamps = seq(from=as.Date('1940-01-01'), by='1 month', len=50*12)
albedo_df = read_greb(file=scenario_file, tstamps=tstamps, varname='albedo', 
                      ivar=5, nvar=5, nlon=96, nlat=48) 


# plot september albedo
coast = ne_coastline() %>% fortify 
albedo_df = albedo_df %>% mutate(lon=ifelse(lon<180, lon, lon-360))
ggplot(filter(albedo_df, month(time)==9) %>% mutate(year=year(time))) +
  geom_raster(aes(x=lon, y=lat, fill=albedo)) + facet_wrap(~year) +
  geom_path(data=coast, mapping=aes(x=long, y=lat, group=id))


# plot global average albedo over time
albedo_df %>% group_by(time) %>% summarise(albedo=mean(albedo)) %>%
ggplot() + geom_line(aes(x=time, y=albedo))
