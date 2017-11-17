library(tidyverse)
library(stringr)
library(lubridate)

# define lons and lats
nlon   = 96
dlon   = 360 / nlon
lon    = seq(-180+dlon/2, 180-dlon/2, len=nlon)
nlat   = 48
dlat   = 180 / nlat
lat    = seq(-90+dlat/2, 90-dlat/2, len=nlat)
# 2d fields are saved as arrays in row-major order, hence expand.grid(lon, lat)
lonlat = expand.grid(lon=lon, lat=lat) %>% as_data_frame

# read topography file and create data frame
topo = readBin('../input/topography', numeric(), n=nlon*nlat, size=4)
topo_df = lonlat %>% mutate(topo=topo, type = ifelse(topo<0, 'ocean', 'land'))

# plot land/sea mask
ggplot(topo_df) + geom_raster(aes(x=lon, y=lat, fill=type))


# read wind fields (global field, twice daily for one year)
nstep_yr = 365 * 2 
ufield = readBin('../input/zonal.wind', numeric(), n=nlon*nlat*nstep_yr, size=4)
vfield = readBin('../input/meridional.wind', numeric(), n=nlon*nlat*nstep_yr, size=4)
time_df = data_frame(
  time = seq(from = as.POSIXlt('1900-01-01 00:00'), 
             to   = as.POSIXlt('1900-12-31 12:00'), 
             len  = nstep_yr))
wind_df = bind_cols(
  time_df %>% slice(rep(1:n(), each=nlon*nlat)),
  lonlat  %>% slice(rep(1:n(), nstep_yr)),
  u = ufield,
  v = vfield)

# plot only the first wind field
wind_df1 = wind_df %>% filter(yday(time)==1)
ggplot() + 
  geom_raster(data=topo_df, mapping=aes(x=lon, y=lat, fill=type), show.legend=FALSE) +
  geom_segment(data=wind_df1, 
               mapping=aes(x=lon-u/4, y=lat-v/4, xend=lon+u/4, yend=lat+v/4), 
               arrow=arrow(length=unit(0.1,'cm')))

# calculate monthly averages per grid point
wind_monavg_df = wind_df %>% 
  mutate(month=month(time)) %>%
  group_by(month, lon, lat) %>% 
  summarise(u=mean(u), v=mean(v))

# plot monthly average wind fields
ggplot(wind_monavg_df) + facet_wrap(~month, ncol=3) +
geom_raster(data=topo_df, mapping=aes(x=lon, y=lat, fill=type), show.legend=FALSE) +
geom_segment(data=wind_monavg_df, mapping=aes(x=lon-u/2, y=lat-v/2, xend=lon+u/2, yend=lat+v/2))



