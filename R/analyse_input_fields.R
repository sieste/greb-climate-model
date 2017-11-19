source('functions.R')

# read topography file (stationary global field)
topo = read_greb('../input/topography', tstamps=as.Date('1900-01-01'), 
                 varname='topo', ivar=1, nvar=1)
topo = topo %>% mutate(type = ifelse(topo<0, 'ocean', 'land'))

# plot land/sea mask
ggplot(topo) + geom_raster(aes(x=lon, y=lat, fill=type))


# read wind fields (global field, twice daily for one year)
tstamps = seq(from=as.POSIXlt('1900-01-01 00:00'), to=as.POSIXlt('1900-12-31 12:00'), len=2*365)
ufield  = read_greb(file='../input/zonal.wind', tstamps=tstamps, ivar=1, nvar=1, varname='u')
vfield  = read_greb(file='../input/meridional.wind', tstamps=tstamps, ivar=1, nvar=1, varname='v')
wind    = ufield %>% mutate(v = vfield$v)

# plot only the first wind field
ggplot() + 
  geom_raster(data=topo, mapping=aes(x=lon, y=lat, fill=type), show.legend=FALSE) +
  geom_segment(data=filter(wind, yday(time)==1), 
               mapping=aes(x=lon-u/4, y=lat-v/4, xend=lon+u/4, yend=lat+v/4), 
               arrow=arrow(length=unit(0.1,'cm')))

# calculate monthly wind speeds per grid point
wind_monavg = 
  wind %>% mutate(month=month(time)) %>%
  group_by(month, lon, lat) %>% summarise(u=mean(u), v=mean(v))

# plot monthly average wind fields
ggplot(wind_monavg) + facet_wrap(~month, ncol=3) +
  geom_raster(data=topo, mapping=aes(x=lon, y=lat, fill=type), show.legend=FALSE) +
  geom_segment(data=wind_monavg, mapping=aes(x=lon-u/2, y=lat-v/2, xend=lon+u/2, yend=lat+v/2))


