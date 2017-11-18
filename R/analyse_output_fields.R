library(tidyverse)
library(stringr)
library(lubridate)
library(rnaturalearth)

# define lons and lats
nlon   = 96
dlon   = 360 / nlon
lon    = seq(dlon/2, 360-dlon/2, len=nlon)
nlat   = 48
dlat   = 180 / nlat
lat    = seq(-90+dlat/2, 90-dlat/2, len=nlat)
# 2d fields are saved as arrays in row-major order, hence expand.grid(lon, lat)
lonlat = expand.grid(lon=lon, lat=lat) %>% as_data_frame
n_grid = nlon * nlat



# Read albedo from scenario run
#
# Currently 50 years of monthly mean fields for 5 variables (Tsurf, Tair, Tocn,
# qatm, albedo) are written in order [year, month, variable, lon, lat]. Each
# saved number uses 4 bytes of memory. So output fields are written in chunks
# of size (4*n_grid) bytes. Albedo is the 5th variable, so the first albedo
# field (for year 1, month 1) starts at position (n_grid * 4), the second
# albedo field (year 1, month 2) starts at (n_grid * 4 + 5 * n_grid * 4), and
# so on. So we jump 4 * 4 * n_grid bytes using seek(), and then read 4 * n_grid
# bytes using readBin(). The data is iteratively appended to a list which
# avoids making deep copies, as opposed to appending to a vector.
scenario_file = '../output/scenario'
n_yr          = 50
scenario_con  = file(scenario_file, open='rb')
albedo        = list()
offset        = 0
for (ii in seq_len(n_yr * 12)) {
  seek(scenario_con, where = (ii-1) * 5 * 4 * n_grid + 4 * 4 * n_grid, origin='start')
  albedo[offset + 1:n_grid] = readBin(con=scenario_con, what=numeric(), n=n_grid, size=4)
  offset = offset + n_grid
}
close(scenario_con)
albedo = unlist(albedo)

# create data frame
tstamps = seq(from=as.POSIXct('1940-01-01'), by='1 month', len=n_yr*12)
albedo_df =
bind_cols(data_frame(time = tstamps) %>% slice(rep(1:n(), each=n_grid)),
          lonlat                     %>% slice(rep(1:n(), n_yr*12)),
          data_frame(albedo=albedo))
#save(file='../output/albedo.Rdata', list='albedo_df')

# plot september albedo
coast = ne_coastline() %>% fortify 
albedo_df = albedo_df %>% mutate(lon=ifelse(lon<180, lon, lon-360))
ggplot(filter(albedo_df, month(time)==9) %>% mutate(year=year(time))) +
  geom_raster(aes(x=lon, y=lat, fill=albedo)) + facet_wrap(~year) +
  geom_path(data=coast, mapping=aes(x=long, y=lat, group=id))


# plot global average albedo over time
albedo_df %>% group_by(time) %>% summarise(albedo=mean(albedo)) %>%
ggplot() + geom_line(aes(x=time, y=albedo))
