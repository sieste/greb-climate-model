library(tidyverse)
library(stringr)
library(lubridate)
library(rnaturalearth)

coast = ne_coastline() %>% fortify 


#' Read a variable from greb model output
#'
#' @param file filename
#' @param tstamps time stamps of the outputs
#' @param varname variable name
#' @param ivar variable index if there are several variables
#' @param nvar number of variables in output file
#' @param nlon number of longitudes
#' @param nlat number of latitudes
#' @param nbyte number of bytes per datum
#' @return a tidy data_frame with columns time, lon, lat, and the requested variable
#'
#' @details
#' 2d fields are saved as arrays in row-major order, so we use
#' `expand.grid(lon, lat)` instead of `(lat, lon)`.
#'
#' In the default greb output, 50 years of monthly mean fields for 5 variables
#' (Tsurf, Tair, Tocn, qatm, albedo) are written in order [year, month,
#' variable, lon, lat]. Each saved number uses 4 bytes of memory. So output
#' fields are written in chunks of size (4*n_grid) bytes. For example, albedo is
#' the 5th variable, so the first albedo field (for year 1, month 1) starts at
#' position (n_grid * 4), the second albedo field (year 1, month 2) starts at
#' (n_grid * 4 + 5 * n_grid * 4), and so on. So we jump 4 * 4 * n_grid bytes
#' using seek(), and then read 4 * n_grid bytes using readBin(). The data is
#' iteratively appended to a list which avoids making deep copies, as opposed to
#' appending to a vector.
# 
read_greb = function(file, tstamps, 
                     varname=str_c('variable_',ivar), 
                     ivar=1L, nvar=5L, 
                     nlon=96, nlat=48, nbyte=4) {

  stopifnot(file.exists(file))
  stopifnot(file.size(file) == nlon * nlat * nvar * nbyte * length(tstamps))
  # latitude/longitude grid (longitude varies fastest)
  ngrid = nlon * nlat
  dlon   = 360 / nlon
  dlat   = 180 / nlat
  lon    = seq(dlon/2, 360-dlon/2, len=nlon)
  lat    = seq(-90+dlat/2, 90-dlat/2, len=nlat)
  lonlat = expand.grid(lon=lon, lat=lat) %>% as_data_frame

  ntime = length(tstamps)

  # read requested data from file
  con  = file(file, open='rb')
  out  = list()
  nout = 0
  for (ii in seq_len(ntime)) {
    seek(con, where = nbyte * ngrid * ((ii-1) * nvar + (ivar - 1)))
    out[nout + 1:ngrid] = readBin(con=con, what=numeric(), n=ngrid, size=nbyte)
    nout = nout + ngrid
  }
  close(con)
  out = list(unlist(out)) %>% setNames(varname) %>% as_data_frame

  # generate data frame with time, lon, lat, and value
  out_df = bind_cols(
    data_frame(time = tstamps) %>% slice(rep(1:ntime, each=ngrid)),
    lonlat                     %>% slice(rep(1:ngrid, ntime)),
    out)

  return(out_df) 
} 


#' Transform longitudes between ranges [0,360] and [-180, 180]
#'
#' @param df data frame with column `lon` or `long`
#' @param how character to which range to transform ('0_360' or '-180_180')
#' @return the original data frame with transformed longitude column
wrap_lon = function(df, how=c('0_360', '-180_180')) {
  how = match.arg(how)
  is_long = FALSE
  if ('long' %in% names(df)) {
    names(df) = names(df) %>% str_replace('^long$', 'lon')
    is_long = TRUE
  }
  if (how == '0_360') {
    df = df %>% mutate(lon = ifelse(lon < 0, lon + 360, lon))
  }
  if (how == '-180_180') {
    df = df %>% mutate(lon = ifelse(lon > 180, lon - 360, lon))
  }
  if (is_long) {
    names(df) = names(df) %>% str_replace('^lon$', 'long')
  }
  return(df)
}
