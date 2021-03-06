# Namelist parameters

Parameters can be defined in a namelist file. The namelist filename can be
provided as an argument to the `greb` executable. If not provided, the default
namelist filename is `namelist`.

The following parameters can be set, default values are shown:


## PHYSICS

```
&PHYSICS_PAR
  pi        = 3.1416   
  sig       = 5.6704e-8      ! stefan-boltzmann constant [W/m^2/K^4]
  rho_ocean = 999.1          ! density of water at T=15C [kg/m^2]
  rho_land  = 2600.          ! density of solid rock [kg/m^2]
  rho_air   = 1.2            ! density of air at 20C at NN 
  cp_ocean  = 4186.          ! specific heat capacity of water at T=15C [J/kg/K]
  cp_land   = 926.222        ! specific heat capacity of dry land [J/kg/K] 
  cp_air    = 1005.          ! specific heat capacity of air      [J/kg/K]
  eps       = 1.             ! emissivity for IR
  d_ocean   = 50.                      ! depth of ocean column [m]  
  d_land    = 2.                       ! depth of land column  [m]
  d_air     = 5000.                    ! depth of air column   [m]
  ct_sens   = 22.5                     ! coupling for sensible heat
  da_ice    = 0.25                     ! albedo diff for ice covered points
  a_no_ice  = 0.1                      ! albedo for non-ice covered points
  a_cloud   = 0.35                     ! albedo for clouds
  Tl_ice1   = 273.15-10.               ! temperature range of land snow-albedo feedback
  Tl_ice2   = 273.15                   ! temperature range of land snow-albedo feedback
  To_ice1   = 273.15-7.                ! temperature range of ocean ice-albedo feedback
  To_ice2   = 273.15-1.7               ! temperature range of ocean ice-albedo feedback 
  co_turb   = 5.0                      ! turbulent mixing to deep ocean [W/K/m^2]
  kappa     = 8e5                      ! atmos. diffusion coefficient [m^2/s]
  ce        = 2e-3                     ! latent heat transfer coefficient for ocean
  cq_latent = 2.257e6                  ! latent heat of condensation/evapoartion f water [J/kg]
  cq_rain   = -0.1/24./3600.           ! decrease in air water vapor due to rain [1/s]
  z_air     = 8400.                    ! scaling height atmos. heat, CO2
  z_vapor   = 5000.                    ! scaling height atmos. water vapor diffusion
  r_qviwv   = 2.6736e3                 ! regres. factor between viwv and q_air  [kg/m^3]
  p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028, 0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)
/
```

## NUMERICS 

```
&NUMERICS_PAR
  ipx        = 1                ! point to use for console print out
  ipy        = 1                ! point to use for console print out
  time_flux  = 0                ! length of integration for flux correction [yrs]
  time_scnr  = 0                ! length of integration for scenario run [yrs]
  year0      = 1940             ! start year of the simulation
/
```

## DIAGNOSTICS

```
&DIAGNOSTICS_PAR
  output_file = "output/scenario"  ! output file name
  ens_id      = ""                 ! a suffix appended to output_file
/
```

## CO2

The CO2 concentration pathway (`co2_ppm`) can be specified as an array of
annual CO2 concentrations in ppm, for example `co2_ppm = 680, 681, 682`. If the
length of the array `co2_ppm` is less than `time_scnr`, it is padded by
repeating the last value.

```
&CO2_PAR
  co2_flux = 298   ! co2 concentration used for flux correction
  co2_ppm  = 680   ! co2 concentration for scenario
/
```

