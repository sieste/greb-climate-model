!
!----------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model 
!----------------------------------------------------------
!
!   Authors; Dietmar Dommenget and Janine Flöter 
!            with numerical opitmizations by Micheal Rezny
!
!   Reference: Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model
!            by Dietmar Dommenget and Janine Flöter, submitted to Climate Dynamics 2010.
!
! 
!  input fields: The GREB model needs the following fields to be specified before 
!                the main subroutine greb_model is called:
!
!   z_topo(xdim,ydim):            topography (<0 are ocean points) [m]
!  glacier(xdim,ydim):            glacier mask ( >0.5 are glacier points )
!    Tclim(xdim,ydim,nstep_yr):   mean Tsurf                       [K]
!    uclim(xdim,ydim,nstep_yr):   mean zonal wind speed            [m/s]
!    vclim(xdim,ydim,nstep_yr):   mean meridional wind speed       [m/s]
!    qclim(xdim,ydim,nstep_yr):   mean atmospheric humidity        [kg/kg]
!  mldclim(xdim,ydim,nstep_yr):   mean ocean mixed layer depth     [m]
!   Toclim(xdim,ydim,nstep_yr):   mean deep ocean temperature      [K]
! swetclim(xdim,ydim,nstep_yr):   soil wetnees, fraction of total  [0-1]
! sw_solar(ydim,nstep_yr):        24hrs mean solar radiation       [W/m^2]
!
!
!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: dt        = 12*3600           ! time step [s]
  integer, parameter :: dt_crcl   = 0.5*3600          ! time step circulation [s]  
  integer, parameter :: ndt_days  = 24*3600/dt        ! number of timesteps per day
  integer, parameter :: nstep_yr  = ndays_yr*ndt_days ! number of timesteps per year
  integer            :: time_flux = 0                 ! length of integration for flux correction [yrs]
  integer            :: time_ctrl = 0                 ! length of integration for control run  [yrs]
  integer            :: time_scnr = 0                 ! length of integration for scenario run [yrs]
  integer            :: ipx       = 1                 ! points for diagonstic print outs
  integer            :: ipy       = 1                 ! points for diagonstic print outs
  integer, parameter, dimension(12) :: jday_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per 
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat

  integer            :: ireal     = 4         ! record length for IO (machine dependent)
! 												ireal = 4 for Mac Book Pro 

  namelist / numerics / time_flux, time_ctrl, time_scnr

end module mo_numerics

!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  integer  :: log_exp = 0                ! process control logics for sens. exp.

! physical parameter (natural constants)
  parameter( pi        = 3.1416 )  
  parameter( sig       = 5.6704e-8 )     ! stefan-boltzmann constant [W/m^2/K^4]
  parameter( rho_ocean = 999.1 )         ! density of water at T=15C [kg/m^2]
  parameter( rho_land  = 2600. )         ! density of solid rock [kg/m^2]
  parameter( rho_air   = 1.2 )           ! density of air at 20C at NN 
  parameter( cp_ocean  = 4186. )         ! specific heat capacity of water at T=15C [J/kg/K]
  parameter( cp_land   = cp_ocean/4.5 )  ! specific heat capacity of dry land [J/kg/K]
  parameter( cp_air    = 1005. )         ! specific heat capacity of air      [J/kg/K]
  parameter( eps       = 1. )            ! emissivity for IR

! physical parameter (model values)
  parameter( d_ocean   = 50. )                     ! depth of ocean column [m]  
  parameter( d_land    = 2. )                      ! depth of land column  [m]
  parameter( d_air     = 5000. )                   ! depth of air column   [m]
  parameter( cap_ocean = cp_ocean*rho_ocean )      ! heat capacity 1m ocean  [J/K/m^2] 
  parameter( cap_land  = cp_land*rho_land*d_land ) ! heat capacity land   [J/K/m^2]
  parameter( cap_air   = cp_air*rho_air*d_air )    ! heat capacity air    [J/K/m^2]
  parameter( ct_sens   = 22.5 )                    ! coupling for sensible heat
  parameter( da_ice    = 0.25 )                    ! albedo diff for ice covered points
  parameter( a_no_ice  = 0.1 )                     ! albedo for non-ice covered points
  parameter( a_cloud   = 0.35 )                     ! albedo for clouds
  parameter( Tl_ice1   = 273.15-10. )              ! temperature range of land snow-albedo feedback
  parameter( Tl_ice2   = 273.15  )                 ! temperature range of land snow-albedo feedback
  parameter( To_ice1   = 273.15-7. )               ! temperature range of ocean ice-albedo feedback
  parameter( To_ice2   = 273.15-1.7 )              ! temperature range of ocean ice-albedo feedback 
  parameter( co_turb   = 5.0 )                     ! turbolent mixing to deep ocean [W/K/m^2]
  parameter( kappa     = 8e5 )                     ! atmos. diffusion coefficient [m^2/s]
  parameter( ce        = 2e-3  )                   ! laten heat transfer coefficient for ocean
  parameter( cq_latent = 2.257e6 )                 ! latent heat of condensation/evapoartion f water [J/kg]
  parameter( cq_rain   = -0.1/24./3600. )          ! decrease in air water vapor due to rain [1/s]
  parameter( z_air     = 8400. )                   ! scaling height atmos. heat, CO2
  parameter( z_vapor   = 5000. )                   ! scaling height atmos. water vapor diffusion
  parameter( r_qviwv   = 2.6736e3)                 ! regres. factor between viwv and q_air  [kg/m^3]

! parameter emisivity
  real, parameter, dimension(10) :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
&                                             0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)

! declare climate fields
  real, dimension(xdim,ydim)          ::  z_topo, glacier,z_ocean
  real, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, qclim, mldclim, Toclim, cldclim
  real, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, swetclim, dTrad
  real, dimension(ydim,nstep_yr)      ::  sw_solar

! declare constant fields
  real, dimension(xdim,ydim)          ::  cap_surf
  integer jday, ityr

! declare some program constants
  real, dimension(xdim, ydim)         :: wz_air, wz_vapor
  real, dimension(xdim,ydim,nstep_yr) :: uclim_m, uclim_p 
  real, dimension(xdim,ydim,nstep_yr) :: vclim_m, vclim_p 

  real :: t0, t1, t2

  namelist / physics / log_exp
  
end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim

 ! declare diagnostic fields
  real, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
&                                        ftmn, fqmn, amn, Tomn

! declare output fields
  real, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, apmm

end module mo_diagnostics

!+++++++++++++++++++++++++++++++++++++++
subroutine greb_model
!+++++++++++++++++++++++++++++++++++++++
!   climate model main loop

  use mo_numerics
  use mo_physics
  use mo_diagnostics
 
! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1,       &
&                               ts_ini, ta_ini, q_ini, to_ini       

  open(21,file='output/control',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(22,file='output/scenario',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

  dTrad = -0.16*Tclim -5. ! offset Tatmos-rad

! set ocean depth
  z_ocean=0
  do i=1,nstep_yr
     where(mldclim(:,:,i).gt.z_ocean) z_ocean = mldclim(:,:,i)
  end do
  z_ocean = 3.0*z_ocean

  if (log_exp ==  1) where(z_topo > 1.) z_topo = 1.0      ! sens. exp. constant topo
  if (log_exp <=  2) cldclim = 0.7                        ! sens. exp. constant cloud cover
  if (log_exp <=  3) qclim   = 0.0052                     ! sens. exp. constant water vapor
  if (log_exp <=  9) mldclim = d_ocean                    ! sens. exp. no deep ocean
  if (log_exp == 11) mldclim = d_ocean                    ! sens. exp. no deep ocean

! heat capacity global [J/K/m^2]
  where (z_topo  > 0.) cap_surf = cap_land
  where (z_topo <= 0.) cap_surf = cap_ocean*mldclim(:,:,1)

! initialize fields
  Ts_ini   = Tclim(:,:,nstep_yr)                          ! initial value temp. surf
  Ta_ini   = Ts_ini                                       ! initial value atm. temp.
  To_ini   = Toclim(:,:,nstep_yr)                         ! initial value temp. surf
  q_ini    = qclim(:,:,nstep_yr)                          ! initial value atmos water vapor

  CO2_ctrl = 340.
  if (log_exp == 12 .or. log_exp == 13 ) CO2_ctrl = 298.  ! A1B scenario

  ! define some program constants
  wz_air   = exp(-z_topo/z_air)
  wz_vapor = exp(-z_topo/z_vapor)
  where (uclim(:,:,:) >= 0.0) 
     uclim_m = uclim
     uclim_p = 0.0
  elsewhere
     uclim_m = 0.0
     uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0) 
     vclim_m = vclim
     vclim_p = 0.0
  elsewhere
     vclim_m = 0.0
     vclim_p = vclim
  end where
  
! compute Q-flux corrections
  print*,'% flux correction ', CO2_ctrl
  call qflux_correction(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini)

! test write qflux
  do irec=1, nstep_yr
     write(21,rec=irec)  TF_correct(:,:,irec)
  end do

! control run
  print*,'% CONTROL RUN CO2=',CO2_ctrl,'  time=', time_ctrl,'yr'
  Ts1 = Ts_ini; Ta1 = Ta_ini; To1 = To_ini; q1 = q_ini;                   ! initialize fields
  mon=1; year=1970; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;  
  do it=1, time_ctrl*nstep_yr                                             ! main time loop
     call time_loop(it, isrec, year, CO2_ctrl, irec, mon, 21, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0 ) 
    Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0    
  end do

! scenario run
  print*,'% SCENARIO EXP: ',log_exp,'  time=', time_scnr,'yr'
  Ts1 = Ts_ini; Ta1 = Ta_ini; q1 = q_ini; To1 = To_ini                     ! initialize fields
  year=1940; CO2=280.0; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.; 
  do it=1, time_scnr*nstep_yr                                              ! main time loop
     call co2_level(it, year, CO2)

     ! sens. exp. SST+1
     if(log_exp >= 14 .and. log_exp <= 16) CO2 = CO2_ctrl
     if(log_exp >= 14 .and. log_exp <= 16) where (z_topo < 0.0) Ts1 = Tclim(:,:,ityr)+1.0 

     call time_loop(it,isrec, year, CO2, irec, mon, 22, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0 ) 
     Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0      
     if (mod(it,nstep_yr) == 0) year=year+1
  end do

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine time_loop(it, isrec, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0)
!+++++++++++++++++++++++++++++++++++++++
! main time loop

  use mo_numerics
  use mo_physics

  real, dimension(xdim,ydim):: Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, sw,       &
&                              albedo, Q_sens, Q_lat, Q_lat_air, dq_eva,      &
&                              dq_rain, dTa_crcl, dq_crcl, dq, dT_ocean, dTo, &
&                              LW_surf, LWair_down, LWair_up, em

  jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
  ityr = mod((it-1),nstep_yr)+1           ! time step in year

  call tendencies(CO2, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat,   &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl,       &
&                    dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)
  ! surface temperature
  Ts0  = Ts1  +dT_ocean +dt*( SW +LW_surf -LWair_down +Q_lat +Q_sens +TF_correct(:,:,ityr)) / cap_surf 
  ! air temperature
  Ta0  = Ta1 +dTa_crcl +dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens )/cap_air
  ! deep ocean temperature
  To0  = To1 +dTo +ToF_correct(:,:,ityr)
  ! air water vapor
  dq = dt*(dq_eva+dq_rain) +dq_crcl + qF_correct(:,:,ityr)
  where(dq .le. -q1 ) dq = -0.9*q1 ! no negative q;  numerical stability
  q0 = q1 + dq
  ! sea ice heat capacity
  call seaice(Ts0)
  ! write output
  call output(it, ionum, irec, mon, ts0, ta0, to0, q0, albedo)
  ! diagnostics: annual means plots
  call diagonstics(it, year, CO2, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)

end subroutine time_loop

!+++++++++++++++++++++++++++++++++++++++
subroutine tendencies(CO2, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat, Q_sens, Q_lat_air, dq_eva,   & 
&                     dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts1, Ta1, To1, q1, albedo, sw, LWair_up,       &
&                               LWair_down, em, Q_sens, Q_lat, Q_lat_air,      &
&                               dq_eva, dq_rain, dTa_crcl, dq_crcl, LW_surf,   &
&                               dT_ocean, dTo

    ! SW radiation model
    call SWradiation(Ts1, sw, albedo)
    ! LW radiation model
    call LWradiation(Ts1, Ta1, q1, CO2, LW_surf, LWair_up, LWair_down, em) 
    ! sensible heat flux
    Q_sens = ct_sens*(Ta1-Ts1)
    ! hydro. model
    call hydro(Ts1, q1, Q_lat, Q_lat_air, dq_eva, dq_rain)
    ! atmos. circulation
!$omp parallel sections
!$omp section
    call circulation(Ta1, dTa_crcl, z_air, wz_air)       ! air temp 
!$omp section
    call circulation( q1,  dq_crcl, z_vapor, wz_vapor)   ! atmos water vapor
!$omp end parallel sections
    ! deep ocean interaction
    call deep_ocean(Ts1, To1, dT_ocean, dTo)

end subroutine tendencies

!+++++++++++++++++++++++++++++++++++++++
subroutine  qflux_correction(CO2_ctrl, Ts1, Ta1, q1, To1)
!+++++++++++++++++++++++++++++++++++++++
!              compute heat flux correction values

  USE mo_numerics
  USE mo_physics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1, sw, albedo,    &
&                               Q_sens, Q_lat, Q_lat_air, dq_eva, dq_rain, LW_surf,  &
&                               LWair_down, LWair_up, em, dTa_crcl, dq_crcl, dTs,    &
&                               dTa, dq, T_error, dT_ocean, dTo

! time loop
  do it=1, time_flux*ndt_days*ndays_yr
     jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
     ityr = mod((it-1),nstep_yr)+1           ! time step in year
     call tendencies(CO2_ctrl, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat,  &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl, dTa_crcl,    &
&                    dT_ocean, dTo, LWair_down, LWair_up, em)

    ! surface temperature without heat flux correction
    dTs = dt*( sw +LW_surf -LWair_down +Q_lat +Q_sens) / cap_surf
    Ts0  = Ts1 +dTs +dT_ocean
    ! air temperature
    dTa = dt*( LWair_up +LWair_down -em*LW_surf +Q_lat_air -Q_sens)/cap_air
    Ta0  = Ta1 + dTa +dTa_crcl
    ! deep ocean temperature without heat flux correction
    To0  = To1 +dTo 
   ! air water vapor without flux correction
    dq = dt*(dq_eva+dq_rain) 
    q0 = q1 +dq +dq_crcl
   ! heat flux correction Tsurf
    T_error              = Tclim(:,:,ityr) -Ts0 ! error relative to Tclim
    TF_correct(:,:,ityr) = T_error*cap_surf/dt  ! heat flux in [W/m^2]
    ! surface temperature with heat flux correction
    Ts0  = Ts1 +dTs +dT_ocean +TF_correct(:,:,ityr)*dt/ cap_surf
   ! heat flux correction deep ocean
    ToF_correct(:,:,ityr) = Toclim(:,:,ityr) -To0  ! heat flux in [K/dt]
    ! deep ocean temperature with heat flux correction
    To0  = To1 +dTo +ToF_correct(:,:,ityr)
    ! water vapor flux correction
    qF_correct(:,:,ityr) = qclim(:,:,ityr) -q0
    ! air water vapor with flux correction
    q0 = q1 + dq +dq_crcl + qF_correct(:,:,ityr)
    ! sea ice heat capacity
    call seaice(Ts0)
    ! diagnostics: annual means plots
    call diagonstics(it, 0.0, CO2_ctrl, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)
    ! memory
    Ts1=Ts0; Ta1=Ta0; q1=q0;  To1=To0; 
  end do

end subroutine qflux_correction

!+++++++++++++++++++++++++++++++++++++++
subroutine SWradiation(Tsurf, sw, albedo)
!+++++++++++++++++++++++++++++++++++++++
!    SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, sw_solar,da_ice, a_no_ice, a_cloud, z_topo  &
&                         , Tl_ice1, Tl_ice2, To_ice1, To_ice2, glacier       &
&                         , cldclim, log_exp

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, sw, albedo, a_surf, a_atmos 

! atmos albedo
  a_atmos=cldclim(:,:,ityr)*a_cloud

! surface albedo
! Land:  ice -> albedo linear function of T_surf
   where(z_topo >= 0. .and. Tsurf <= Tl_ice1) a_surf = a_no_ice+da_ice   ! ice
   where(z_topo >= 0. .and. Tsurf >= Tl_ice2) a_surf = a_no_ice          ! no ice
   where(z_topo >= 0. .and. Tsurf > Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&       a_surf = a_no_ice +da_ice*(1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))
! Ocean: ice -> albedo/heat capacity linear function of T_surf
  where(z_topo < 0. .and. Tsurf <= To_ice1) a_surf = a_no_ice+da_ice      ! ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) a_surf = a_no_ice             ! no ice
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       a_surf = a_no_ice+da_ice*(1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

! glacier -> no albedo changes
  where(glacier > 0.5) a_surf = a_no_ice+da_ice

  if (log_exp <= 5) a_surf = a_no_ice

! SW flux
  albedo=a_surf+a_atmos-a_surf*a_atmos
  forall (i=1:xdim) 
     sw(i,:)=SW_solar(:,ityr)*(1-albedo(i,:)) 
  end forall

end subroutine SWradiation


!+++++++++++++++++++++++++++++++++++++++
subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
!+++++++++++++++++++++++++++++++++++++++
! new approach with LW atmos

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: sig, eps, qclim, cldclim, z_topo, jday, ityr,         &
&                           r_qviwv, z_air, z_vapor, dTrad, p_emi, log_exp

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tair, q, LWsurf, LWair, e_co2, e_cloud,   &
&                                LWair_up, LWair_down, e_vapor, em

 
  e_co2   = exp(-z_topo/z_air)*CO2         ! CO2
  e_vapor = exp(-z_topo/z_air)*r_qviwv*q   ! water vapor
  e_cloud = cldclim(:,:,ityr)              ! clouds

  if(log_exp == 11) e_vapor = exp(-z_topo/z_air)*r_qviwv*qclim(:,:,ityr)     ! sens. exp. linear-function

! total  
  em      = p_emi(4)*log( p_emi(1)*e_co2 +p_emi(2)*e_vapor +p_emi(3) ) +p_emi(7)   &
&          +p_emi(5)*log( p_emi(1)*e_co2   +p_emi(3) )                             &
&          +p_emi(6)*log( p_emi(2)*e_vapor +p_emi(3) )
  em      = (p_emi(8)-e_cloud)/p_emi(9)*(em-p_emi(10))+p_emi(10) 
  if(log_exp == 11)  em = em +0.022/(0.15*24.)*r_qviwv*(q-qclim(:,:,ityr)) ! sens. exp. linear-function
   
  LWsurf      = -sig*Tsurf**4
  LWair_down  = -em*sig*(Tair+dTrad(:,:,ityr))**4 
  LWair_up    = LWair_down 

end subroutine LWradiation


!+++++++++++++++++++++++++++++++++++++++
subroutine hydro(Tsurf, q, Qlat, Qlat_air, dq_eva, dq_rain)
!+++++++++++++++++++++++++++++++++++++++
!    hydrological model for latent heat and water vapor

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: rho_air, uclim, vclim, z_topo, swetclim, ityr,   &
&                           ce, cq_latent, cq_rain, z_air, r_qviwv, log_exp

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, q, Qlat, Qlat_air, qs, dq_eva,        &
&                                dq_rain, abswind

  Qlat=0.; Qlat_air=0.; dq_eva=0.; dq_rain=0.
  if(log_exp <=  6 .or. log_exp == 13 .or. log_exp == 15) return

  abswind = sqrt(uclim(:,:,ityr)**2 +vclim(:,:,ityr)**2) 
  where(z_topo > 0. ) abswind = sqrt(abswind**2 +2.0**2) ! land
  where(z_topo < 0. ) abswind = sqrt(abswind**2 +3.0**2) ! ocean

! saturated humiditiy (max. air water vapor)
  qs = 3.75e-3*exp(17.08085*(Tsurf-273.15)/(Tsurf-273.15+234.175));
  qs = qs*exp(-z_topo/z_air) ! scale qs by topography
! latent heat flux surface
  Qlat    = (q-qs)*abswind*cq_latent*rho_air*ce*swetclim(:,:,ityr) 

! change in water vapor
  dq_eva  = -Qlat/cq_latent/r_qviwv  ! evaporation
  dq_rain = cq_rain*q                ! rain

! latent heat flux atmos
  Qlat_air = -dq_rain*cq_latent*r_qviwv 

end subroutine hydro

!+++++++++++++++++++++++++++++++++++++++
subroutine seaice(Tsurf)
!+++++++++++++++++++++++++++++++++++++++
!              SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, z_topo, cap_surf, cap_land, cap_ocean, &
&                           log_exp, To_ice1, To_ice2, glacier, mldclim

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf

  where(z_topo < 0. .and. Tsurf <= To_ice1) cap_surf = cap_land                    ! sea ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) cap_surf = cap_ocean*mldclim(:,:,ityr) ! open ocean
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       cap_surf = cap_land + (cap_ocean*mldclim(:,:,ityr)-cap_land)     &
&                            /(To_ice2-To_ice1)*(Tsurf-To_ice1)

  if( log_exp <= 5 ) then
     where(z_topo > 0. ) cap_surf = cap_land                     ! sea ice
     where(z_topo < 0. ) cap_surf = cap_ocean*mldclim(:,:,ityr)  ! open ocean 
  end if

! glacier -> no sea ice change
  where(glacier > 0.5) cap_surf = cap_land                       ! ice sheet

end subroutine seaice

!+++++++++++++++++++++++++++++++++++++++
subroutine deep_ocean(Ts, To, dT_ocean, dTo)
!+++++++++++++++++++++++++++++++++++++++
!              deep ocean model

  USE mo_numerics,    ONLY: xdim, ydim, nstep_yr, dt 
  USE mo_physics,     ONLY: ityr, z_topo, mldclim, log_exp, To_ice2,     &
&                           cap_ocean, co_turb, z_ocean

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts, To, dT_ocean, dTo, dmld, Tx
  dT_ocean = 0.0;  dTo     = 0.0
  if ( log_exp <= 9 .or. log_exp == 11 )   return
  if ( log_exp >= 14 .and. log_exp <= 16 ) return

  if (ityr >  1) dmld = mldclim(:,:,ityr)-mldclim(:,:,ityr-1)
  if (ityr == 1) dmld = mldclim(:,:,ityr)-mldclim(:,:,nstep_yr)

! entrainment & detrainment
  where ( z_topo < 0 .and. Ts >= To_ice2 .and. dmld < 0)     &
&       dTo      = -dmld/(z_ocean-mldclim(:,:,ityr))*(Ts-To)
  where ( z_topo < 0 .and. Ts >= To_ice2 .and. dmld > 0)     &
&       dT_ocean =  dmld/mldclim(:,:,ityr)*(To-Ts)

 c_effmix = 0.5
 dTo      = c_effmix*dTo
 dT_ocean = c_effmix*dT_ocean 

! turbulent mixing
  Tx = max(To_ice2,Ts)   
  dTo      = dTo      + dt*co_turb*(Tx-To)/(cap_ocean*(z_ocean-mldclim(:,:,ityr)))
  dT_ocean = dT_ocean + dt*co_turb*(To-Tx)/(cap_ocean*mldclim(:,:,ityr)) 

end subroutine deep_ocean

!+++++++++++++++++++++++++++++++++++++++
subroutine circulation(X_in, dX_crcl, h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
! circulation with shorter time step

  USE mo_numerics,  ONLY: xdim, ydim, dt, dt_crcl
  USE mo_physics,   ONLY: log_exp, z_vapor
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: X_in, wz
  real,                       intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_crcl

  real, dimension(xdim,ydim) :: X, dx_diffuse, dx_advec
  integer time, tt

  if(log_exp  <=  4 ) return
  if(log_exp .eq.  7 .and. h_scl .eq. z_vapor) return
  if(log_exp .eq. 16 .and. h_scl .eq. z_vapor) return

  time=max(1,nint(float(dt)/dt_crcl))

  X = X_in;
  if(log_exp .eq. 8 .and. h_scl .eq. z_vapor) then
     do tt=1, time   ! time loop circulation
        call diffusion(X, dx_diffuse, h_scl, wz)
        X = X + dx_diffuse
     end do           ! time loop
  else
     do tt=1, time   ! time loop circulation
        call diffusion(X, dx_diffuse, h_scl, wz)
        call advection(X, dx_advec, h_scl, wz)
        X = X + dx_diffuse + dx_advec
     end do           ! time loop
  end if
  dX_crcl = X - X_in

end subroutine circulation

!+++++++++++++++++++++++++++++++++++++++
subroutine diffusion(T1, dX_diffuse,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    diffusion

  USE mo_numerics,   ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,    ONLY: pi, z_topo, log_exp, kappa, z_vapor
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_diffuse

  integer :: i
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/) 
  real, dimension(ydim)      :: lat, dxlat, ccx
  real, dimension(xdim)      :: T1h, dTxh
  real, dimension(xdim,ydim) :: dTx, dTy
  
  real    :: deg, dd, dx, dy, dyy, ccy, ccx2
  integer :: j, k, km1, kp1, jm1, jm2, jm3, jp1, jp2, jp3
  integer :: time2, dtdff2, tt2

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m] 
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=kappa*dt_crcl/dyy**2
  ccx=kappa*dt_crcl/dxlat**2

     ! latitudinal   
     do k=1, ydim
        km1=k-1;  kp1=k+1
        if ( k>=2 .and. k<=ydim-1)   dTy(:,k)=ccy*(                                      & 
&                         wz(:,km1)*(T1(:,km1)-T1(:,k)) +wz(:,kp1)*(T1(:,kp1)-T1(:,k)) )
        if ( k==1 )                  dTy(:,k)=ccy*wz(:,kp1)*(-T1(:,k)+T1(:,kp1))
        if ( k==ydim )               dTy(:,k)=ccy*wz(:,km1)*(T1(:,km1)-T1(:,k))
        ! longitudinal
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 2
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = 3
           jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           do j=4, xdim-3              ! longitudinal
              jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
              ! 3.order solution: stable unitl 84degree (dx=2.5degree, a=5e5)
              dTx(j,k)=ccx(k)*(                                                           &
&               10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&               +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&               +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&               +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&               +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           end do
           j = xdim-2
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1;
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
           j = xdim
           jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
           dTx(j,k)=ccx(k)*(                                                           &
&            10*( wz(jm1,k)*(T1(jm1,k)-T1(j,k))   +wz(jp1,k)*(T1(jp1,k) -T1(j,k))   )  &
&            +4*( wz(jm2,k)*(T1(jm2,k)-T1(jm1,k)) +wz(jm1,k)*(T1(j,k)   -T1(jm1,k)) )  &
&            +4*( wz(jp1,k)*(T1(j,k)  -T1(jp1,k)) +wz(jp2,k)*(T1(jp2,k) -T1(jp1,k)) )  &
&            +1*( wz(jm3,k)*(T1(jm3,k)-T1(jm2,k)) +wz(jm2,k)*(T1(jm1,k) -T1(jm2,k)) )  &
&            +1*( wz(jp2,k)*(T1(jp1,k)-T1(jp2,k)) +wz(jp3,k)*(T1(jp3,k) -T1(jp2,k)) ) )/20.
        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(1.*dxlat(k)**2/kappa))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=kappa*dtdff2/dxlat(k)**2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
              j = 1
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = xdim; jm2 = xdim-1; jm3 = xdim-2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 2
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = xdim; jm3 = xdim-1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = 3
              jp1 = j+1; jp2 = j+2; jp3 = j+3; jm1 = j-1; jm2 = j-2; jm3 = xdim;
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)=ccx2*(                                                           &
&                      10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                      +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                      +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                      +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                      +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.

                end do           ! longitudinal
              j = xdim-2
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = j+2; jp3 = 1
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim-1
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = j+1; jp2 = 1; jp3 = 2
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
              j = xdim
              jm1 = j-1; jm2 = j-2; jm3 = j-3; jp1 = 1; jp2 = 2; jp3 = 3
              dTxh(j) = ccx2*(                                                         &
&                10*( wz(jm1,k)*(T1h(jm1)-T1h(j))   +wz(jp1,k)*(T1h(jp1) -T1h(j))   )  &
&                +4*( wz(jm2,k)*(T1h(jm2)-T1h(jm1)) +wz(jm1,k)*(T1h(j)   -T1h(jm1)) )  &
&                +4*( wz(jp1,k)*(T1h(j)  -T1h(jp1)) +wz(jp2,k)*(T1h(jp2) -T1h(jp1)) )  &
&                +1*( wz(jm3,k)*(T1h(jm3)-T1h(jm2)) +wz(jm2,k)*(T1h(jm1) -T1h(jm2)) )  &
&                +1*( wz(jp2,k)*(T1h(jp1)-T1h(jp2)) +wz(jp3,k)*(T1h(jp3) -T1h(jp2)) ) )/20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability
                T1h=T1h+dTxh
            end do               ! additional time loop
            dTx(:,k)=T1h-T1(:,k)
        end if
    end do          ! y-loop
    dX_diffuse = wz * (dTx + dTy);

end subroutine diffusion

!+++++++++++++++++++++++++++++++++++++++
subroutine advection(T1, dX_advec,h_scl, wz)
!+++++++++++++++++++++++++++++++++++++++
!    advection after DD

  USE mo_numerics, ONLY: xdim, ydim, dt, dlon, dlat, dt_crcl
  USE mo_physics,  ONLY: pi, z_topo, uclim, vclim, ityr, z_vapor, log_exp
  USE mo_physics,  ONLY: uclim_m, uclim_p, vclim_m, vclim_p
  implicit none

  real, dimension(xdim,ydim), intent(in)  :: T1, wz
  real                      , intent(in)  :: h_scl
  real, dimension(xdim,ydim), intent(out) :: dX_advec

  integer :: i
  integer, dimension(ydim):: ilat = (/(i,i=1,ydim)/)
  real, dimension(ydim) :: lat, dxlat, ccx
  real, dimension(xdim) :: T1h, dTxh
  real, dimension(xdim,ydim) :: ddx, T, dTx, dTy
  integer time2, dtdff2, tt2

  real    :: deg, dx, dy, dd, dyy, ccy, ccx2
  integer :: j, k, km1, km2, kp1, kp2, jm1, jm2, jm3, jp1, jp2, jp3
  
  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m] 
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)
  ccy=dt_crcl/dyy/2.
  ccx=dt_crcl/dxlat/2.
  
     ! latitudinal   
     k=1
     kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                     vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
     end do
     k=2
     km1=k-1; kp1=k+1; kp2=k+2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1)))          &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))           &
&                                        +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) )/3. )
     end do
     do k=3, ydim-2
        km1=k-1; kp1=k+1; km2=k-2; kp2=k+2
        do j = 1, xdim
           dTy(j,k) = ccy * (                                                     &
&                       -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))        &
&                                           +wz(j,km2)*(T1(j,k)-T1(j,km2)) )      &
&                      + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1))        &
&                                           +wz(j,kp2)*(T1(j,k)-T1(j,kp2)) ) )/3.
        end do
     end do
     k=ydim-1
     km1=k-1; kp1=k+1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) )/3.      &
&                   + vclim_p(j,k,ityr)*( wz(j,kp1)*(T1(j,k)-T1(j,kp1)) ) )
     end do
     k=ydim
     km1=k-1; km2=k-2
     do j = 1, xdim
        dTy(j,k) = ccy * (                                                        &
&                    -vclim_m(j,k,ityr)*( wz(j,km1)*(T1(j,k)-T1(j,km1))           &
&                                        +wz(j,km2)*(T1(j,k)-T1(j,km2)) ) )/3.
     end do

     ! longitudinal
     do k=1, ydim
        if ( dxlat(k) > 2.5e5) then  ! unitl 25degree
           j = 1
           jm1 = xdim; jm2 = xdim-1; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = 2
           jm1 = j-1; jm2 = xdim; jp1 = j+1; jp2 = j+2
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           do j=3, xdim-2              ! longitudinal
                jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2
                dTx(j,k)= ccx(k) * (                                                  &
&                           -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))        &
&                                               +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )      &
&                          + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))        &
&                                               +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           end do
           j = xdim-1
           jm1 = j-1; jm2 = j-2; jp1 = j+1; jp2 = 1
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.
           j = xdim
           jm1 = j-1; jm2 = j-2; jp1 = 1; jp2 = 2 
           dTx(j,k)= ccx(k) * (                                                      &
&                      -uclim_m(j,k,ityr)*( wz(jm1,k)*(T1(j,k)-T1(jm1,k))            &
&                                          +wz(jm2,k)*(T1(j,k)-T1(jm2,k)) )          &
&                     + uclim_p(j,k,ityr)*( wz(jp1,k)*(T1(j,k)-T1(jp1,k))            &
&                                          +wz(jp2,k)*(T1(j,k)-T1(jp2,k)) ) )/3.

        else  ! high resolution -> smaller time steps
            dd=max(1,nint(dt_crcl/(dxlat(k)/10.0/1.))); dtdff2=dt_crcl/dd
            time2=max(1,nint(float(dt_crcl)/float(dtdff2)))
            ccx2=dtdff2/dxlat(k)/2
            T1h=T1(:,k)
            do tt2=1, time2      ! additional time loop
                j = 1
                jm1=xdim; jm2=xdim-1; jm3=xdim-2; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 2
                jm1=j-1; jm2=xdim; jm3=xdim-1; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = 3
                jm1=j-1; jm2=j-2; jm3=xdim; jp1=j+1; jp2=j+2; jp3=j+3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                do j=4, xdim-3     ! longitudinal
                    jm1=j-1; jp1=j+1; jm2=j-2; jp2=j+2; jm3=j-3; jp3=j+3
                    dTxh(j)= ccx2 * (                                                          &
&                            -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )          &
&                                                 +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )          & 
&                                                 +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )        &
&                           + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )          &
&                                                 +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )          & 
&                                                 +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                end do           ! longitudinal
                j = xdim-2
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim-1; jp2=xdim-1; jp3=1
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim-1
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=xdim; jp2=1; jp3=2
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                j = xdim
                jm1=j-1; jm2=j-2; jm3=j-3; jp1=1; jp2=2; jp3=3
                dTxh(j)= ccx2 * (                                                              &
&                          -uclim_m(j,k,ityr)*( 10*wz(jm1,k)*(T1h(j)   - T1h(jm1) )            &
&                                               +4*wz(jm2,k)*(T1h(jm1) - T1h(jm2) )            & 
&                                               +1*wz(jm3,k)*(T1h(jm2) - T1h(jm3) ) )          &
&                         + uclim_p(j,k,ityr)*( 10*wz(jp1,k)*(T1h(j)   - T1h(jp1) )            &
&                                               +4*wz(jp2,k)*(T1h(jp1) - T1h(jp2) )            & 
&                                               +1*wz(jp3,k)*(T1h(jp2) - T1h(jp3) ) ) ) /20.
                where(dTxh .le. -T1h ) dTxh = -0.9*T1h ! no negative q;  numerical stability                
                T1h = T1h + dTxh
            end do               ! additional time loop
            dTx(:,k) = T1h - T1(:,k)
        end if
    end do          ! y-loop
    dX_advec = dTx + dTy;

end subroutine advection

!+++++++++++++++++++++++++++++++++++++++
subroutine co2_level(it, year, CO2)
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: ndays_yr, ndt_days
  USE mo_physics,     ONLY: log_exp

  CO2 = 680.
  if( log_exp .eq. 12 .or.  log_exp .eq. 13 ) then
     CO2_1950=310.;  CO2_2000=370.;  CO2_2050=520.
     if (year <= 2000.) CO2=CO2_1950 + 60./50.*(year-1950.)
     if (year  > 2000. .and. year <= 2050.) CO2=CO2_2000 + 150./50.*(year-2000.)
     if (year  > 2050. .and. year <= 2100.) CO2=CO2_2050 + 180./50.*(year-2050.)
  end if

end subroutine co2_level

!+++++++++++++++++++++++++++++++++++++++
subroutine diagonstics(it, year, CO2, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)
!+++++++++++++++++++++++++++++++++++++++
!    diagonstics plots

  USE mo_numerics,    ONLY: ndays_yr, xdim, ydim, ipx ,ipy, ndt_days, nstep_yr
  USE mo_physics,     ONLY: ityr, TF_correct, qF_correct, cap_surf, Tclim
  use mo_diagnostics

! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, sw, albedo, Q_sens, Q_lat,  LW_surf

  ! diagnostics: annual means
  tsmn=tsmn+Ts0; tamn=tamn+ta0; tomn=tomn+to0; qmn=qmn+q0; amn=amn+albedo
  swmn=swmn+sw;  lwmn=lwmn+LW_surf; qlatmn=qlatmn+q_lat; qsensmn=qsensmn+Q_sens;
  ftmn=ftmn+TF_correct(:,:,ityr); fqmn=fqmn+qF_correct(:,:,ityr);
  if ( ityr == nstep_yr ) then
     tsmn    = tsmn/nstep_yr;      tamn = tamn/nstep_yr;    tomn = tomn/nstep_yr;
     qmn     = qmn/nstep_yr;
     amn     = amn/nstep_yr;       swmn = swmn/nstep_yr;    lwmn = lwmn/nstep_yr;
     qlatmn  = qlatmn/nstep_yr; qsensmn = qsensmn/nstep_yr; ftmn = ftmn/nstep_yr;
     fqmn    = fqmn/nstep_yr;
     print *, year, sum(tsmn)/(96*48)-273.15, tsmn(48,24+3)-273.15, tsmn(16,24+14)-273.15
     tsmn=0.; tamn=0.; qmn=0.; amn=0.; swmn=0.;        ! reset annual mean values
     lwmn=0.; qlatmn=0.; qsensmn=0.; ftmn=0.; fqmn=0.; ! reset annual mean values
  end if

end subroutine diagonstics

!+++++++++++++++++++++++++++++++++++++++
subroutine output(it, iunit, irec, mon, ts0, ta0, to0, q0, albedo)
!+++++++++++++++++++++++++++++++++++++++
!    write output

  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days
  USE mo_physics,      ONLY: jday
  use mo_diagnostics,  ONLY: Tmm, Tamm, Tomm, qmm, apmm
 
  ! declare temporary fields
  real, dimension(xdim,ydim)  :: Ts0, Ta0, To0, q0, albedo

  ! diagnostics: monthly means
  Tmm=Tmm+Ts0; Tamm=Tamm+ta0; Tomm=Tomm+to0; qmm=qmm+q0; apmm=apmm+albedo
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) ) then
     ndm=jday_mon(mon)*ndt_days
     irec=irec+1; write(iunit,rec=irec)  Tmm/ndm
     irec=irec+1; write(iunit,rec=irec)  Tamm/ndm
     irec=irec+1; write(iunit,rec=irec)  Tomm/ndm
     irec=irec+1; write(iunit,rec=irec)   qmm/ndm
     irec=irec+1; write(iunit,rec=irec)  apmm/ndm
     Tmm=0.; Tamm=0.;Tomm=0.; qmm=0.; apmm=0.; 
     mon=mon+1; if (mon==13) mon=1
  end if

end subroutine output
