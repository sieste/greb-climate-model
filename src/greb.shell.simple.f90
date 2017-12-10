  program  time_ex

  USE mo_numerics
  USE mo_physics



! declare output fields
  real, dimension(xdim,ydim,ndays_yr) :: Tc1, Ta1, q1, ap1
  real, dimension(xdim,ydim,ndays_yr) :: Tc2, Ta2, q2, ap2

  integer, dimension(ndays_yr)::  t = (/(i,i=1,ndays_yr)/) ! jday index



100 FORMAT('climate: ',F9.2, 5E12.4)

  print*,'% start climate shell'

  ipx=46; ipy=24+8
  print*,'% diagonstic point lat/lon: ',3.75*ipy-90, 3.75*ipx
   
  open(11,file='input/tsurf',           ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(12,file='input/vapor',           ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(13,file='input/topography',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(14,file='input/soil.moisture',   ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(15,file='input/solar.radiation', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*ydim*nstep_yr)
  open(16,file='input/zonal.wind',      ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(17,file='input/meridional.wind', ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(18,file='input/ocean.mld',       ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(19,file='input/cloud.cover',     ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(20,file='input/glacier.masks',   ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)

  ! initialise modules 
  open(10,file='namelist_simple')
  call init_default_mo_physics
  call init_default_mo_numerics
  call namelist_mo_physics
  call namelist_mo_numerics


! read fix data
  read(13,rec=1)  z_topo
  read(15,rec=1)  sw_solar
  read(20,rec=1)  glacier

  do n=1,nstep_yr
     read(11,rec=n) tclim(:,:,n)
     read(12,rec=n) qclim(:,:,n)
     read(14,rec=n) swetclim(:,:,n)
     read(16,rec=n) uclim(:,:,n)
     read(17,rec=n) vclim(:,:,n)
     read(18,rec=n) mldclim(:,:,n)
     read(19,rec=n) cldclim(:,:,n)
  end do

! define deep ocean temp. as min of Tsurf but > 3.0 Celcius
  forall (i=1:xdim, j=1:ydim)
     Toclim(i,j,:) = minval(Tclim(i,j,:))
  end forall
  where (Toclim(:,:,1)-273.15 < -1.7) Toclim(:,:,1) = -1.7+273.15
  forall (i=1:xdim, j=1:ydim)
     Toclim(i,j,:) = Toclim(i,j,1)
  end forall

  print*,'% time flux/control/scenario: ', time_flux, time_ctrl, time_scnr  


  call greb_model
  
  END
