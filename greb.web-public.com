#!/bin/csh 
#
# author: Dietmar Dommenget


# compile GREB model
\rm -f greb.x *.mod
# ifort -assume byterecl -O3 -xhost -align all -fno-alias greb.model.f90 greb.shell.web-public.f90 -o greb.x 
g95 greb.model.f90 greb.shell.web-public.f90 -o greb.x 

#  namelist
cat >namelist <<EOF
&NUMERICS
time_flux = 3  ! length of flux corrections run [yrs]
time_ctrl = 3  ! length of control run [yrs
time_scnr = 50 ! length of scenariorun [yrs
/
&PHYSICS
 log_exp = 10 ! complete GREB model; 2xCO2 forcing
/
EOF

# run model
./greb.x

exit