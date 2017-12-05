### Subroutines in greb.model.f90

- subroutine greb_model
- subroutine time_loop(it, isrec, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0)
- subroutine tendencies(CO2, Ts1, Ta1, To1, q1, albedo, SW, LW_surf, Q_lat, Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em)
- subroutine qflux_correction(CO2_ctrl, Ts1, Ta1, q1, To1)
- subroutine SWradiation(Tsurf, sw, albedo)
- subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
- subroutine hydro(Tsurf, q, Qlat, Qlat_air, dq_eva, dq_rain)
- subroutine seaice(Tsurf)
- subroutine deep_ocean(Ts, To, dT_ocean, dTo)
- subroutine circulation(X_in, dX_crcl, h_scl, wz)
- subroutine diffusion(T1, dX_diffuse,h_scl, wz)
- subroutine advection(T1, dX_advec,h_scl, wz)
- subroutine co2_level(it, year, CO2)
- subroutine diagnostics(it, year, CO2, ts0, ta0, to0, q0, albedo, sw, lw_surf, q_lat, q_sens)
- subroutine output(it, iunit, irec, mon, ts0, ta0, to0, q0, albedo)


