module parameter_inc
 
implicit none
save

integer, parameter :: NRD=801, NRW=128, NCONST=256, NARRX=99, NSBMX=60, &
    NSDELOUT=39, NCNBM=60, NCNBTM=25000, NRDX=500, NTVAR=250000, &
    NTARR=250000, NEQNS=19, NBDMAX=256, NBDTMAX=1500, plot_modes=9

end module parameter_inc

module timeoutput_inc

use parameter_inc, only: NRW

implicit none

integer, parameter :: NTIMES=1024

! TOUT   - Time variables output array
! TTOUT  - time-coordinate array for time output [s] TTOUT(1:LTOUT<=NTIMES)
double precision TTOUT(NTIMES), TOUT(NTIMES, NRW), TPOUT

end module timeoutput_inc

!--------------------------------
module expdat

use parameter_inc, only: NTVAR

implicit none

double precision :: VARDAT(3,NTVAR)
integer :: INDVAR(NTVAR), IVAR

end module expdat

!--------------------------------
module ac_neg1

implicit none

integer :: NUM(4), JMIN, JMAX, NKL1, NKL2, MODK(2)

end module ac_neg1

!--------------------------------
module plasma_state ! for flight simulator plasma yes/no

implicit none

integer plasma_up, plasma_trig

end module plasma_state

!--------------------------------
module flight_sim_geometrics ! for flight simulator diagnostics

implicit none

double precision, dimension(300) :: geom1d

end module flight_sim_geometrics

!--------------------------------
module fenix_params ! for flight simulator parameters

implicit none

integer :: lhmodel, btipdirec, pr_clamp, reinitcirc, reinitpsi, &
    ipsmk2, eq_cmd, use_zlim_pot, cmnd_dioh2s, cmnd_dioh2u, nequiz, &
    s_adapt, yesfitcc, s_fazt, kastr2, com_solver, &
    isafazt, ispid_contour, &
    max_max_iteri, max_max_iterb, max_max_iterj, res_trigts06

double precision :: wallpos, d_j_m, dw_j_m, dw2_j_m, &
    vsoldiv, solwidth, &
    recycl_wall, boron_wall, predep_w, &
    transp_variance,  &
    hmodetransp, lmodetransp, &
    solmod_dt, tctr_dt, equi_dt, &
    saves_dt, savep_dt, neocl_dt, trmod_ty, nbieqmix_dt, &
    torba_dt, simdtmultip, ped_width, chie_chii, &
    D_chie, D_ped_mult, gs2d_tmin_multip, &
    timecirc, timepsi, diohdt, diohdtthreshold, dteqz2, &
    zibkdw, zifbey, dt_adapt, VV, dt_fazt, ipl_bf_bkdw, &
    resres_oh6, alp0, alpnew0, rx00, zx00
double precision, dimension(500, 2) :: psitok

end module fenix_params
