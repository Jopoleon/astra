subroutine EQCTST

use astra2fbe, only: tau_circuit_feqis, tau_gseq_feqis, time_astra, solve_fix, &
    activate_coil_feqis, current_limit_feqis, force_coil, refit_mode, fast_mode, &
    use_reduce_circuit, use_limiter_astra, reconnect_circuits, execute_plasma, &
    n_equivalence, n_fourier_restab_boundary, n_of_newton_iterations, &
    new_equivalence, new_resistance, resistance_change, &
    raxis_astra, zaxis_astra, psi0_astra, psib_astra, &
    sigma_b, sigma_axis, sigma_coils, sigma_energy, &
    cur_init, dr_factor_init_astra, dz_factor_init_astra, psplex_from_fbe
use const_inc, only: tau, time, RTOR, shift, psiax, psibo, iplx, taumin, ncnb, NA1
use status_inc, only: MU, SHIV
use outcmn_inc, only: machine, ccoil, ncnbt

implicit none

logical :: file_existence
integer :: ictrl, i, j, diagzz, ii_pl, jj_pl, n_p_grid, n_limz, &
    i_lim_min, j_counta, j_init, j_currr, j_rampdown, j_readgrid, &
    ni_p, nj_p

real*8, dimension(257) :: r_c_e, z_c_e 
real*8, dimension(257, 257, 130) :: m_c_e

double precision :: r_mag, z_mag, icur, dteqz, &
    r00, r01, r10, r11, z00, z01, z10, z11, &
    i00, i11, dt_eqnew, dtmaxe, ztt0, ztt1, dztt1, &
    rtt0, rtt1, drtt1, err_given, rtt2, ztt2, drtt2, &
    dztt2, err_checka, ioh1, ioh2, li3r, psiprima, &
    a_min_b, a_ratio, r_pl, z_pl, tbkdw, LL, RR, &
    ipl_threshold, dt_afazt, err_checkrzp, psplexold, &
    dist_lim_min, v_95_pos
double precision, dimension(5) :: yfircs, jumbo
double precision, dimension(12) :: icc1, icc2, dccdt
double precision, dimension(30) :: t_p_grid, rho_p_grid
double precision, dimension(130) :: conduc_cur
double precision, dimension(200) :: r_lim, z_lim, dist_lim
double precision, dimension(30, 30) :: r_p_grid, z_p_grid, ds_p_grid
double precision, dimension(257, 257) :: hoop_force, dhb_f, psi_c_e, &
    b_r_e, b_z_e
character(len=80) :: fname

data r00/0./
data j_init/0/
data r01/0./
data dt_eqnew/1.e-6/
data ioh2/0./
data psiprima/0./
data j_rampdown/0/
data j_readgrid/0/
data j_counta/0/

save r00,r01,r10,r11,a_min_b,j_init
save z10,dt_eqnew,dt_afazt
save i00,i11,psplexold
save ztt0,ztt1,dztt1,ztt2,dztt2
save rtt0,rtt1,drtt1,rtt2,drtt2
save ioh1,ioh2,icc1,icc2,dccdt
save psiprima,r_pl,z_pl
save tbkdw,j_currr,j_rampdown,j_readgrid
save m_c_e,r_c_e,z_c_e,r_lim,z_lim
save n_limz,j_counta

tau_circuit_feqis = tau
tau_gseq_feqis    = tau
time_astra     = time
activate_coil_feqis = 1 ! if 0, coil is disconnected if use_reduce_circuit and reconnect_circuits is used, otherwise just sets the current to zero (you will get a different result!)
current_limit_feqis(:, 1) =  1.e6 ! 1 is upper, 2 is lower
current_limit_feqis(:, 2) = -1.e6 ! 1 is upper, 2 is lower
force_coil = 0. ! where it is 1, forces coil i,i to current of i,j --> better use reconnect circuits.

use_reduce_circuit = 0 ! run with reconnected circuits, does not reset the matrix
reconnect_circuits = 0 ! this resets the circuit matrix, to change connections
n_equivalence = 0
new_equivalence = 0
!new_equivalence(1:14,1) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
resistance_change = 0 ! this resets the circuit matrix, to change resistances (diagonals)
new_resistance=0.
!	new_resistance(1:12) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) !diagonal resistance in microOhm

raxis_astra = RTOR + SHIFT
zaxis_astra = SHIV(1)
psi0_astra = PSIAX
psib_astra = PSIBO
n_fourier_restab_boundary = 5
use_limiter_astra = 1   ! do not use limiter for DEMO
refit_mode = 0   ! if -1 - 1 turn only, 0 - stab method, if 1 - restab with prescribed axis , 3 - full fit like spider but only for eddy currents, 101 - only Z stab
sigma_B = 1.
sigma_axis = 50000.
sigma_coils = 1.e5
sigma_energy = 1

solve_fix = 0   ! if 0 - solve full fix boundary problem, >0 - N pass only, -2 - uses fbe solution 
execute_plasma = 1   ! if 0 - only circuit equations, if 1 - solve plasma gseq too
n_of_newton_iterations = 150

!factors of dr and dz for initial iterations
dr_factor_init_astra = 1.
dz_factor_init_astra = 1.
fast_mode = 0
psplex_from_fbe = 0
ipl_threshold = IPLX !Ip at which astra starts to work

dt_eqnew = taumin
dteqz  = dt_eqnew ! spider GS solver time step

execute_plasma = 1
tau_gseq_feqis = dteqz

if (MACHINE(1:3) == 'aug') then
    cur_init( 1:12) = CCOIL(1:12)/1.e3
    cur_init(13:52) = 0.
endif

if (TIME > 2.52) fast_mode = 1

write(*,*) 'eqtime', time, fast_mode, execute_plasma, tau_gseq_feqis, &
    ncnb, ncnbt, cur_init(1:12)
write(*,*) V_95_POS(1./mu(1:na1))

return
end subroutine eqctst
