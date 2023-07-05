subroutine a2cdf

use parameters_a2spider, only: equil_now
use parameter_inc
use const_inc
use status_inc
use outcmn_inc, only: AWD, exp_file, equ_file
use netcdf

implicit none

! NetCDF variables

logical, parameter :: verbose=.False.
integer, parameter :: l_name=8, l_unit=25, l_desc=50

integer :: ncid, j_call=1

character(len = *), parameter :: UNIT='units'
character(len = *), parameter :: DESC='long_name'

character(len = l_name), parameter :: t_lbl='TIME', r_lbl='XRHO', rh_lbl='RHO_SURF', th_lbl='THETA'
character(len = l_unit), parameter :: t_unit='s', r_unit ='-', rh_unit ='-', th_unit='rad'
character(len = l_desc), parameter :: t_desc='Time', r_desc='rho toroidal', rh_desc='rho toroidal', th_desc='Pol. angle'

integer :: n_t, n_r, n_rh, n_th, ios, j, jid, jrho, jthe
integer :: varid(1000), t_id, r_id, rh_id, th_id, nrho_surf, nthe_surf
integer :: n_devar, n_devarx, n_const, n_delout, n_int2, n_prof, n_profx, neq_1d, neq_2d
double precision, dimension(:, :), allocatable :: tmp
character(len=120) :: f_var, f_varx, f_const, f_intern, f_intern2, f_prof, f_profx, feq_1d, feq_2d, netcdf_out

save j_call

nrho_surf = SIZE(equil_now%profiles_1d%rho_tor)
nthe_surf = SIZE(equil_now%coord_sys%position%teta2d)

allocate(tmp(nrho_surf, nthe_surf))

! NetCDF output

f_var     = TRIM(AWD) // '/main/variables.txt'
f_varx    = TRIM(AWD) // '/main/variables_x.txt'
f_const   = TRIM(AWD) // '/main/constants.txt'
f_intern  = TRIM(AWD) // '/main/internal.txt'
f_intern2 = TRIM(AWD) // '/main/intern2.txt'
f_prof    = TRIM(AWD) // '/main/profiles.txt'
f_profx   = TRIM(AWD) // '/main/profiles_x.txt'
feq_1d    = TRIM(AWD) // '/main/equil_1d.txt'
feq_2d    = TRIM(AWD) // '/main/equil_2d.txt'

write(netcdf_out, '(4A, i0, A)') TRIM(awd), '/.res/ncdf/', TRIM(exp_file), TRIM(equ_file), j_call, '.cdf'

call nfcheck( nf90_create(netcdf_out, nf90_clobber, ncid) )

! Coordinate variables (time, space)

if (verbose) then
    write(6, *) '   Dimensions...'
    write(6, *) '      rho:     ', NA1
endif
write(6, *) 'NetCDF  equil2D: ', nrho_surf, nthe_surf 

call nfcheck( nf90_def_dim(ncid, t_lbl, 1  , n_t) )
call nfcheck( nf90_def_dim(ncid, r_lbl, NA1, n_r) )
call nfcheck( nf90_def_dim(ncid, rh_lbl, nrho_surf, n_rh) )
call nfcheck( nf90_def_dim(ncid, th_lbl, nthe_surf, n_th) )
  
if (verbose) then
    write(6, *) '   Defining coordinate variables...'
endif
call nfcheck( nf90_def_var(ncid,  t_lbl, NF90_DOUBLE, (/n_t /),  t_id) )
call nfcheck( nf90_def_var(ncid,  r_lbl, NF90_DOUBLE, (/n_r /),  r_id) )
call nfcheck( nf90_def_var(ncid, rh_lbl, NF90_DOUBLE, (/n_rh/), rh_id) )
call nfcheck( nf90_def_var(ncid, th_lbl, NF90_DOUBLE, (/n_th/), th_id) )

if (verbose) then
    write(6, *) '   Assigning attributes to coordinate variables...'
endif
call nfcheck( NF90_PUT_ATT(ncid,  t_id, UNIT,  t_unit) )
call nfcheck( NF90_PUT_ATT(ncid,  r_id, UNIT,  r_unit) )
call nfcheck( NF90_PUT_ATT(ncid, rh_id, UNIT, rh_unit) )
call nfcheck( NF90_PUT_ATT(ncid, th_id, UNIT, th_unit) )
call nfcheck( NF90_PUT_ATT(ncid,  t_id, DESC,  t_desc) )
call nfcheck( NF90_PUT_ATT(ncid,  r_id, DESC,  r_desc) )
call nfcheck( NF90_PUT_ATT(ncid, rh_id, DESC, rh_desc) )
call nfcheck( NF90_PUT_ATT(ncid, th_id, DESC, th_desc) )

!---------------------------------------
! Define NetCDF variables and attributes
!---------------------------------------

! Time traces

if (verbose) write(6, *) '   Defining scalar variables...'
jid = 1
call nf90_set(ncid, jid, 1, (/n_t/), f_var, n_devar, varid)
jid = jid + n_devar
call nf90_set(ncid, jid, 1, (/n_t/), f_varx, n_devarx, varid)
jid = jid + n_devarx
call nf90_set(ncid, jid, 1, (/n_t/), f_const, n_const, varid)
jid = jid + n_const
call nf90_set(ncid, jid, 1, (/n_t/), f_intern, n_delout, varid)
jid = jid + n_delout
call nf90_set(ncid, jid, 1, (/n_t/), f_intern2, n_int2, varid)

! Profiles

if (verbose) write(6, *) '   Defining profileX variables...'

jid = jid + n_int2
call nf90_set(ncid, jid, 1, (/n_r/), f_profx, n_profx, varid)
if (verbose) write(6, *) '   Defining profile variables...'
jid = jid + n_profx
call nf90_set(ncid, jid, 1, (/n_r/), f_prof, n_prof, varid)

! Equilibrium

jid = jid + n_prof
if (verbose) then
    write(6, *) '   Assigning attributes to 1d variables...'
endif
call nf90_set(ncid, jid, 1, (/n_rh/), feq_1d, neq_1d, varid)

if (verbose) then
    write(6, *) '   Assigning attributes to 2d variables...'
endif
jid = jid + neq_1d
call nf90_set(ncid, jid, 2, (/n_th, n_rh/), feq_2d, neq_2d, varid)

call nfcheck( nf90_enddef(ncid) ) ! End define mode

!--------------------
! Writing NetCDF data
!--------------------

if (verbose) write(6, *) '   Writing grid data...'
call nfcheck( nf90_put_var(ncid,  t_id, (/TIME/)) )
call nfcheck( nf90_put_var(ncid,  r_id, XRHO(1:NA1)) )
call nfcheck( nf90_put_var(ncid, rh_id, equil_now%profiles_1d%rho_tor(:)) )
call nfcheck( nf90_put_var(ncid, th_id, equil_now%coord_sys%position%teta2d(:)) )

!------------------------------------------------
if (verbose) write(6, *) '   Writing time traces'
!------------------------------------------------
jid = 0
do j = 1, n_devar
    jid = jid + 1
    call nfcheck( nf90_put_var(ncid, varid(jid), DEVAR(j)) )
enddo
do j = 1, n_devarx
    jid = jid + 1
    call nfcheck( nf90_put_var(ncid, varid(jid), DEVARX(j)) )
enddo
do j = 1, n_const
    jid = jid + 1
    call nfcheck( nf90_put_var(ncid, varid(jid), CONSTF(j)) )
enddo
do j = 1, n_delout + 1 ! + 1 because of TIME, which is then skipped
    if (j .NE. 4) then ! Skip TIME, a coordinate
        jid = jid + 1
        call nfcheck(nf90_put_var(ncid, varid(jid), DELOUT(j)))
    endif
enddo
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TSTART) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TAU   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TAUPRP) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), HRO   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), HROX  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ALBPL ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NNCX  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QETB  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF0B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF1B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF2B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF3B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF4B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF5B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF6B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF7B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF8B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QFF9B ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QITB  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QNNB  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FTO   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FTN   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BTN   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), IPLN  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), IFBEY ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROB   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROWALL) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROC   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROCO  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RON   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROE   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROI   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO0   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO1   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO2   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO3   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO4   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO5   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO6   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO7   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO8   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RO9   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ROU   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VOLUME) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIAX ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIBO ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DFPDRB) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), GP    ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), GP2   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIFB ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIFBO) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIEXT) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSPLEX) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSIEXO) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSPLXO) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), IPLFBE) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ATREQ ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PTREQ ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RBDOT ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BBDOT ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA    ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA1   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA1N  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA1E  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA1I  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA1U  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA10  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA11  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA12  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA13  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA14  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA15  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA16  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA17  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA18  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NA19  ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NAB   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ITREQ ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NITOT ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NSTEPS) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QBEAM ) )

!-----------------------------------------
if (verbose) write(6, *) '   Writing profiles...'
!-----------------------------------------

do j = 1, n_profx
    jid = jid + 1
    call nfcheck(nf90_put_var(ncid, varid(jid), EXT(1: NA1, j)))
enddo
do j = 1, 64
    jid = jid + 1
    call nfcheck(nf90_put_var(ncid, varid(jid), CAR(1:NA1,j)))
enddo

jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), AMAIN (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), AMETR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), AREAT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), B0DB2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BDB0  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BDB02 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BMAXT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), BMINT (1:NA1)) )
jid=jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CC    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CD    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CE    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CNPAD (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CNPAP (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CNPAR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CU    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUBM  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUBS  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUECR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUFI  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUFW  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUICR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CULH  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CUTOR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), CV    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DC    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DDNEO (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DDNEOD(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DIMP1 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DIMP2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DIMP3 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DLNEO (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DLNEOD(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), DRODA (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ELON  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), EQFF  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), EQPF  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ER    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F0    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F0O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F1    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F1O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F2    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F2O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F3    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F3O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F4    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F4O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F5    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F5O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F6    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F6O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F7    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F7O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F8    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F8O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F9    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), F9O   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FOFB  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FP    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FPO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), FV    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G11   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G22   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G22E  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G33   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G33E  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G41   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G42   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G43   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G44   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), G45   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), GN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), GRADRO(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), HC(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), HE(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), IPOL  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), MRHO  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), MU    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), MV    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NALF  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NDEUT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NE    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NEO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NHE3  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NHYDR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NIBM  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NIO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NIZ1  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NIZ2  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NIZ3  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NMAIN (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NNBM1 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NNBM2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NNBM3 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), NTRIT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBEAM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBLON (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBOL1 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBOL2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBOL3 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PBPER (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PDE   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PDI   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PE    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEBM  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEECR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEFW  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEICR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEIQI (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PELH  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PELON (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PEPER (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PERIM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PETOT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PFAST (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PIBM  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PIFW  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PIICR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PITOT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PRAD  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PRES  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSXR1 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSXR2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), PSXR3 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QE    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF0   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF1   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF2   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF3   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF4   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF5   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF6   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF7   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF8   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QF9   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), QU    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RHO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RUPAR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RUPFR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), RUPYR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SCUBM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD0   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD1   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD2   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD3   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD4   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD5   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD6   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD7   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD8   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SD9   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SDN   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF0TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF1TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF2TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF3TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF4TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF5TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF6TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF7TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF8TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SF9TOT(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SGNEO (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SGNEOD(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SHEAR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SHIF  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SHIV  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SLAT  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNEBM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNIBM1(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNIBM2(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNIBM3(1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNNBM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SNTOT (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SQEPS (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SRHO  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), SXHO  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TE    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TEO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TIO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TN    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TRIA  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TTRQ  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), TTRQI (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ULON  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPAR  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPARO (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPL   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS0  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS0O (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS1  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS1O (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS2  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), UPS2O (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VIMP1 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VIMP2 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VIMP3 (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VOLUM (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VP    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VPFP  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VPOL  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VR    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VRO   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VRS   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), VTOR  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), XC    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), XI    (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), XUPAD (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), XUPAP (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), XUPAR (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZEF   (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZEF1  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZEF2  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZEF3  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZIM1  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZIM2  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZIM3  (1:NA1)) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), ZMAIN (1:NA1)) )

! Equilibrium quantities
!------------------------------------------------
if (verbose) write(6, *) '   Writing equilibrium quantities'
!------------------------------------------------

jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), equil_now%profiles_1d%phi   ) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), equil_now%profiles_1d%pprime) )
jid = jid + 1
call nfcheck( nf90_put_var(ncid, varid(jid), equil_now%profiles_1d%ffprime) )
jid = jid + 1
tmp = TRANSPOSE(equil_now%coord_sys%position%r)
if (verbose) write(*, *) 'R2D', SHAPE(tmp)
call nfcheck( nf90_put_var(ncid, varid(jid), tmp) ) 
jid = jid + 1
tmp = TRANSPOSE(equil_now%coord_sys%position%z)
call nfcheck( nf90_put_var(ncid, varid(jid), tmp) )
jid = jid + 1
tmp = TRANSPOSE(equil_now%coord_sys%position%psirz)
call nfcheck( nf90_put_var(ncid, varid(jid), tmp) )

! Close the NetCDF file.
call nfcheck( nf90_close(ncid) )
write(6, *) '   Written file ' // netcdf_out

j_call = j_call + 1

return
end subroutine a2cdf

!---------------------------------------------------------------
subroutine nf90_set(ncid, jid_in, nlen, dimid, file_in, nvars, varid)

use netcdf

implicit none

integer, parameter :: nunit=31, l_name=8, l_unit=25, l_desc=50
character(len=*), parameter :: UNIT='units'
character(len=*), parameter :: DESC='long_name'

integer, intent(in) :: ncid, jid_in, nlen
integer, intent(in), dimension(nlen) :: dimid
character(len=*), intent(in) :: file_in
integer, intent(out) :: nvars
integer, dimension(1000), intent(inout) :: varid

integer :: jid, ios
character(len=l_name) :: s_name
character(len=l_desc) :: s_desc
character(len=l_unit) :: s_unit

jid = jid_in

nvars = -1 ! Error

open(nunit, file=TRIM(file_in), iostat=ios)
if (ios /= 0) then
    write(*, *) 'File ' // TRIM(file_in) // ' not found'
    write(*, *) 'Exiting subroutine netcdf'
    return
endif

read(nunit, '(/A)')
do
    read(nunit, '(A8, X, A25, X, A50)', iostat=ios) s_name, s_unit, s_desc
    if (ios /= 0) EXIT
    if (TRIM(s_name) /= 'TIME' .and. TRIM(s_name) /= 'XRHO') then
!        write(6, '( 3(A, "%") )') s_name, s_desc, s_unit
        call nfcheck( nf90_def_var(ncid, s_name, NF90_DOUBLE, dimid, varid(jid)) )
        call nfcheck( nf90_put_att(ncid, varid(jid), UNIT, s_unit) )
        call nfcheck( nf90_put_att(ncid, varid(jid), DESC, s_desc) )
        jid = jid + 1
    endif
enddo
close(nunit)

nvars = jid - jid_in

return
end subroutine nf90_set

!---------------------------------------------------------------
SUBROUTINE NFCHECK(status)

implicit none

integer, intent(in) :: status

! if (status /= NF90_NOERR) then
if (status /= 0) then
    write(6, *) 'Error', status
    CALL exit(2)
endif

return
end subroutine NFCHECK
