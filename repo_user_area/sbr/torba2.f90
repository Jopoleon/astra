! http://www2.ipp.mpg.de/~emp/index.php?page=input
!======================================================================|
subroutine TORBA2

!----------------------------------------------------------------------|
! File input for: beam curvature and width, 
! tor and pol angle
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use const_inc, only: NA1, RTOR, BTOR, TIME, ROC, SGNIP, SGNBT
use status_inc, only: TE, NE, FP, XRHO, ZEF, MU, ELON, SHif , IPOL, &
   AMETR, VOLUM, PEECR, CUECR, AREAT, rho_pol, FP_NORM
use outcmn_inc, only: AWD, nml_file
use numerical_tools, only: qinterp, integr, deriv_cde
use parameters_a2equil, only : equil_now, GP2

implicit none

logical, parameter :: dump_flag=.FALSE.
integer, parameter :: n_gy_max=30, maxint=50, maxflt=50, &
     mmax=150, nmax=150, prdim = 2*mmax+2*nmax, ndat=100000, &
     npnt=5000, ianexp=2, &
     nprofvw=25, n_interp=150

logical, dimension(n_gy_max) :: beam_on
integer :: jr, jz, i, jrho, jgy, j, n_ne, n_te, ios, n_rho
integer :: npow, ncd, ndns, nte, nshot, nastra, n_gyro
integer :: ncdroutine, lfd
integer :: n_Rrect, n_Zrect, eqdim
integer, dimension(maxint) :: intinbeam
integer, dimension(n_gy_max) :: nmod
integer :: noout, iend, kend, icnt, ibgout, neq_block
integer :: extrap_coef_cPoints , extrap_coef_cFreq

double precision :: xrtol, xatol, xstep, xtbeg, xtend, xpw0, xrmaj,  &
    xrmin, xb0, xdns, edgdns, xe1, xe2,  xte0, xteedg, xe1t, xe2t,  &
    xdel0, xdeled, xelo0, xeloed, xq0, xqedg
double precision :: xpoldeg, xtordeg, alpha, beta
double precision :: rhoresult(0: 19)
double precision, dimension(n_gy_max) :: power_gyro, freq_n, &
    xryyb, xrzzb, xwyyb, xwzzb, theta_n, phi_n, theta_t, phi_t, RR_n, ZZ_n
double precision :: floatinbeam(maxflt)
double precision, dimension(:), allocatable :: eqdata, dum_r, pfmr, dum_z, pfmz, B_Z, B_T
double precision :: prdata(prdim)
double precision :: volprofw(2*nprofvw)
double precision, dimension(6*ndat) :: t1data, t1tdata
double precision, dimension(5*ndat) :: t2data
double precision, dimension(3*npnt) :: t2ndata
double precision, dimension(NA1) :: ECR, CCD, total_int
double precision, dimension(npnt) :: ctorb, rtorb, ptorb
double precision :: Rmin, Rmax, zmin, zmax, dr, dz, drho_eq, drho_interp
double precision, dimension(:, :), allocatable :: r_surf, z_surf, B_R
double precision, dimension(:), allocatable :: rho_interp, te_interp, ne_interp

double precision :: ecrh_int, eccd_int, psi_sep, psi_axis
double precision, external :: VINT, IINT

character(len=120) :: fort_name, as_nml, time_str, pecr_file, phi_file, theta_file

NAMELIST / torbeam / n_gyro, pecr_file, theta_file, phi_file, &
           beam_on, nmod, freq_n,  &
           theta_n, phi_n, RR_n, ZZ_n, xryyb, xrzzb, xwyyb, xwzzb

!----------------------------------------------------------------------|

write(6, *) 'Calling TORBEAM...'

! Initialise TORBEAM input arrays
intinbeam   = 0
floatinbeam = 0.d0
prdata = 0.d0
extrap_coef_cPoints = 3 ! skips 3 steps
extrap_coef_cFreq   = 6 ! skips 6 steps to make it faster

n_Rrect = SIZE(equil_now%eqgeometry%rectgrid%r2d)
n_Zrect = SIZE(equil_now%eqgeometry%rectgrid%z2d)

neq_block = n_Rrect*n_Zrect
eqdim = 1 + n_Rrect + n_Zrect + 4*n_Rrect*n_Zrect
if (.not. allocated(eqdata)) then
    allocate(pfmr(n_Rrect), dum_r(n_Rrect))
    allocate(pfmz(n_Zrect), dum_z(n_Zrect))
    allocate(B_R(n_Rrect, n_Zrect))
    allocate(B_T(n_Rrect))
    allocate(B_Z(n_Rrect))
    allocate(eqdata(eqdim))
endif
! Read geometry and settings

as_nml = TRIM(awd) // TRIM(nml_file)
write(*, *) 'Reading namelist ', TRIM(as_nml)

open(57, FILE=TRIM(as_nml), delim='apostrophe')
read(57, nml=torbeam, iostat=ios)
close(57)

xrmaj = RTOR*100.

eqdata = 0.d0
eqdata(1) = FP(NA1)

write(*, '(A)') 'Acquiring eqdata'
write(*, '(A, f9.4, f9.4)') 'Sign of Ip, Bt', SGNIP, SGNBT

eqdata(2: n_Rrect+1) = equil_now%eqgeometry%rectgrid%r2d
eqdata(n_Rrect+2: n_Rrect+n_Zrect+1) = equil_now%eqgeometry%rectgrid%z2d

do jr=1, n_Rrect
    pfmz = equil_now%eqgeometry%rectgrid%psirz2d(jr, :)
    call DERIV_CDE(equil_now%eqgeometry%rectgrid%z2d, 3, pfmz, dum_z, n_Zrect)
    B_R(jr, :) = dum_z/(GP2*equil_now%eqgeometry%rectgrid%r2d(jr))
enddo

do jz=1, n_Zrect
    pfmr = equil_now%eqgeometry%rectgrid%psirz2d(:, jz)
    call DERIV_CDE(equil_now%eqgeometry%rectgrid%r2d, 3, pfmr, dum_r, n_Rrect)
    B_Z = -dum_r/(GP2*equil_now%eqgeometry%rectgrid%r2d)
    B_T = equil_now%eqgeometry%rectgrid%fdia2d(:, jz)/equil_now%eqgeometry%rectgrid%r2d
    do jr=1, n_Rrect
        j = n_Rrect + n_Zrect + jr + (jz - 1)*n_Rrect + 1
        eqdata(j)               = SGNIP*B_R(jr, jz)
        eqdata(j +   neq_block) = SGNBT*B_T(jr)
        eqdata(j + 2*neq_block) = SGNIP*B_Z(jr)
        eqdata(j + 3*neq_block) = equil_now%eqgeometry%rectgrid%psirz2d(jr, jz)
    enddo
enddo

! Profiles
write(6, *) 'Acquiring profiles', NA1

n_rho = min(n_interp, NA1)
allocate(rho_interp(n_rho))
allocate( te_interp(n_rho))
allocate( ne_interp(n_rho))
if (n_interp < NA1) then ! Interpolate on reduced space grid
    drho_interp = 1./(n_rho - 1.d0)
    rho_interp = (/ (drho_interp*(i - 1.d0), i=1, n_rho) /)
    call qinterp(rho_pol(1: NA1), NE(1: NA1), NA1, rho_interp, ne_interp, n_rho)
    call qinterp(rho_pol(1: NA1), TE(1: NA1), NA1, rho_interp, te_interp, n_rho)
else ! use original profiles
    rho_interp = rho_pol(1:n_rho)
    ne_interp  = NE(1:n_rho)
    te_interp  = TE(1:n_rho)
endif

n_ne = n_rho
prdata(1: n_ne) = rho_interp(1: n_ne)
prdata(n_ne+1: 2*n_ne) = ne_interp(1: n_ne)

n_te = n_rho
prdata(2*n_ne+1: 2*n_ne+n_te) = rho_interp(1: n_te)
prdata(2*n_ne+n_te+1: 2*n_ne+2*n_te) = te_interp(1: n_te)

! Defaults
! Beam params: ---------------------------------------------------------

ndns = 2         ! analytic (1) or interpolated (2) ne profile
nte  = 2         ! analytic (1) or interpolated (2) Te profile
npow = 1         ! power absorption on(1)/off(0)
ncd  = 1         ! current drive calc. on(1)/off(0)

ncdroutine = 2   ! 0->Curba, 1->Lin-Liu, 2->Lin-Liu+momentum conservation
noout = 0        ! 1 does not write any output
nshot = 1        ! shot number
nastra = 1       ! For ASTRA <j.B>/B0 /= <j.B>/B

intinbeam(1)  = ndns
intinbeam(2)  = nte
intinbeam(4)  = npow
intinbeam(5)  = ncd
intinbeam(6)  = ianexp
intinbeam(7)  = ncdroutine
intinbeam(8)  = nprofvw
intinbeam(9)  = noout
intinbeam(10) = 1 ! 0=weakly relativistic, 1=fully
intinbeam(11) = 3 ! #harmonics
intinbeam(12) = 1 ! git, default=1 (Farina), 0<->Westerhof
intinbeam(13) = nastra
intinbeam(14) = 0 ! 1 = Fix CD in the center, not well tested
! ODE solver params: ---------------------------------------------------

xdns   = 1.02e+14  ! electron density [cm**(-3)]
edgdns = 1.e+13    ! electron density [cm**(-3)]
xe1 = 2.000        ! exponent in the density profile
xe2 = 0.500        ! exponent in the density profile
xte0 = 25.00       ! electron temperature [keV]
xteedg = 1.000     ! electron temperature [keV]
xe1t = 2.0000      ! exponent in the temperature profile
xe2t = 2.0000      ! exponent in the temperature profile
xrtol = 1.e-07     ! required rel. error
xatol = 1.e-07     ! required abs. error
xstep = 1.         ! integration step, def = 1
xtbeg = 2.         ! tbegin
xtend = 2.         ! tende
xpw0 = 1.00        ! initial beam power [MW]

! Plasma params: -------------------------------------------------

xrmin  = AMETR(NA1)*100. ! minor radius
xb0    = BTOR            ! central toroidal field
xdel0  = SHif (1)  *100. ! Shafranov shif t on the axis
xdeled = SHif (NA1)*100. ! Shafranov shif t on the edge
xelo0  = 1.0             ! elongation on the axis
xeloed = ELON(NA1)       ! elongation on the edge
xq0    = 1./MU(2)        ! safety factor on the axis
xqedg  = 1./MU(NA1-1)    ! safety factor on the edge

floatinbeam(7)  = xdns
floatinbeam(8)  = edgdns
floatinbeam(9)  = xe1
floatinbeam(10) = xe2
floatinbeam(11) = xte0
floatinbeam(12) = xteedg
floatinbeam(13) = xe1t
floatinbeam(14) = xe2t
floatinbeam(15) = xrtol
floatinbeam(16) = xatol
floatinbeam(17) = xstep
floatinbeam(18) = xtbeg
floatinbeam(19) = xtend
floatinbeam(24) = xpw0
floatinbeam(25) = xrmaj
floatinbeam(26) = xrmin
floatinbeam(27) = xb0
floatinbeam(28) = xdel0
floatinbeam(29) = xdeled
floatinbeam(30) = xelo0
floatinbeam(31) = xeloed
floatinbeam(32) = xq0
floatinbeam(33) = xqedg
floatinbeam(34) = 1.
floatinbeam(35) = ZEF(1)

pecr_file  = TRIM(awd) // TRIM(pecr_file)
theta_file = TRIM(awd) // TRIM(theta_file)
phi_file   = TRIM(awd) // TRIM(phi_file)

call uf2dr(pecr_file, TIME, power_gyro(1:n_gyro))
power_gyro(1:n_gyro) = 1d-6*power_gyro(1:n_gyro)
call uf2dr(theta_file, TIME, theta_t(1:n_gyro))
call uf2dr(phi_file, TIME, phi_t(1:n_gyro))

if (dump_flag) then
    write(*, *) 'Writing fort.61, 62 files for TORBEAM standalone'
    write(time_str, '(f6.4)') TIME

    fort_name = 'tb_eqdata_t' // TRIM(time_str) // 's.dat'
    open(61, file=TRIM(fort_name))
    write(61, '(e13.5)') eqdata
    close(61)

    fort_name = 'tb_prdata_t' // TRIM(time_str) // 's.dat'
    open(62, file=TRIM(fort_name))
    write(62, '(e13.5)') prdata(1:2*n_ne+2*n_te)
    close(62)

    fort_name = 'tb_magn_t' // TRIM(time_str) // 's.dat'
    open(62, file=TRIM(fort_name))
    write(62, '(e13.5)') equil_now%profiles_1d%rho_tor
    write(62, '(e13.5)') equil_now%profiles_1d%F_dia
    write(62, '(e13.5)') equil_now%profiles_1d%psi
    close(62)
endif 

PEECR = 0.d0
CUECR = 0.d0

gyro_loop: do jgy=1, n_gyro

! TORBEAM only for gyrotrons with finite power 
    if (beam_on(jgy) .and. (power_gyro(jgy) > 0.02)) then
! ITER         alpha = gpi*theta_n(jgy)/180.d0
! ITER         beta = gpi*phi_n(jgy)/180.d0
! ITER         xpoldeg =  asin(cos(beta)*sin(alpha))*180.d0/gpi
! ITER         xtordeg =  -atan2(tan(beta), cos(alpha))*180.d0/gpi
        if (jgy <= 4) then
            xpoldeg = theta_n(jgy)
            xtordeg = phi_n(jgy)
        else
            xpoldeg = theta_t(jgy)
            xtordeg = phi_t(jgy)
        endif 

! Initialise TORBEAM input arrays intinbeam, floatinbeam

        intinbeam(3) = nmod(jgy)

        floatinbeam(1)  = freq_n(jgy)
        floatinbeam(2)  = xtordeg
        floatinbeam(3)  = xpoldeg
        floatinbeam(4)  = RR_n(jgy) ! beam launching R
        floatinbeam(5)  = 0.d0      ! beam launching "Y"
        floatinbeam(6)  = ZZ_n(jgy) ! beam launching Z
        floatinbeam(20) = xryyb(jgy)
        floatinbeam(21) = xrzzb(jgy)
        floatinbeam(22) = xwyyb(jgy)
        floatinbeam(23) = xwzzb(jgy)
        floatinbeam(36) = 0.d0

        if (dump_flag) then
            write(*, *) 'Writing fort. files for TORBEAM standalone'
            write(fort_name, '(A, i1, 3A)') 'tb_int_gy', jgy, '_t', TRIM(time_str), 's.dat'
            open(63, file=TRIM(fort_name))
            write(63, *) intinbeam
            close(63)

            write(fort_name, '(A, i1, 3A)') 'tb_flt_gy', jgy, '_t', TRIM(time_str), 's.dat'
            open(64, file=TRIM(fort_name))
            write(64, *) floatinbeam
            close(64)
        endif 

! Call TORBEAM library:

        write(*, '(/A, i)') 'TORBEAM call for gyrotron ', jgy
        write(6, '(A, 2f9.4)') 'Phi, theta', xtordeg, xpoldeg

        iend = 0
        kend = 0
        icnt = 0
        ibgout = 0
        t1data = 0.d0
        t2data = 0.d0
        t1tdata = 0.d0
        t2ndata = 0.d0
        rhoresult = 0.d0

! use libtorbeamB
        call beam(intinbeam, floatinbeam, n_Rrect, n_Zrect, &
            eqdata, n_ne, n_te, prdata, &
            rhoresult, iend, t1data, t1tdata, kend, t2data, t2ndata, &
            icnt, ibgout, nprofvw, volprofw)

        if (rhoresult(19) /= 0.0) then
            write(6, *) 'Error on exit', rhoresult(19)
            if (dump_flag) then
                write(fort_name, '(A, i1, 3A)') 'tb_err_gy', jgy, '_t', TRIM(time_str), 's.dat'
                open(63, file=TRIM(fort_name))
                write(63, *) rhoresult(19)
                close(63)
            endif 
        endif 

! Trajectories

        do lfd=1, npnt
            rtorb(lfd) = t2ndata(lfd)
            ptorb(lfd) = t2ndata(npnt+lfd)
            ctorb(lfd) = t2ndata(2*npnt+lfd)
        enddo
        if (rhoresult(19) /= 0.0) then
            do lfd=1, npnt
                if (ISNAN(ptorb(lfd))) then
                    write(6, *) 'P isnan at j=', lfd, rtorb(lfd)
                endif 
                if (ISNAN(ctorb(lfd))) then
                    write(6, *) 'CU isnan at j=', lfd, rtorb(lfd)
                endif 
            enddo
            write(6, *) 'P  max', maxval(abs(ptorb))
            write(6, *) 'CU max', maxval(abs(ctorb))
        endif

        ECR = 0.d0
        CCD = 0.d0

        call qinterp(rtorb, ptorb, npnt, rho_pol(1: NA1), ECR, NA1)
        call qinterp(rtorb, ctorb, npnt, rho_pol(1: NA1), CCD, NA1)

! Profiles are equidistand in rho_tor, even though they are expressed
! as function of (irregular) rho_pol, because this rho_pol grid is
! equidistand in rho_tor
        ecrh_int = VINT(ECR, ROC)
        eccd_int = IINT(CCD, ROC)
        write(6, *) 'ecr int', ecrh_int
        if (ecrh_int > 1.d-6) then
            ECR = ECR/ecrh_int
            CCD = CCD/eccd_int
        endif 
        write(*, *) 'P_gyro=', power_gyro(jgy)
        write(*, *) 'Absorption per injected MW', rhoresult(13)
        write(*, *) 'Total driven current MA per MW / total MA', &
            1.e-3*rhoresult(12), &
            1.e-3*rhoresult(13)*SGNIP*power_gyro(jgy)
        PEECR(1: NA1) = PEECR(1: NA1) + rhoresult(13)*power_gyro(jgy)*ECR
        CUECR(1: NA1) = CUECR(1: NA1) + 1.e-3*rhoresult(12)*SGNIP*power_gyro(jgy)*CCD

    endif
enddo gyro_loop

deallocate(rho_interp)
deallocate( te_interp)
deallocate( ne_interp)

write(*, *) 'Exiting torba'

return
end subroutine torba2
