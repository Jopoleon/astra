module dimensions_ef_parameters

implicit none

integer, parameter :: i_dim2=300

end module dimensions_ef_parameters

!-------------------------------------------
module pi_grec_vars

implicit none

double precision, parameter :: GPI=3.141592653589793, &
    GPI2=2.*GPI, GPI4=GPI2**2.0, mu0=0.4*GPI, muvac=4.e-7*GPI

end module pi_grec_vars

!-------------------------------------------
module ef_circuit       ! declaration of minimal CPOs

use dimensions_ef_parameters, only: i_dim2

implicit none 

!generic
character*120 :: data_dir
integer :: max_iter

! coordinates:
! r, z --> rectangular grid in meters
! rho --> distance of point from magnetic axis in meters
! teta --> angle from LFS midplane which is 0, coincides with Zmag plane
! psigrid --> in poloidal flux equispaced

!grids
integer :: nr, nz, nrho, nteta, nr2, nz2

double precision :: rmin, rmax, zmin, zmax
double precision :: rho(i_dim2, i_dim2), teta(i_dim2) ! rho is defined as actual distance in meters as in astra
double precision :: rpol(i_dim2, i_dim2), zpol(i_dim2, i_dim2) ! R, Z in polar coordinates
double precision :: psia_2d(i_dim2), ffp_2d(i_dim2), ppp_2d(i_dim2), ipol_2d(i_dim2), pres_2d(i_dim2)

! potential
double precision :: psirhoteta(i_dim2, i_dim2)

! boundary and axis PBE
double precision :: rbndp(i_dim2), zbndp(i_dim2), psibndp, psiaxisp
double precision :: raxp, zaxp
double precision :: rexp(i_dim2), zexp(i_dim2)

! plasma parameters
double precision :: iplasma, btor0, rgeom0, psplex, li3, betapol
double precision :: pprime(i_dim2), ffprime(i_dim2), pressure(i_dim2), psigrida(i_dim2)     !these 3 come from astra, psi is FP of astra
! current density
double precision :: jrhoteta(i_dim2, i_dim2), ipol(i_dim2)

end module ef_circuit

!------------------------------------
module metric_coefficients_pbe

use dimensions_ef_parameters, only: i_dim2

implicit none

! metric coefficients in polar coordinates 
double precision, dimension(i_dim2, i_dim2) :: lambda2d, lambda2dp

end module metric_coefficients_pbe

!------------------------------------
module exchange_with_astra       ! declaration of minimal CPOs

integer :: nonegcurr ! nonegcurr = 0 --> no negative current allowed in plasma
double precision :: raxis_astra, zaxis_astra, psi0_astra, psib_astra

end module exchange_with_astra

!----------------------------------------
module parameters_gsef

implicit none

integer, parameter :: diagnostic_gsef=0

character(len=10), parameter :: file_fields='bfields.wr'
character(len=9), parameter :: file_eqout='outp.gsef'

integer :: max_iter
character(len=120) :: name_gsefdir = 'exp/equ/aug/' ! working directory path

end module

!-----------------
module debugger_ef

implicit none

integer :: debug
character(len=132) :: last_mark, sec_last_mark

contains

   subroutine markloc_ef(str_in, debug_lev)

   character(len=*), intent(in) :: str_in
   integer, optional, intent(in) :: debug_lev

   integer :: verbose

   if (PRESENT(debug_lev)) then
      verbose = debug_lev
   else
      verbose = debug
   endif

   SELECT CASE(verbose)
   CASE(1)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '(4X, A)') TRIM(str_in)
   CASE(2)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '(4X, A)') TRIM(last_mark)
      write(*, '(8X, A)') TRIM(str_in)
   CASE(3)
      write(*, '(/A)') 'EQEF-Trackback:'
      write(*, '( 4X, A)') TRIM(sec_last_mark)
      write(*, '( 8X, A)') TRIM(last_mark)
      write(*, '(12X, A)') TRIM(str_in)
   END SELECT

   sec_last_mark = TRIM(last_mark)
   last_mark = TRIM(str_in)

   return
   end subroutine markloc_ef

end module debugger_ef
