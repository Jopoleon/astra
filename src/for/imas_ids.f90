module imas_ids

INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10)
integer, parameter, private :: DP=kind(1.0D0)

type type_coreprofile  !    Structure for core plasma profile; Time-dependent
    real(DP),pointer :: value(:) => null()     ! /value - Signal value; Time-dependent; Vector (nrho)
    character(len=132), dimension(:), pointer ::source => null()       ! /source - Source of the profile (any comment describing the origin of the profile : code, path to diagnostic s
endtype type_coreprofile

type type_profiles_1d
    real(DP),pointer :: rho_tor_norm(:) => null()     ! /profiles_1d/rho_tor - Normalised toroidal flux coordinate
    real(DP),pointer :: rho_tor(:) => null()     ! /profiles_1d/rho_tor - Toroidal flux coordinate [m]. Defined as sqrt(phi/pi/B0), where B0 = equil
    real(DP),pointer :: phi(:) => null()     ! /profiles_1d/phi - toroidal flux [Wb]; Time-dependent; Vector (npsi)
    real(DP),pointer :: psi(:) => null()     ! /profiles_1d/psi - Poloidal flux [Wb], without 1/2pi and such that Bp=|grad psi| /R/2/pi. Time-dependent; Vector (npsi)
    real(DP),pointer :: pressure(:) => null()     ! /profiles_1d/pressure - pressure profile as a function of the poloidal flux [Pa]; Time-dependent; Vector (npsi)
    real(DP),pointer :: F_dia(:) => null()     ! /profiles_1d/F_dia - diamagnetic profile (R B_phi) [T m]; Time-dependent; Vector (npsi)
    real(DP),pointer :: pprime(:) => null()     ! /profiles_1d/pprime - psi derivative of the pressure profile [Pa/Wb]; Time-dependent; Vector (npsi)
    real(DP),pointer :: ffprime(:) => null()     ! /profiles_1d/ffprime - psi derivative of F_dia multiplied with F_dia [T^2 m^2/Wb]; Time-dependent; Vector (npsi)
    real(DP),pointer :: q(:) => null()     ! /profiles_1d/q - Safety factor = dphi/dpsi [-]; Time-dependent; Vector (npsi)
    real(DP),pointer :: jparallel(:) => null()     ! /profiles_1d/jparallel - flux surface averaged parallel current density = average(j.B) / B0, where B0 = equilibrium/global_param/toroid_field/b0 ; [A/m^2];

    real(DP),pointer :: r_inboard(:) => null()     ! /profiles_1d/r_inboard - radial coordinate (major radius) at the height and on the left of the magnetic axis [m]; Time-dependent; Vector (npsi)
    real(DP),pointer :: r_outboard(:) => null()     ! /profiles_1d/r_outboard - radial coordinate (major radius) at the height and on the right of the magnetic axis [m]; Time-dependent; Vector (npsi)
    real(DP),pointer :: elongation(:) => null()     ! /profiles_1d/elongation - Elongation; Time-dependent; Vector (npsi)
    real(DP),pointer :: tria_upper(:) => null()     ! /profiles_1d/tria_upper - Upper triangularity profile; Time-dependent; Vector (npsi)
    real(DP),pointer :: tria_lower(:) => null()     ! /profiles_1d/tria_lower - Lower triangularity profile; Time-dependent; Vector (npsi)
    real(DP),pointer :: squareness(:) => null()     ! /profiles_1d/squareness - squareness; Time-dependent; Vector (npsi)
    real(DP),pointer :: volume(:) => null()     ! /profiles_1d/volume - Volume enclosed in the flux surface [m^3]; Time-dependent; Vector (npsi)
    !real(DP),pointer :: ftrap(:) => null()     ! /profiles_1d/ftrap - Trapped particle fraction; Time-dependent; Vector (npsi)
    real(DP),pointer :: fofb(:) => null() ! /profiles_1d/fofb -average(B0/B^2(l-sqrt(1-B/Bmax)(1+0.5B/Bmax)) [1]; Time-dependent; Vector (npsi)
    real(DP),pointer :: areat(:) => null() ! /profiles_1d/areat -toroidal area [m**2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: perim(:) => null() ! /profiles_1d/perim -perimeter [m]; Time-dependent; Vector (npsi)
    real(DP),pointer :: gm1(:) => null()     ! /profiles_1d/gm1 - average(1/R^2); Time-dependent; Vector (npsi)
    real(DP),pointer :: gm4(:) => null()     ! /profiles_1d/gm4 - average(1/B^2) [T^-2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: gm5(:) => null()     ! /profiles_1d/gm5 - average(B^2) [T^2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: gm41(:) => null()     ! /profiles_1d/gm41 - average(R^2) [m^2]; Time-dependent; Vector (npsi)

    real(DP),pointer :: rbp_b2(:) => null()     ! /profiles_1d/<R^2 Bp^2/B^2> - Time-dependent; Vector (npsi)
    real(DP),pointer :: bplfs(:) => null()     ! /profiles_1d/B_p on low field side - Time-dependent; Vector (npsi)

    real(DP),pointer :: acosB2a(:,:) => null()     ! /profiles_1d/<B^2 cos nt>  - Time-dependent; Vector (npsi,n)
    real(DP),pointer :: asinB2a(:,:) => null()     ! /profiles_1d/<B^2 sin nt> - Time-dependent; Vector (npsi,n)
    real(DP),pointer :: acosBlnBa(:,:) => null()     ! /profiles_1d/<BlnB cos nt> - Time-dependent; Vector (npsi,n)
    real(DP),pointer :: asinBlnBa(:,:) => null()     ! /profiles_1d/<BlnB sin nt> - Time-dependent; Vector (npsi,n)
    
    real(DP),pointer :: g1(:) => null()     ! /profiles_1d/g1 - average(grad_V^2) [m^4]; Time-dependent; Vector (npsi)
    real(DP),pointer :: g2(:) => null()     ! /profiles_1d/g2 - average(grad_V^2/R^2) [m^2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: g2int(:) => null()     ! /profiles_1d/g2int - int(grad_V^2/R^2 dV) [m^5]; Time-dependent; Vector (npsi)
    real(DP),pointer :: dPSIdV(:) => null()     ! /profiles_1d/dpsidv - dpsidv [Wb/m^3]; Time-dependent; Vector (npsi)
    real(DP),pointer :: ggradro(:) => null() ! /profiles_1d/ggradro - average(grad_V) [m^2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: gdroda(:) => null() ! /profiles_1d/gdroda - dV/dr [m^2]; Time-dependent; Vector (npsi)
    real(DP),pointer :: bmaxt(:) => null() ! /profiles_1d/bmaxt - max |B| [T]; Time-dependent; Vector (npsi)
    real(DP),pointer :: bmint(:) => null() ! /profiles_1d/bmint - min |B| [T]; Time-dependent; Vector (npsi)
    real(DP),pointer :: bdb0(:) => null() ! /profiles_1d/bdb0 - B [T]; Time-dependent; Vector (npsi)
    real(DP),pointer :: shif(:) => null()     ! /profiles_1d/shif - Shafranov shift; Time-dependent; Vector (npsi)
    real(DP),pointer :: shiv(:) => null()     ! /profiles_1d/shiv - z-center; Time-dependent; Vector (npsi)
    real(DP),pointer :: surface(:) => null()     ! /profiles_1d/surface - lateral area (m^2)
    
    type (type_coreprofile) :: jni  ! /coreprof/profiles1d/jni - [A/m^2] non-inductive parallel current density = average(jni.B) / B0, where B0 = coreprof/toroid_field/b0 
    type (type_coreprofile) :: sigmapar  ! /coreprof/profiles1d/sigmapar - Parallel conductivity [ohm^-1.m^-1].  Time-dependent. 
    type (type_coreprofile) :: te  ! /coreprof/te - Electron temperature [eV]; (source term in [W.m^-3]). Time-dependent;

endtype type_profiles_1d

type type_b0r0  !    Characteristics of the vacuum toroidal field, redundant with the toroidfield CPO, normalisation used by the ETS
    real(DP) :: r0=-9.0D40       ! /r0 - Characteristic major radius of the device (used in publications, usually middle of the vessel at the equatorial midplane) [m]. Sca
    real(DP) :: b0=-9.0D40       ! /b0 - Vacuum field at r0 [T]; Positive sign means anti-clockwise when viewed from above. Scalar. Time-dependent. 
endtype type_b0r0

type type_global_param
    type (type_b0r0) :: toroid_field  ! /global_param/toroid_field - Characteristics of the vacuum toroidal field, redundant with the toroidfield CPO, to be used by the ETS
    real(DP) :: i_plasma=-9.0D40       ! /global_param/i_plasma - total toroidal plasma current [A]; Positive sign means anti-clockwise when viewed from above. Time-dependent; Scalar
    real(DP) :: Zcurr=-9.0D40       ! /global_param/Zcurr - current centroid, time dep, scalar [m]
    real(DP) :: Rcurr=-9.0D40       ! /global_param/Zcurr - current centroid, time dep, scalar [m]
    real(DP) :: psibound=-9.0D40       ! /global_param/psibound - psi at plasma boundary alpsep [Wb]
    real(DP) :: psibound100=-9.0D40       ! /global_param/psibound - psi at plasma boundary 100% [Wb]
    real(DP) :: psiaxis=-9.0D40       ! /global_param/psiaxis - psi at mag axis [Wb]
    real(DP) :: psplex=-9.0D40       ! /global_param/psplex - avg integral at plasma boundary [?]
    real(DP) :: psiext=-9.0D40       ! /global_param/psiext - avg integral of psi external at plasma boundary [?]
    real(DP) :: Vloop=-9.0D40       ! /global_param/Vloop - loop voltage at plasma boundary [V]
    real(DP) :: li3=-9.0D40       ! /global_param/li3 - l_i defined as in ITER li3: 2 V <Bp^2>_volavg / (Rgeo*(mu0 Ip)^2)
    real(DP) :: lext=-9.0D40       ! /global_param/lext - l_ext defined as 2 Wb_ext / Ip^2
    real(DP) :: li_aug=-9.0D40       ! /global_param/li_aug - l_i defined as in AUG:  s^2 /(mu0*Ip)^2 * <Bp^2>_volavg = li3 * s^2 R0 / (2V), with s the boundary perimeter.
    real(DP) :: betpol=-9.0D40       ! /global_param/betpol - beta poloidal: Wp/Wm  , with Wp=volint(pressure) and Wm=volint(Bpol^2/(2mu0))
    real(DP) :: betpol_iter=-9.0D40       ! /global_param/betpol - beta poloidal iter definition: 4 Wp / (mu0 Rgeo Ip^2)
    real(DP) :: wkin=-9.0D40       ! /global_param/wkin - pressure energy
    real(DP) :: bpkin=-9.0D40       ! /global_param/bpkin - poloidal energy
endtype type_global_param

type type_rz2D  !    Structure for list of R,Z positions (2D)
    real(DP),pointer :: r(:,:) => null()     ! /r - Major radius [m]
    real(DP),pointer :: z(:,:) => null()     ! /z - Altitude [m]
    real(DP),pointer :: theta2d(:) => null()     ! /theta [rad]
    real(DP),pointer :: rmin(:,:) => null()     ! /r minor local [m]
    real(DP),pointer :: psirz(:,:) => null()     ! /r minor local [m]
!  real(DP),pointer :: Epol_eta(:,:) => null()     ! /poloidal electric field PS over eta [V/m * sigma]
endtype type_rz2D

type type_coord_sys
     type (type_rz2D) :: position  ! /coord_sys/position - R and Z position of grid points; Time-dependent; Matrix (ndim1, ndim2)
     
     real(DP),pointer :: gradvcell(:,:) => null() ! /coord_sys/gradvcell - |grad_V| [m^2] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
     real(DP),pointer :: bpcell(:,:) => null() ! /coord_sys/bpcell - Bp [T] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
     real(DP),pointer :: bcell(:,:) => null() ! /coord_sys/bcell - B [T] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
     real(DP),pointer :: rcell(:,:) => null() ! /coord_sys/rcell - R [m] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
     real(DP),pointer :: darea(:,:) => null() ! /coord_sys/darea - dA [m^2] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
     real(DP),pointer :: jphi(:,:) => null() ! /coord_sys/jphi - j toroidal [MA/m^2] in cells; Time-dependent; Vector (ndim1-1,ndim2-1)
      
endtype type_coord_sys

type type_rz1D_npoints  !    Structure for list of R,Z positions (1D)
    real(DP),pointer :: r(:) => null()     ! /r - Major radius [m]. Vector(max_npoints). Time-dependent
    real(DP),pointer :: z(:) => null()     ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    integer :: npoints=-999999999       ! /npoints - Number of meaningful points in the above vectors at a given time slice. Time-dependent
endtype

type type_rect_npoints  !    Structure for list of R,Z positions (1D)
    real(DP),pointer :: r2d(:) => null()     ! /r - Major radius [m]. Vector(max_npoints). Time-dependent
    real(DP),pointer :: z2d(:) => null()     ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    real(DP),pointer :: psirz2d(:,:) => null()     ! /z - Altitude [m]. Vector(max_npoints). Time-dependent
    real(DP),pointer :: fdia2d(:,:) => null()     ! F function
    integer :: npointsr=-999999999        ! /npoints - Number of meaningful points in the above vectors at a given time slice. Time-dependent
    integer :: npointsz=-999999999        ! /npoints - Number of meaningful points in the above vectors at a given time slice. Time-dependent
    double precision :: psi_axis=0.       !value of psi on axis
    double precision :: psi_boundary=0.   !value of psi on plasma boundary, alpsep
    double precision :: psi_boundary100=0.   !value of psi on plasma boundary 100%
endtype type_rect_npoints

type type_eqgeometry
    type (type_rz1D_npoints) :: boundary  ! /eqgeometry/boundary - RZ description of the plasma boundary; Time-dependent;
    type (type_rect_npoints) :: rectgrid  ! /eqgeometry/rectgrid - RZ description of the equilibrium; Time-dependent;
endtype type_eqgeometry

type type_equilibrium
    type (type_profiles_1d) :: profiles_1d  ! /equilibrium/profiles_1d
    type (type_eqgeometry) :: eqgeometry  ! /equilibrium/eqgeometry - 
    type (type_coord_sys) :: coord_sys  ! /equilibrium/coord_sys
    type (type_global_param) :: global_param  ! /equilibrium/global_param - 
    type (type_coreprofile) :: coreprofile  ! /equilibrium/global_param - 
endtype type_equilibrium

end module
