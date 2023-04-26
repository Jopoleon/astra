subroutine definitions_feqis(equil_in, j_call, ifplasma)

use pi_vars, only: GPI2, mu0
use imas_ids, only: type_equilibrium
use ef_circuit, only: data_dir, nteta, nrho, &
    psigrida, Rgeom0, teta, &
    Rexp, Zexp, Rbndp, Zbndp, &
    btor0, iplasma, pressure, pprime, ffprime, ipol, &
    psia_2d, ffp_2d, ppp_2d, ipol_2d, pres_2d

implicit none

integer, intent(in) :: j_call, ifplasma
type(type_equilibrium), intent(in) :: equil_in

integer :: i

if (j_call == 0) then
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho  = SIZE(equil_in%profiles_1d%pressure)
!normalized psi from 0 axis to 1 edge,  equispaced
    do i=1, nrho
        psia_2d(i) = (i - 1.)/(nrho - 1.)
    enddo
!constants

    Rgeom0 = equil_in%global_param%toroid_field%r0

! teta for polar grid,  goes from 0 to 2*pi-dteta,  but point nt+1 is the periodic one

    do i=1, nteta+1
        teta(i) = GPI2*(i - 1.)/(nteta + 0.)
    enddo
endif

if (ifplasma == 1) then
    btor0   = equil_in%global_param%toroid_field%b0
    iplasma = equil_in%global_param%i_plasma/1.e6
    pressure(1: nrho) = equil_in%profiles_1d%pressure(1: nrho)
    pprime(  1: nrho) = equil_in%profiles_1d%pprime(1: nrho)
    ffprime( 1: nrho) = equil_in%profiles_1d%ffprime(1: nrho)
    psigrida(1: nrho) = equil_in%profiles_1d%psi(1: nrho) !unnormalized
    psigrida(1: nrho) = (psigrida(1: nrho) - psigrida(1)) / &
                        (psigrida(nrho)    - psigrida(1)) ! normalized: 0 axis,  1 sep
    call linterp_feqis(psigrida(1: nrho), ffprime(1: nrho), nrho, &
        psia_2d(1: nrho), ffp_2d(1: nrho), nrho)
    call linterp_feqis(psigrida(1: nrho), pprime(1: nrho), nrho, &
        psia_2d(1: nrho), ppp_2d(1: nrho), nrho)
    call linterp_feqis(psigrida(1: nrho), IPOL(1: nrho), nrho, &
        psia_2d(1: nrho), ipol_2d(1: nrho), nrho)
    call linterp_feqis(psigrida(1: nrho), pressure(1: nrho), nrho, &
        psia_2d(1: nrho), pres_2d(1: nrho), nrho)
    ffp_2d = -GPI2/mu0*ffp_2d
    ppp_2d = -GPI2*1.e-6*ppp_2d

    ipol(1: nrho) = equil_in%profiles_1d%F_dia(1: nrho)

    rexp(1: nteta) = equil_in%eqgeometry%boundary%r(1: nteta)
    zexp(1: nteta) = equil_in%eqgeometry%boundary%z(1: nteta)
    rexp(nteta+1)  = rexp(1)
    zexp(nteta+1)  = zexp(1)
    rbndp(1: nteta+1) = rexp(1: nteta+1)
    zbndp(1: nteta+1) = zexp(1: nteta+1)

endif

return
end subroutine definitions_feqis
