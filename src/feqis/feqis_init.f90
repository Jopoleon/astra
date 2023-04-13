subroutine definitions_feqis(equil_in, params, j_call, ifplasma)

use pi_vars, only: GPI2, mu0
use parameters_a2spider, only: type_parameters
use imas_ids, only: type_equilibrium
use ef_circuit, only: data_dir, nteta, nrho, &
    psigrida, Rgeom0, max_iter, teta, &
    raxp, zaxp, Rexp, Zexp, Rbndp, Zbndp, &
    btor0, iplasma, pressure, pprime, ffprime, ipol, &
    psia_2d, ffp_2d, ppp_2d, ipol_2d, pres_2d
use exchange_with_astra, only: psi0_astra, psib_astra, &
    raxis_astra, zaxis_astra, nonegcurr
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_call, ifplasma
type(type_parameters) , intent(in) :: params
type(type_equilibrium), intent(in) :: equil_in

integer :: i

call markloc_ef('definitions_ef_equil', debug_lev=debug)

if (j_call == 0) then
    data_dir = params%prename(1: params%kname)
    nteta = equil_in%eqgeometry%boundary%npoints
    nrho = params%neql
!normalized psi from 0 axis to 1 edge,  equispaced
    do i=1, nrho
        psia_2d(i) = (i - 1.)/(nrho - 1.)
    enddo
!constants

    Rgeom0 = equil_in%global_param%toroid_field%r0
    max_iter = 520 !hardwired

! teta for polar grid,  goes from 0 to 2*pi-dteta,  but point nt+1 is the periodic one

    do i=1, nteta+1
        teta(i) = GPI2*(i - 1.)/(nteta + 0.)
    enddo
    psi0_astra = equil_in%profiles_1d%psi(1)
    psib_astra = equil_in%profiles_1d%psi(nrho)
    raxp = raxis_astra
    zaxp = zaxis_astra
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
    if (nonegcurr == 1) then
        do i = 1, nrho
            ffp_2d(i) = max(0., ffp_2d(i))
            ppp_2d(i) = max(0., ppp_2d(i))
        enddo
    endif

    ipol(1: nrho) = equil_in%profiles_1d%F_dia(1: nrho)

    if (params%k_fixfree == 0) then !if 1,  comes from free boundary
        rexp(1: nteta) = equil_in%eqgeometry%boundary%r(1: nteta)
        zexp(1: nteta) = equil_in%eqgeometry%boundary%z(1: nteta)
        rexp(nteta+1)  = rexp(1)
        zexp(nteta+1)  = zexp(1)
        rbndp(1: nteta+1) = rexp(1: nteta+1)
        zbndp(1: nteta+1) = zexp(1: nteta+1)
    endif
endif

return
end subroutine definitions_feqis
