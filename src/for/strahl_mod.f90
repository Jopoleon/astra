module strahl_mod

use parameter_inc, only: NRD

implicit none

double precision, dimension(NRD) :: prad_tot, ne_source, nneut_imp
double precision, dimension(NRD, 11) :: prad_strahl, nimp_strahl, nesrc_strahl

contains

!----------------------------------------------------------------------
subroutine profiles_file_write_strahl(strahl_dir, rho_coord, rhopol, &
    ngrid, ne, te, ti, ne_decayl, te_decayl, ti_decayl, time, &
    rhopolg, teg, neg, tig, ngmax)

implicit none

integer, parameter :: nch_w2=42
character(len=14), parameter :: strahl_prof_in='nete/pp11111.1'

integer, intent(in) :: ngrid, ngmax
double precision, intent(in) :: time,  ne_decayl, te_decayl, ti_decayl
double precision, intent(in), dimension(ngrid) :: ne, te, ti, rhopol
double precision, intent(in), dimension(ngmax) :: rhopolg, teg, tig, neg
character(len=20), intent(in) :: rho_coord
character(len=80), intent(in) :: strahl_dir

integer :: i

!Start with main parameter file
open(nch_w2, file=TRIM(strahl_dir)//strahl_prof_in)
write(nch_w2, '(A)') &
    '          ******************** ', &
    '          **** from ASTRA **** ', &
    '          ******************** ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    '      ' // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ELECTRON DENSITY (cm**-3)'
if (ngmax < ngrid) then
    write(nch_w2, 107) neg(1)*1.E+13
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) neg(i)/neg(1)
    enddo
else
    write(nch_w2, 107) ne(1)*1.E+13
    do i=1, ngrid
        write(nch_w2, 101) ne(i)/ne(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ne decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', ne_decayl
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    "      " // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ELECTRON TEMPERATURE (eV)'
if (ngmax < ngrid) then
    write(nch_w2, 107) teg(1)*1.E+3
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) teg(i)/teg(1)
    enddo
else
    write(nch_w2, 107) te(1)*1.E+3
    do i=1, ngrid
        write(nch_w2, 101) te(i)/te(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv te decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', te_decayl
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv    time-vector  ', &
    '   1'
write(nch_w2, 103) '   ', time
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv      Function', &
    "      'interp'   ", &
    '   ', &
    '   ', &
    'cv    x-coordinate', &
    "      " // rho_coord, &
    '   ', &
    '   ', &
    'cv   # of radial points'
write(nch_w2, *) '       ', min(ngmax, ngrid)
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv   x-grid'
do i=1, min(ngmax, ngrid)
    if (ngmax < ngrid) then
        write(nch_w2, 101) rhopolg(i)
    else
        write(nch_w2, 101) rhopol(i)
    endif
enddo
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ION TEMPERATURE (eV)'
if (ngmax < ngrid) then
    write(nch_w2, 107) tig(1)*1.E+3
    do i=1, min(ngmax, ngrid)
        write(nch_w2, 101) tig(i)/tig(1)
    enddo
else
    write(nch_w2, 107) ti(1)*1.E+3
    do i=1, ngrid
        write(nch_w2, 101) ti(i)/ti(1)
    enddo
endif
write(nch_w2, '(A)') &
    '   ', &
    '   ', &
    'cv ti decay length [cm] in rho_volume'
write(nch_w2, 105) '   ', ti_decayl
close(nch_w2)

101 format(F12.4)
103 format(A, F12.4, A, F12.4, A, F12.4, A, F12.4)
105 format(A, F12.4)
107 format(E16.4)

return
end subroutine profiles_file_write_strahl

!----------------------------------------------------------------------
subroutine grid_write_strahl(strahl_dir, nfour_c, Raxis, &
  Rvoltot, Vloop, time, machine)

use parameters_a2spider, only: equil_now

implicit none

integer, parameter :: nch_r3=33, nch_w3=43, nr_max=100
character(len=17) :: strahl_grid_in='nete/grid_11111.1'

integer, intent(in) :: nfour_c
double precision, intent(in) :: Raxis, Rvoltot, Vloop, time
character(len=80), intent(in) :: strahl_dir
character(len=4) , intent(in) :: machine

integer :: i, j, nequil
double precision :: rpol, rvol, sigma, dum1, bdb0, bdb02, bmaxt, btor, fc

nequil = SIZE(equil_now%profiles_1d%volume)

open(nch_w3, file=TRIM(strahl_dir)//strahl_grid_in)
write(nch_w3, '(A)') &
    '   ', &
    'cv  rho volume(LCFS)[cm]  R_axis[cm]   U_loop[V]    time[s] '
write(nch_w3, 103) '  ', Rvoltot*100., '  ', Raxis*100., '   ', Vloop, '  ', time
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  number of grid points  points up to separtrix  fourier coefficients'
write(nch_w3, *) ' ', nequil, '  ', nequil, '  ', nfour_c
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  sqrt( (Psi-Psi_ax) / (Psi_sep - Psi_ax) )   '
do i=1, nequil
    rpol = ((equil_now%profiles_1d%psi(i)      - equil_now%profiles_1d%psi(1)) / &
            (equil_now%profiles_1d%psi(nequil) - equil_now%profiles_1d%psi(1)))**0.5
    write(nch_w3, 101) rpol
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv     rho volume / rho_volume(LCFS)    '

do i=1, nequil
    rvol = sqrt(equil_now%profiles_1d%volume(i)/equil_now%profiles_1d%volume(nequil))
    write(nch_w3, 101) rvol
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv   large radius low field side / R_axis '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%r_outboard(i)/equil_now%profiles_1d%r_outboard(1)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv   large radius high field side / R_axis '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%r_inboard(i)/equil_now%profiles_1d%r_inboard(1)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  safety factor '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%q(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  fraction of circulating particles  '

btor = equil_now%global_param%toroid_field%b0
do i=1, nequil
    BDB0  = equil_now%profiles_1d%bdb0(i)
    BMAXT = equil_now%profiles_1d%bmaxt(i)
    BDB02 = equil_now%profiles_1d%gm4(i)
    dum1  = min(.999999d0, BDB0/BMAXT)
    sigma = 1. - BDB02/BDB0**2 * (1. - SQRT(1. - dum1)*(1. + 0.5*dum1))
    dum1  = 1. - (BDB02/btor**2)*equil_now%profiles_1d%fofb(i)
    fc = 1. - 0.75*sigma - 0.25*dum1
    if (i == 1) fc = 1.
    write(nch_w3, 101) fc
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  Integral( dl_p / B_p) [m/T]  '
do i=1, nequil
    write(nch_w3, 101) 1./equil_now%profiles_1d%dPSIdV(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < B_total > [T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%bdb0(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < B_total**2 > [T**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm4(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < 1./B_total**2 > [1/T**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm5(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < R B_T > [m*T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%F_dia(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < R**2 B_p**2/B**2 > [m**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%rbp_b2(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  < 1/R**2 > [1/m**2]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%gm1(i)
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <cos (m theta) B_total**2> [T**2]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%acosB2a(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <sin (m theta) B_total**2> [T**2]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%asinB2a(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <cos (m theta) B_total ln(B_total)> [Tln(T)]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%acosBlnBa(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  <sin (m theta) B_total ln(B_total)> [Tln(T)]'
do j=1, nfour_c
    do i=1, nequil
        write(nch_w3, 101) equil_now%profiles_1d%asinBlnBa(i, j)
    enddo
enddo
write(nch_w3, '(A)') &
    '   ', &
    '   ', &
    'cv  Bp at LFS equator [T]  '
do i=1, nequil
    write(nch_w3, 101) equil_now%profiles_1d%bplfs(i)
enddo
close(nch_w3)

101 format(F12.4)
103 format(A, F12.4, A, F12.4, A, F12.4, A, F12.4)

return
end subroutine grid_write_strahl

end module strahl_mod
