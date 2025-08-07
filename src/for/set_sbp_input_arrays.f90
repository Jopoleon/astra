subroutine set_sbp_input_arrays

use const_inc, only: NA1
use ipc_mod, only: n_sbp_arr_in  
use status_inc, only: aimpt, amain, ametr, cu, elon, er, fp, g11, ipol, mu, nalf, ndeut, ne, nhe3, nhydr, ni, nibm, nimpt, niz1, niz2, niz3, nn, ntrit, pblon, pbper, pfast, rho, shear, shif, te, ti, tn, tria, upl, vpol, vr, vrs, vtor, zef, zim1, zim2, zim3, zimpt, zmain

implicit none

double precision, allocatable, dimension(:, :) :: sbp_in

allocate(sbp_in(NA1, n_sbp_arr_in))
sbp_in = 0.d0

sbp_in(:,  1) = ne(1: NA1)
sbp_in(:,  2) = te(1: NA1)
sbp_in(:,  3) = ni(1: NA1)
sbp_in(:,  4) = ndeut(1: NA1)
sbp_in(:,  5) = ntrit(1: NA1)
sbp_in(:,  6) = niz1(1: NA1)
sbp_in(:,  7) = niz2(1: NA1)
sbp_in(:,  8) = ti(1: NA1)
sbp_in(:,  9) = zef(1: NA1)
sbp_in(:, 10) = zim1(1: NA1)
sbp_in(:, 11) = amain(1: NA1)
sbp_in(:, 12) = mu(1: NA1)
sbp_in(:, 13) = rho(1: NA1)
sbp_in(:, 14) = ametr(1: NA1)
sbp_in(:, 15) = shif(1: NA1)
sbp_in(:, 16) = elon(1: NA1)
sbp_in(:, 17) = tria(1: NA1)
sbp_in(:, 18) = er(1: NA1)
sbp_in(:, 19) = nibm(1: NA1)
sbp_in(:, 20) = g11(1: NA1)
sbp_in(:, 21) = vpol(1: NA1)
sbp_in(:, 22) = vrs(1: NA1)
sbp_in(:, 23) = vtor(1: NA1)
sbp_in(:, 24) = shear(1: NA1)
sbp_in(:, 25) = pblon(1: NA1)
sbp_in(:, 26) = pbper(1: NA1)
sbp_in(:, 27) = pfast(1: NA1)
sbp_in(:, 28) = niz3(1: NA1)
sbp_in(:, 29) = zim2(1: NA1)
sbp_in(:, 30) = zim3(1: NA1)
sbp_in(:, 31) = zimpt(1: NA1)
sbp_in(:, 32) = nimpt(1: NA1)
sbp_in(:, 33) = aimpt(1: NA1)

call fill_arr2shm(sbp_in)

return
end subroutine set_sbp_input_arrays
