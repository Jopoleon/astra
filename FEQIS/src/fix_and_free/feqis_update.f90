subroutine circuit_eq_advance_ef(j_init)
!circuit eq advance

use pi_grec_vars, only: gpi2
use ef_circuit, only: tau_new, tau_old, nconduc, cur_con_old, &
    curconduc, dpc, psiplasmatoconduc, psi_cur_old, &
    indconduc, resconduc, voltage
use exchange_with_astra, only: tau_circuit_ef, tau_gseq_ef, &
    activate_coil_ef, reconnect_circuits, force_coil, i_dim1, &
    n_equivalence, new_equivalence, current_limit_ef, use_reduce_circuit
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: j_init

integer :: i, ic, j, k, iii, jjj, invertcommand, i_equivalence, &
    i_cnew, firstcall
integer, dimension(200) :: rem_coils
double precision, dimension(200, 200) :: restemp, indtemp
double precision, dimension(i_dim1) :: dpctemp, vtemp, curotemp, curtemp

data firstcall/0/
save restemp, indtemp, rem_coils, i_cnew, firstcall

call markloc_ef('circuit_eq_advance_ef', debug_lev=debug)

tau_new = tau_circuit_ef
invertcommand = 0

if (tau_new /= tau_old) then
    invertcommand = 1
    tau_old = tau_new
endif
if (j_init == 0 .or. firstcall == 0) then
    invertcommand = 1
    tau_old = tau_new
    rem_coils = 0
endif

ic = nconduc

if (j_init == 0) cur_con_old = curconduc
if (j_init == 0) i_cnew = ic
if (j_init == 0) return

firstcall = 1

dpc(1: ic) = GPI2*(psiplasmatoconduc(1: ic) - psi_cur_old(1: ic))/tau_gseq_ef !plasma contribution
if (j_init == 0) dpc = 0.

do j = 1, nconduc
    if (activate_coil_ef(j) == 0) cur_con_old(j) = 0.
    if (activate_coil_ef(j) == 0) dpc(j) = 0.
enddo

if (reconnect_circuits == 1) then

    invertcommand = 1
    i_cnew = ic
    indtemp(1: ic, 1: ic) = indconduc(1: ic, 1: ic)
    restemp(1: ic, 1: ic) = resconduc(1: ic, 1: ic)

    do k=1, n_equivalence
        i_equivalence = 1
        do i=1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do j=i + 1, nconduc
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        indtemp(i, i) = indconduc(i, i) + indconduc(j, j) + indconduc(i, j) + indconduc(j, i)
                        restemp(i, i) = resconduc(i, i) + resconduc(j, j) + resconduc(i, j) + resconduc(j, i)
                        do iii=1, nconduc
                            if (iii /= i .and. iii /= j) then
                                indtemp(i, iii) = indconduc(i, iii) + indconduc(j, iii)
                                restemp(i, iii) = resconduc(i, iii) + resconduc(j, iii)
                                indtemp(iii, i) = indconduc(iii, i) + indconduc(iii, j)
                                restemp(iii, i) = resconduc(iii, i) + resconduc(iii, j)
                            endif
                        enddo
                        rem_coils(j) = new_equivalence(j, k)
                        ic = ic - 1
                    endif            
                enddo    
! Remove obsolete columns and rows
                do jjj=1, nconduc
                    j = nconduc - jjj + 1
                    if (rem_coils(j) > 0) then
                        indtemp(1: nconduc, j: nconduc) = indtemp(1: nconduc, j + 1: nconduc + 1)
                        restemp(1: nconduc, j: nconduc) = restemp(1: nconduc, j + 1: nconduc + 1)
                        indtemp(j: nconduc, 1: nconduc) = indtemp(j + 1: nconduc + 1, 1: nconduc)
                        restemp(j: nconduc, 1: nconduc) = restemp(j + 1: nconduc + 1, 1: nconduc)
                    endif
                enddo
            endif
        enddo    
    enddo

    i_cnew = ic
endif

if (use_reduce_circuit == 1) then
    ic = nconduc
    dpctemp(1: ic) = dpc(1: ic)
    vtemp(1: ic) = voltage(1: ic)
    curotemp(1: ic) = cur_con_old(1: ic)
    do k = 1, n_equivalence
        i_equivalence = 1
        do i = 1, nconduc
            if (new_equivalence(i, k) > 0 .and. i_equivalence == 1) then
                i_equivalence = 0
                do jjj = i + 1, nconduc
                    j = nconduc - jjj + i + 1
                    if (new_equivalence(j, k) == new_equivalence(i, k)) then
                        dpctemp(i) = dpctemp(i) + dpctemp(j)
                    ic = ic - 1
                    endif            
                enddo    
            endif
        enddo    
    enddo
    do jjj = 1, nconduc
        j = nconduc - jjj + 1
        if (rem_coils(j) > 0) then
            dpctemp(j: nconduc) = dpctemp(j + 1: nconduc + 1)    
            vtemp(j: nconduc) = vtemp(j + 1: nconduc + 1)    
            curotemp(j: nconduc) = curotemp(j + 1: nconduc + 1)
        endif
    enddo
endif

i = i_cnew

write(*, *) 'mmm', use_reduce_circuit, i, indconduc(5, 5),  & 
    resconduc(5, 5), cur_con_old(5), curconduc(5), voltage(5), &
    dpc(5), tau_new, invertcommand

if (use_reduce_circuit == 0) then
    call solve_circuit_equations(i, & 
        indconduc(1: i, 1: i), resconduc(1: i, 1: i),  & 
        cur_con_old(1: i), curconduc(1: i), & 
        voltage(1: i), dpc(1: i), tau_new, invertcommand)
else
    call solve_circuit_equations(i, & 
        indtemp(1: i, 1: i), restemp(1: i, 1: i),  & 
        curotemp(1: i), curtemp(1: i), & 
        vtemp(1: i), dpctemp(1: i), tau_new, invertcommand)    
! Readapt currents
    do i = 1, nconduc
        if (rem_coils(i) == 0) then
            curconduc(i) = curtemp(i)
        endif
        if (rem_coils(i) > 0) then
            curconduc(i) = curtemp(rem_coils(i))
            curtemp(i + 1: i_cnew + 1) = curtemp(i: i_cnew)
        endif
    enddo
endif

write(*, *) 'cur ', cur_con_old(5), curconduc(5), dpc(5), tau_new

ic = nconduc
i = ic
do j = 1, ic
    if (activate_coil_ef(j) == 0) curconduc(j) = 0.
    curconduc(j) = max(curconduc(j), current_limit_ef(j, 2))
    curconduc(j) = min(curconduc(j), current_limit_ef(j, 1))
    if (sum(force_coil(j, 1: i)) > 0.5) then
        do k = 1, i
            if (force_coil(j, k) == 1.) then
                curconduc(j) = curconduc(k)
            endif
        enddo
    endif
enddo

reconnect_circuits = 0

return
end subroutine circuit_eq_advance_ef

!---------------------------------------------------------
subroutine solve_circuit_equations(nc, im, rm, I0, I1, & 
    V, dpc, tau, invertcommand)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: nc, invertcommand
double precision, intent(in) :: tau
double precision, intent(in) , dimension(nc) :: i0, v, dpc
double precision, intent(in) , dimension(nc, nc) :: im, rm
double precision, intent(out), dimension(nc) :: i1

integer :: i
double precision, dimension(nc) :: b(nc)
double precision, dimension(nc, nc) :: matrix 
double precision, dimension(200, 200) :: invmatrix

save invmatrix

call markloc_ef('solve_circuit_equations', debug_lev=debug)

!equation is im*(i1 - i0)/tau  +  rm*i1 = v - dpc

do i = 1, nc
    b(i) = v(i) - dpc(i) + sum(im(i, 1: nc)*i0(1: nc))/tau
enddo

if (invertcommand == 1) then
    matrix(1: nc, 1: nc) = im(1: nc, 1: nc)/tau + rm(1: nc, 1: nc)
    call inverse_matrix_feqis(matrix, invmatrix(1: nc, 1: nc), nc)
endif

write(*, *) 'ii', i0(12: 15)
do i = 1, nc
    i1(i) = sum(invmatrix(i, 1: nc)*b(1: nc))
enddo
write(*, *) 'ii', i1(12: 15)

return
end subroutine solve_circuit_equations
