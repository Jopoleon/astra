subroutine find_new_X0Y0(Nr, Nt, PSI, X, Y, X0, Y0, psiax, iax, jax)

use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Nr, Nt
double precision, intent(in), dimension(Nr, Nt) :: PSI, X, Y

integer, intent(out) :: iax, jax
double precision, intent(out) :: X0, Y0, psiax

integer :: jt, info, jmin(2)
double precision :: g3(2*Nt+1)
double precision :: work(2*(2*Nt+1)*6), matrix(2*Nt+1, 6)

call markloc_ef('find_new_X0Y0', debug_lev=debug)

jmin = minloc(PSI)
iax = jmin(1)
jax = jmin(2)

if (iax /= 1) then
    x0 = x(iax, jax)
    y0 = y(iax, jax)
    psiax = PSI(iax, jax)
    return
endif

g3(1)            = PSI(1, 1)
g3(2: Nt+1)      = PSI(2, 1: Nt)
g3(2+Nt: 2*Nt+1) = PSI(3, 1: Nt)
matrix(1, 1) = x(1, 1)**2.   
matrix(1, 2) = x(1, 1)   
matrix(1, 3) = 1.   
matrix(1, 4) = y(1, 1)**2.   
matrix(1, 5) = y(1, 1)   
matrix(1, 6) = x(1, 1)*y(1, 1)   
do jt = 1, Nt
    matrix(jt + 1, 1) = x(2, jt)**2.
    matrix(jt + 1, 2) = x(2, jt)   
    matrix(jt + 1, 3) = 1.   
    matrix(jt + 1, 4) = y(2, jt)**2.   
    matrix(jt + 1, 5) = y(2, jt)   
    matrix(jt + 1, 6) = x(2, jt)*y(2, jt)   
enddo 
do jt = 1, Nt
    matrix(jt + 1 + Nt, 1) = x(3, jt)**2.   
    matrix(jt + 1 + Nt, 2) = x(3, jt)   
    matrix(jt + 1 + Nt, 3) = 1.   
    matrix(jt + 1 + Nt, 4) = y(3, jt)**2.   
    matrix(jt + 1 + Nt, 5) = y(3, jt)   
    matrix(jt + 1 + Nt, 6) = x(3, jt)*y(3, jt)   
enddo 

! Lapack DGELS
call dgels('N', 2*Nt + 1, 6, 1, matrix, 2*Nt + 1, g3, 2*Nt + 1, &
    WORK, 2*(2*Nt + 1)*6, INFO )    

x0 = (g3(6)*g3(5) - 2*g3(4)*g3(2))/(4*g3(1)*g3(4) - g3(6)**2.)
y0 = (g3(6)*g3(2) - 2*g3(1)*g3(5))/(4*g3(1)*g3(4) - g3(6)**2.)
psiax = g3(1)*x0**2. + g3(2)*x0 + g3(3) +  &
        g3(4)*y0**2. + g3(5)*y0 + g3(6)*x0*y0

return
end subroutine find_new_X0Y0
