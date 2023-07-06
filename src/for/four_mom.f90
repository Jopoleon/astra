subroutine four_mom(n_rho, n_mom, rcos, rsin, zcos, zsin)

use const_inc, only: GP
use parameters_a2equil, only: equil_now

implicit none

integer, intent(in) :: n_rho, n_mom
double precision, intent(out), dimension(n_rho, 0:n_mom) :: rcos, rsin, zcos, zsin

integer :: jrho, jmom, jthe, n_theta
double precision :: dtet, tet0, tet1, tetc, drc, dzc, sint0, sint1, cost0, cost1, &
    r0, r1, z0, z1, rc, zc, rcos1, rcos2, rsin1, rsin2, zcos1, zcos2, zsin1, zsin2

rcos = 0.
rsin = 0.
zcos = 0.
zsin = 0.
n_theta = SIZE(equil_now%coord_sys%position%teta2d)

do jrho=1, n_rho
    do jthe=2, n_theta-1
        rc = 0.5d0*(equil_now%coord_sys%position%r(jrho, jthe) + equil_now%coord_sys%position%r(jrho, jthe+1))              
        zc = 0.5d0*(equil_now%coord_sys%position%z(jrho, jthe) + equil_now%coord_sys%position%z(jrho, jthe+1))              
        dtet = 0.5d0*(equil_now%coord_sys%position%teta2d(jthe+1) - equil_now%coord_sys%position%teta2d(jthe))/GP
        rcos(jrho, 0) = rcos(jrho, 0) + rc*dtet
        zcos(jrho, 0) = zcos(jrho, 0) + zc*dtet
    enddo

    if (jrho > 1) then
        do jmom=1, n_mom
            do jthe=2, n_theta-1
                tet0 = equil_now%coord_sys%position%teta2d(jthe  )         
                tet1 = equil_now%coord_sys%position%teta2d(jthe+1)
                r0 = equil_now%coord_sys%position%r(jrho, jthe)
                r1 = equil_now%coord_sys%position%r(jrho, jthe+1)
                z0 = equil_now%coord_sys%position%z(jrho, jthe)
                z1 = equil_now%coord_sys%position%z(jrho, jthe+1)
                rc = 0.5d0*(r0 + r1)   
                zc = 0.5d0*(z0 + z1) 
                tetc = 0.5d0*(tet0 + tet1)*jmom         
                dtet = (tet1 - tet0)/GP

                drc = (r1 - r0)/(tet1 - tet0)
                dzc = (z1 - z0)/(tet1 - tet0)
                sint0 = dsin(tet0*jmom)
                sint1 = dsin(tet1*jmom)
                cost0 = dcos(tet0*jmom)
                cost1 = dcos(tet1*jmom)

                rcos1 =  rc*(sint1 - sint0)/jmom/GP
                zcos1 =  zc*(sint1 - sint0)/jmom/GP
                rsin1 = -rc*(cost1 - cost0)/jmom/GP
                zsin1 = -zc*(cost1 - cost0)/jmom/GP

                rcos2 = drc*( 0.5d0*(sint0 + sint1)*dtet + (cost1 - cost0)/GP/jmom )/jmom
                zcos2 = dzc*( 0.5d0*(sint0 + sint1)*dtet + (cost1 - cost0)/GP/jmom )/jmom
                rsin2 = drc*(-0.5d0*(cost0 + cost1)*dtet + (sint1 - sint0)/GP/jmom )/jmom
                zsin2 = dzc*(-0.5d0*(cost0 + cost1)*dtet + (sint1 - sint0)/GP/jmom )/jmom
 
                rcos(jrho, jmom) = rcos(jrho, jmom) + rcos1 + rcos2
                zcos(jrho, jmom) = zcos(jrho, jmom) + zcos1 + zcos2
                rsin(jrho, jmom) = rsin(jrho, jmom) + rsin1 + rsin2
                zsin(jrho, jmom) = zsin(jrho, jmom) + zsin1 + zsin2
            enddo
        enddo
     endif
enddo

return
end subroutine four_mom
