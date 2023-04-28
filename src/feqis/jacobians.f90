subroutine jacobians(nrho, ntheta, Rb, Zb, X0, Y0, lambda2d_in, lambda2dp_in, &
    psig, psiax, psib, j_ok, psigp, PSI, &
! Output
    dArea, X2, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, X, Y, &
    thetap, thetap_i, lambda2d_out, lambda2dp_out)

use interp_mod, only: linterp

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI

integer, intent(in) :: ntheta, nrho, j_ok
double precision, intent(in) :: psiax, psib, X0, Y0
double precision, intent(in), dimension(ntheta) :: Rb, Zb
double precision, intent(in), dimension(nrho) :: psig, psigp
double precision, intent(in), dimension(nrho, ntheta) :: psi, lambda2d_in, lambda2dp_in

double precision, intent(out), dimension(nrho) :: ddr, ddr_i
double precision, intent(out), dimension(ntheta) :: dtp, dtm, dt_i
double precision, intent(out), dimension(ntheta+1) :: thetap, thetap_i
double precision, intent(out), dimension(nrho, ntheta) :: dArea, dArea2, X2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    X, Y, lambda2d_out, lambda2dp_out

integer :: jrho, jthe, jthe_l, jthe_r, j, k
double precision :: drdX, dhdX, drdY, dhdY, Mdet_inv, dpsi, dthe
double precision, dimension(nrho) ::  lambda2dold
double precision, dimension(ntheta) :: dXb0, dXb0_i
double precision, dimension(nrho, ntheta) :: lambda2d, lambda2dp, &
    lambda2di, lambda2dpi, psin, Y2, X_i, Y_i, &
    dXdr2, dYdr2, dXdr, dYdr, dXdri1, dYdri1, &
    dXdh2, dYdh2, dXdh, dYdh, dXdhi1, dYdhi1, &
    Jcbni1, gradhi1, grti1, &
    rr2, r_i, r_i1, X_i1, Y_i1, grt2, rmin, &
    Jcbn, Jcbn2, gradr2

do jthe=1, ntheta
    dXb0(jthe) = sqrt((rb(jthe) - X0)**2 + (zb(jthe) - Y0)**2)
    thetap(jthe) = ATAN2(zb(jthe) - Y0, rb(jthe) - X0)
    if (thetap(jthe) < 0) thetap(jthe) = thetap(jthe) + GPI2
enddo

if (thetap(1) > thetap(ntheta)) then !reorder
    j = 1
    do k=1, ntheta-1
        if (thetap(k+1) < thetap(k)) j = k + 1   ! j is the first teta above 0
    enddo
    if (j > 1) thetap(1: j-1) = thetap(1: j-1) - GPI2
endif
thetap(ntheta+1) = thetap(1) + GPI2

! Now find the intermediate theta grid
do jthe=1, ntheta
    thetap_i(jthe) = 0.5*(thetap(jthe) + thetap(jthe+1))
enddo
thetap_i(ntheta+1) = thetap_i(1) + GPI2

! Compute differentials

do jrho=1, nrho-1
    ddr(jrho)     = psig (jrho+1) - psig (jrho)
    ddr_i(jrho+1) = psigp(jrho+1) - psigp(jrho)
enddo
ddr_i(1) = 0.5*psigp(1)
ddr(nrho)  = 0.

do jthe=1, ntheta
    jthe_l = jthe - 1
    if (jthe == 1) then
        jthe_l = ntheta
    endif
    dtp (jthe) = thetap  (jthe  +1) - thetap  (jthe)
    dtm (jthe) = thetap  (jthe_l+1) - thetap  (jthe_l)
    dt_i(jthe) = thetap_i(jthe_l+1) - thetap_i(jthe_l)
enddo

do jthe=1, ntheta-1
    dXb0_i(jthe) = 0.5*(dXb0(jthe) + dXb0(jthe+1))
enddo
dXb0_i(ntheta) = dXb0_i(1)

lambda2d = lambda2d_in
lambda2dp = lambda2dp_in

if (j_ok == 1) then !relambda
    psin = (psi - psiax)/(psib - psiax)
    psin(   1, :) = 0.
    psin(nrho, :) = 1.
    do jthe=1, ntheta
        lambda2dold = lambda2d(:, jthe)
        call linterp(psin(:, jthe), lambda2dold, nrho, &
            psig (2: nrho-1), lambda2d (2: nrho-1, jthe), nrho-2)
        call linterp(psig, lambda2d(:, jthe), nrho, &
            psigp(1: nrho-1), lambda2dp(1: nrho-1, jthe), nrho-1)
    enddo
endif

do jthe=1, ntheta-1
    lambda2di (:, jthe) = 0.5*(lambda2d (:, jthe) + lambda2d (:, jthe+1))
    lambda2dpi(:, jthe) = 0.5*(lambda2dp(:, jthe) + lambda2dp(:, jthe+1))
enddo
lambda2di (:, ntheta) = 0.5*(lambda2d (:, ntheta) + lambda2d (:, 1))
lambda2dpi(:, ntheta) = 0.5*(lambda2dp(:, ntheta) + lambda2dp(:, 1))

do jthe=1, ntheta
    rmin(:, jthe) = lambda2d (:, jthe)*dXb0(jthe)
    rr2 (:, jthe) = lambda2dp(:, jthe)*dXb0(jthe)
    X   (:, jthe) = X0 + rmin(:, jthe)*cos(thetap(jthe))
    Y   (:, jthe) = Y0 + rmin(:, jthe)*sin(thetap(jthe))
    X2  (:, jthe) = X0 + rr2 (:, jthe)*cos(thetap(jthe))
    Y2  (:, jthe) = Y0 + rr2 (:, jthe)*sin(thetap(jthe))

    r_i1(:, jthe) = lambda2di (:, jthe)*dXb0_i(jthe)
    r_i (:, jthe) = lambda2dpi(:, jthe)*dXb0_i(jthe)
    X_i1(:, jthe) = X0 + r_i1(:, jthe)*cos(thetap_i(jthe))
    Y_i1(:, jthe) = Y0 + r_i1(:, jthe)*sin(thetap_i(jthe))
    X_i (:, jthe) = X0 + r_i (:, jthe)*cos(thetap_i(jthe))
    Y_i (:, jthe) = Y0 + r_i (:, jthe)*sin(thetap_i(jthe))
enddo

! Compute jacobian
! J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r, theta)

do jrho=1, nrho-1
    do jthe=1, ntheta
        if (jthe == 1) then
            jthe_l = ntheta
        else
            jthe_l = jthe - 1
        endif
        dXdr2(jrho, jthe) = (X(jrho+1, jthe) - X(jrho, jthe))/ddr(jrho)
        dYdr2(jrho, jthe) = (Y(jrho+1, jthe) - Y(jrho, jthe))/ddr(jrho)
        dXdh2(jrho, jthe) = (X_i(jrho, jthe) - X_i(jrho, jthe_l))/dt_i(jthe)
        dYdh2(jrho, jthe) = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l))/dt_i(jthe)
    enddo
enddo

Jcbn2 = dXdr2*dYdh2 - dXdh2*dYdr2   ! i+1/2, j

do jrho=2, nrho-1
    do jthe=1, ntheta
        jthe_l = jthe - 1
        jthe_r = jthe + 1
        if (jthe == 1) then
            jthe_l = ntheta
        elseif (jthe == ntheta) then
            jthe_r = 1
        endif
        dXdr  (jrho, jthe) = (X2  (jrho, jthe) - X2 (jrho-1, jthe)) / ddr_i(jrho)
        dYdr  (jrho, jthe) = (Y2  (jrho, jthe) - Y2 (jrho-1, jthe)) / ddr_i(jrho)
        dXdri1(jrho, jthe) = (X_i (jrho, jthe) - X_i(jrho-1, jthe)) / ddr_i(jrho)
        dYdri1(jrho, jthe) = (Y_i (jrho, jthe) - Y_i(jrho-1, jthe)) / ddr_i(jrho)
        dXdh  (jrho, jthe) = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l)) / dt_i(jthe)
        dYdh  (jrho, jthe) = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l)) / dt_i(jthe)
        dXdhi1(jrho, jthe) = (X (jrho, jthe_r) - X(jrho, jthe)) / dtp(jthe)
        dYdhi1(jrho, jthe) = (Y (jrho, jthe_r) - Y(jrho, jthe)) / dtp(jthe)
    enddo
enddo

jrho = 1
do jthe=1, ntheta
    jthe_l = jthe - 1
    jthe_r = jthe + 1
    if (jthe == 1) then
        jthe_l = ntheta
    elseif (jthe == ntheta) then
        jthe_r = 1
    endif
    dXdr  (jrho, jthe) = (X2 (jrho, jthe) - X(jrho, jthe)) / psigp(jrho)
    dYdr  (jrho, jthe) = (Y2 (jrho, jthe) - Y(jrho, jthe)) / psigp(jrho)
    dXdri1(jrho, jthe) = (X_i(jrho, jthe) - X(jrho, jthe)) / psigp(jrho)
    dYdri1(jrho, jthe) = (Y_i(jrho, jthe) - Y(jrho, jthe)) / psigp(jrho)
    dXdh  (jrho, jthe) = (X_i(jrho, jthe) - X_i(jrho, jthe_l)) / dt_i(jthe)
    dYdh  (jrho, jthe) = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l)) / dt_i(jthe)
    dXdhi1(jrho, jthe) = (X2(jrho, jthe_r) - X2(jrho, jthe)) / dtp(jthe)
    dYdhi1(jrho, jthe) = (Y2(jrho, jthe_r) - Y2(jrho, jthe)) / dtp(jthe)
enddo

jrho = nrho
do jthe=1, ntheta
    jthe_l = jthe - 1
    jthe_r = jthe + 1
    if (jthe == 1) then
        jthe_l = ntheta
    elseif (jthe == ntheta) then
        jthe_r = 1
    endif
    dXdr  (jrho, jthe) = (X   (jrho, jthe) - X2 (jrho-1, jthe)) / ddr_i(jrho)
    dYdr  (jrho, jthe) = (Y   (jrho, jthe) - Y2 (jrho-1, jthe)) / ddr_i(jrho)
    dXdri1(jrho, jthe) = (X_i1(jrho, jthe) - X_i(jrho-1, jthe)) / ddr_i(jrho)
    dYdri1(jrho, jthe) = (Y_i1(jrho, jthe) - Y_i(jrho-1, jthe)) / ddr_i(jrho)
    dXdh  (jrho, jthe) = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l)) / dt_i(jthe)
    dYdh  (jrho, jthe) = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l)) / dt_i(jthe)
    dXdhi1(jrho, jthe) = (X(jrho, jthe_r) - X(jrho, jthe)) / dtp(jthe)
    dYdhi1(jrho, jthe) = (Y(jrho, jthe_r) - Y(jrho, jthe)) / dtp(jthe)
enddo

Jcbn   = dXdr  *dYdh   - dXdh  *dYdr    ! i, j
Jcbni1 = dXdri1*dYdhi1 - dXdhi1*dYdri1  ! i, j+1/2

do jthe=1, ntheta
    do jrho=1, nrho-1
        Mdet_inv = 1./(dXdr2(jrho, jthe)*dYdh2(jrho, jthe) - dXdh2(jrho, jthe)*dYdr2(jrho, jthe))
        drdX =  dYdh2(jrho, jthe)*Mdet_inv
        dhdX = -dYdr2(jrho, jthe)*Mdet_inv
        drdY = -dXdh2(jrho, jthe)*Mdet_inv
        dhdY =  dXdr2(jrho, jthe)*Mdet_inv
        gradr2(jrho, jthe) = drdX**2 + drdY**2
        grt2  (jrho, jthe) = drdX*dhdX + drdY*dhdY   ! i+1/2, j

        Mdet_inv = 1./(dXdri1(jrho, jthe)*dYdhi1(jrho, jthe) - dXdhi1(jrho, jthe)*dYdri1(jrho, jthe))
        drdX =  dYdhi1(jrho, jthe)*Mdet_inv
        dhdX = -dYdri1(jrho, jthe)*Mdet_inv
        drdY = -dXdhi1(jrho, jthe)*Mdet_inv
        dhdY =  dXdri1(jrho, jthe)*Mdet_inv
        gradhi1(jrho, jthe) = dhdX**2 + dhdY**2
        grti1  (jrho, jthe) = drdX*dhdX + drdY*dhdY       ! i, j+1/2
    enddo
enddo

gradhi1(1, :) = 0.
grti1  (1, :) = 0.

! Compute differential area, arc lengths

do jrho=1, nrho
    dArea    (jrho, :) = Jcbn  (jrho, :)*ddr_i(jrho)*dt_i(:)
    dArea2   (jrho, :) = Jcbn2 (jrho, :)*ddr  (jrho)*dt_i(:)
    dArc_rp1 (jrho, :) = Jcbn2 (jrho, :)*gradr2 (jrho, :)*dt_i(:)/X2(jrho, :)
    dArc_rpt1(jrho, :) = Jcbn2 (jrho, :)*grt2   (jrho, :)*dt_i(:)/X2(jrho, :)
    dArc_tp1 (jrho, :) = Jcbni1(jrho, :)*gradhi1(jrho, :)*ddr_i(jrho)/X_i1(jrho, :)
    dArc_tpr1(jrho, :) = Jcbni1(jrho, :)*grti1  (jrho, :)*ddr_i(jrho)/X_i1(jrho, :)
enddo

dArc_rm1 (1, :) = 0.   ! i-1/2, j
dArc_rmt1(1, :) = 0.
dArc_tm1 (1, :) = 0.   ! i, j-1/2
dArc_tmr1(1, :) = 0.   ! i, j-1/2
do jrho=2, nrho
    dArc_rm1 (jrho, :) = Jcbn2(jrho-1, :)*gradr2(jrho-1, :)*dt_i(:)/X2(jrho-1, :)
    dArc_rmt1(jrho, :) = Jcbn2(jrho-1, :)*grt2  (jrho-1, :)*dt_i(:)/X2(jrho-1, :)
    do jthe=1, ntheta
        if (jthe == 1) then
            jthe_l = ntheta
        else
            jthe_l = jthe - 1
        endif
        dArc_tm1 (jrho, jthe) = Jcbni1(jrho, jthe_l)*gradhi1(jrho, jthe_l)*ddr_i(jrho)/X_i1(jrho, jthe_l)
        dArc_tmr1(jrho, jthe) = Jcbni1(jrho, jthe_l)*grti1  (jrho, jthe_l)*ddr_i(jrho)/X_i1(jrho, jthe_l)
    enddo
enddo

lambda2d_out  = lambda2d
lambda2dp_out = lambda2dp

return
end subroutine jacobians
