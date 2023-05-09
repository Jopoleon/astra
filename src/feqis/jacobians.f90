subroutine jacobians(nrho, ntheta, Rb, Zb, X0, Y0, lambda2d_in, lambda2dp_in, &
    psig, psiax, psib, relambda_flag, psigp, PSI, &
! Output
    dArea, X2, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, X, Y, &
    thetap, thetap_i, lambda2d, lambda2dp)

use pi_vars, only: GPI2
use numerical_tools, only: linterp

implicit none

logical, intent(in) :: relambda_flag
integer, intent(in) :: ntheta, nrho
double precision, intent(in) :: psiax, psib, X0, Y0
double precision, intent(in), dimension(ntheta) :: Rb, Zb
double precision, intent(in), dimension(nrho) :: psig, psigp
double precision, intent(in), dimension(nrho, ntheta) :: psi, lambda2d_in, lambda2dp_in

double precision, intent(out), dimension(nrho) :: ddr, ddr_i
double precision, intent(out), dimension(ntheta) :: dtp, dtm, dt_i
double precision, intent(out), dimension(ntheta+1) :: thetap, thetap_i
double precision, intent(out), dimension(nrho, ntheta) :: dArea, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    X, Y, X2, lambda2d, lambda2dp

integer :: jrho, jthe, jthe_l, jthe_r, k
double precision :: drdX, dhdX, drdY, dhdY, &
    lambda2di, lambda2dpi, Mdet_inv2, Mdet_invi1, dXb0_i, &
    dXdr2, dYdr2, dXdh2, dYdh2, &
    dXdri1, dYdri1, dXdhi1, dYdhi1, &
    dXdr, dYdr, dXdh, dYdh, &
    Jcbn, Jcbn2, grt2, gradr2, &
    dxcos1, dxcos2, dxsin1, dxsin2 
double precision, dimension(nrho) ::  lambda2dold, psin
double precision, dimension(ntheta) :: dXb0
double precision, dimension(nrho, ntheta) :: &
    Y2, X_i, Y_i, &
    gradhi1, grti1, &
    X_i1, Y_i1, Jcbni1

do jthe=1, ntheta
    dXb0(jthe) = sqrt((rb(jthe) - X0)**2 + (zb(jthe) - Y0)**2)
    thetap(jthe) = ATAN2(zb(jthe) - Y0, rb(jthe) - X0)
    if (thetap(jthe) < 0) thetap(jthe) = thetap(jthe) + GPI2
enddo
if (thetap(1) > thetap(ntheta)) then ! shift some entries
    do k=1, ntheta-1
        if (thetap(k+1) < thetap(k)) then
            thetap(1: k) = thetap(1: k) - GPI2
            EXIT
        endif
    enddo
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
ddr_i(1)  = 0.5*psigp(1)
ddr(nrho) = 0.

lambda2d  = lambda2d_in
lambda2dp = lambda2dp_in

if (relambda_flag) then !relambda
    do jthe=1, ntheta
        psin = (psi(:, jthe) - psiax)/(psib - psiax)
        psin(1) = 0.
        psin(nrho) = 1.
        lambda2dold = lambda2d(:, jthe)
        call linterp(psin, lambda2dold, nrho, &
            psig (2: nrho-1), lambda2d (2: nrho-1, jthe), nrho-2)
        call linterp(psig, lambda2d(:, jthe), nrho, &
            psigp(1: nrho-1), lambda2dp(1: nrho-1, jthe), nrho-1)
    enddo
endif

do jthe=1, ntheta
    if (jthe == ntheta) then
        jthe_r = 1
    else
        jthe_r = jthe + 1
    endif
    dXb0_i = 0.5*(dXb0(jthe) + dXb0(jthe_r))
    dxcos1 = dXb0(jthe)*cos(thetap(jthe))
    dxsin1 = dXb0(jthe)*sin(thetap(jthe))
    dxcos2 = dXb0_i*cos(thetap_i(jthe))
    dxsin2 = dXb0_i*sin(thetap_i(jthe))
    do jrho=1, nrho
        lambda2di  = 0.5*(lambda2d (jrho, jthe) + lambda2d (jrho, jthe_r))
        lambda2dpi = 0.5*(lambda2dp(jrho, jthe) + lambda2dp(jrho, jthe_r))
        X   (jrho, jthe) = X0 + lambda2d (jrho, jthe)*dxcos1
        Y   (jrho, jthe) = Y0 + lambda2d (jrho, jthe)*dxsin1
        X2  (jrho, jthe) = X0 + lambda2dp(jrho, jthe)*dxcos1
        Y2  (jrho, jthe) = Y0 + lambda2dp(jrho, jthe)*dxsin1
        X_i1(jrho, jthe) = X0 + lambda2di *dxcos2
        Y_i1(jrho, jthe) = Y0 + lambda2di *dxsin2
        X_i (jrho, jthe) = X0 + lambda2dpi*dxcos2
        Y_i (jrho, jthe) = Y0 + lambda2dpi*dxsin2
    enddo
enddo

! Compute jacobian
! J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r, theta)

do jthe=1, ntheta
    jthe_l = jthe - 1
    jthe_r = jthe + 1
    if (jthe == 1) then
        jthe_l = ntheta
    elseif (jthe == ntheta) then
        jthe_r = 1
    endif
    dtp (jthe) = thetap  (jthe  +1) - thetap  (jthe)
    dtm (jthe) = thetap  (jthe_l+1) - thetap  (jthe_l)
    dt_i(jthe) = thetap_i(jthe_l+1) - thetap_i(jthe_l)
    do jrho=1, nrho
        if (jrho == 1) then
            dXdr2 = (X(jrho+1, jthe) - X(jrho, jthe))/ddr(jrho)
            dYdr2 = (Y(jrho+1, jthe) - Y(jrho, jthe))/ddr(jrho)
            dXdh2 = (X_i(jrho, jthe) - X_i(jrho, jthe_l))/dt_i(jthe)
            dYdh2 = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l))/dt_i(jthe)
            dXdr   = (X2 (jrho, jthe) - X(jrho, jthe)) / psigp(jrho)
            dYdr   = (Y2 (jrho, jthe) - Y(jrho, jthe)) / psigp(jrho)
            dXdri1 = (X_i(jrho, jthe) - X(jrho, jthe)) / psigp(jrho)
            dYdri1 = (Y_i(jrho, jthe) - Y(jrho, jthe)) / psigp(jrho)
            dXdh   = (X_i(jrho, jthe) - X_i(jrho, jthe_l)) / dt_i(jthe)
            dYdh   = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l)) / dt_i(jthe)
            dXdhi1 = (X2(jrho, jthe_r) - X2(jrho, jthe)) / dtp(jthe)
            dYdhi1 = (Y2(jrho, jthe_r) - Y2(jrho, jthe)) / dtp(jthe)
        else if (jrho == nrho) then
            dXdr   = (X   (jrho, jthe) - X2 (jrho-1, jthe)) / ddr_i(jrho)
            dYdr   = (Y   (jrho, jthe) - Y2 (jrho-1, jthe)) / ddr_i(jrho)
            dXdri1 = (X_i1(jrho, jthe) - X_i(jrho-1, jthe)) / ddr_i(jrho)
            dYdri1 = (Y_i1(jrho, jthe) - Y_i(jrho-1, jthe)) / ddr_i(jrho)
            dXdh   = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l)) / dt_i(jthe)
            dYdh   = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l)) / dt_i(jthe)
            dXdhi1 = (X(jrho, jthe_r) - X(jrho, jthe)) / dtp(jthe)
            dYdhi1 = (Y(jrho, jthe_r) - Y(jrho, jthe)) / dtp(jthe)
        else
            dXdr2 = (X(jrho+1, jthe) - X(jrho, jthe))/ddr(jrho)
            dYdr2 = (Y(jrho+1, jthe) - Y(jrho, jthe))/ddr(jrho)
            dXdh2 = (X_i(jrho, jthe) - X_i(jrho, jthe_l))/dt_i(jthe)
            dYdh2 = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l))/dt_i(jthe)
            dXdr   = (X2  (jrho, jthe) - X2 (jrho-1, jthe)) / ddr_i(jrho)
            dYdr   = (Y2  (jrho, jthe) - Y2 (jrho-1, jthe)) / ddr_i(jrho)
            dXdri1 = (X_i (jrho, jthe) - X_i(jrho-1, jthe)) / ddr_i(jrho)
            dYdri1 = (Y_i (jrho, jthe) - Y_i(jrho-1, jthe)) / ddr_i(jrho)
            dXdh   = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l)) / dt_i(jthe)
            dYdh   = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l)) / dt_i(jthe)
            dXdhi1 = (X (jrho, jthe_r) - X(jrho, jthe)) / dtp(jthe)
            dYdhi1 = (Y (jrho, jthe_r) - Y(jrho, jthe)) / dtp(jthe)
        endif

        Mdet_inv2 = 1./(dXdr2*dYdh2 - dXdh2*dYdr2)
        drdX =  dYdh2*Mdet_inv2
        dhdX = -dYdr2*Mdet_inv2
        drdY = -dXdh2*Mdet_inv2
        dhdY =  dXdr2*Mdet_inv2
        gradr2 = drdX**2 + drdY**2
        grt2   = drdX*dhdX + drdY*dhdY   ! i+1/2, j

        Mdet_invi1 = 1./(dXdri1*dYdhi1 - dXdhi1*dYdri1)
        drdX =  dYdhi1*Mdet_invi1
        dhdX = -dYdri1*Mdet_invi1
        drdY = -dXdhi1*Mdet_invi1
        dhdY =  dXdri1*Mdet_invi1
        if ((jrho == 0) .or. (jrho == nrho)) then
            gradhi1(jrho, jthe) = 0.
            grti1  (jrho, jthe) = 0.
        else
            gradhi1(jrho, jthe) = dhdX**2 + dhdY**2
            grti1  (jrho, jthe) = drdX*dhdX + drdY*dhdY       ! i, j+1/2
        endif
        Jcbn  = dXdr *dYdh  - dXdh *dYdr    ! i, j
        Jcbn2 = dXdr2*dYdh2 - dXdh2*dYdr2   ! i+1/2, j
        Jcbni1(jrho, jthe) = dXdri1*dYdhi1 - dXdhi1*dYdri1  ! i, j+1/2
        dArea    (jrho, jthe) = Jcbn  *ddr_i(jrho)*dt_i(jthe)
        dArea2   (jrho, jthe) = Jcbn2 *ddr  (jrho)*dt_i(jthe)
        dArc_rp1 (jrho, jthe) = Jcbn2 *gradr2*dt_i(jthe)/X2(jrho, jthe)
        dArc_rpt1(jrho, jthe) = Jcbn2 *grt2  *dt_i(jthe)/X2(jrho, jthe)
        dArc_tp1 (jrho, jthe) = Jcbni1(jrho, jthe)*gradhi1(jrho, jthe)*ddr_i(jrho)/X_i1(jrho, jthe)
        dArc_tpr1(jrho, jthe) = Jcbni1(jrho, jthe)*grti1  (jrho, jthe)*ddr_i(jrho)/X_i1(jrho, jthe)
        if (jrho < nrho) then
            dArc_rm1 (jrho+1, jthe) = Jcbn2*gradr2*dt_i(jthe)/X2(jrho, jthe)
            dArc_rmt1(jrho+1, jthe) = Jcbn2*grt2  *dt_i(jthe)/X2(jrho, jthe)
        endif
    enddo
enddo

do jthe=1, ntheta
    jthe_l = jthe - 1
    jthe_r = jthe + 1
    if (jthe == 1) then
        jthe_l = ntheta
    elseif (jthe == ntheta) then
        jthe_r = 1
    endif
    do jrho=1, nrho
        dArc_tm1 (jrho, jthe) = Jcbni1(jrho, jthe_l)*gradhi1(jrho, jthe_l)*ddr_i(jrho)/X_i1(jrho, jthe_l)
        dArc_tmr1(jrho, jthe) = Jcbni1(jrho, jthe_l)*grti1  (jrho, jthe_l)*ddr_i(jrho)/X_i1(jrho, jthe_l)
    enddo
enddo
dArc_rm1 (1, :) = 0.   ! i-1/2, j
dArc_rmt1(1, :) = 0.
dArc_tm1 (1, :) = 0.   ! i, j-1/2
dArc_tmr1(1, :) = 0.   ! i, j-1/2

return
end subroutine jacobians
