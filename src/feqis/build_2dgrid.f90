subroutine build_2dgrid(Nr, Nt, Rb, Zb, X0, Y0, lambda2d, lambda2dp, &
    psig, psiax, psib, j_ok, psigp, PSI, &
    dArea, Rmaj, Rmaj2, Rmaji, Rmaji1, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, Area, X, Y, X_i1, Y_i1, &
    r, thetap, rr2, thetap_i, r_i, r_i1, dphdR, dphdZ, gradr, Jcbn, &
    grt, gradr2, Jcbn2, grt2, ghht2)

use pi_vars, only: GPI2
use debugger_ef, only: markloc_ef, debug

implicit none

integer, intent(in) :: Nt, Nr, j_ok
double precision, intent(in) :: psiax, psib, X0, Y0
double precision, intent(in), dimension(Nt) :: Rb, Zb
double precision, intent(in), dimension(Nr) :: psig, psigp
double precision, intent(in), dimension(Nr, Nt) :: psi

double precision, intent(out) :: Area
double precision, intent(out), dimension(Nt+1) :: thetap, thetap_i
double precision, intent(out), dimension(Nr, Nt) :: dArea, dArea2, &
    Rmaj, Rmaj2, Rmaji, Rmaji1, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, X, Y, X_i1, Y_i1, &
    r, rr2, r_i, r_i1, dphdR, dphdZ, Jcbn, Jcbn2, &
    gradr, gradr2, grt, grt2, ghht2

integer :: jr, jt, jthe_l, jthe_r, j, k
double precision :: drdX, dhdX, drdY, dhdY, Mdet, dpsi, dthe
double precision, dimension(Nt) :: dXb0, dXb0_i
double precision, dimension(Nr, Nt) :: lambda2d, lambda2dold, &
    lambda2dp, lambda2di, lambda2dpi, psin, X2, Y2, X_i, Y_i, &
    dXdr2, dYdr2, dXdr, dYdr, dXdri1, dYdri1, &
    dXdh2, dYdh2, dXdh, dYdh, dXdhi1, dYdhi1, &
    Jcbni1, gradh, gradh2, gradhi1, grti1

call markloc_ef('build_2dgrid', debug_lev=debug)

do jt=1, Nt
    dXb0(jt) = sqrt((rb(jt) - X0)**2.0 + (zb(jt) - Y0)**2.0)
    thetap(jt) = ATAN2(zb(jt) - Y0, rb(jt) - X0)
    if (thetap(jt) < 0) thetap(jt) = thetap(jt) + GPI2
enddo

if (thetap(1) > thetap(Nt)) then !reorder
    j = 1
    do k=1, Nt-1
        if (thetap(k+1) < thetap(k)) j = k + 1   ! j is the first teta above 0
    enddo  
    if (j > 1) thetap(1: j-1) = thetap(1: j-1) - GPI2
endif 

jt = Nt
thetap(jt+1) = GPI2 + thetap(1)

! Now find the intermediate theta grid
do jt=1, Nt
    thetap_i(jt) = (thetap(jt) + thetap(jt+1))/2.0
enddo
thetap_i(Nt+1) = thetap_i(1) + GPI2
 
do jt=1, Nt-1
    dXb0_i(jt) = (dXb0(jt) + dXb0(jt+1))/2.0
enddo
dXb0_i(Nt) = dXb0_i(1)

if (j_ok == 1) then !relambda
!relambda
    lambda2dold(1: Nr, 1: Nt) = lambda2d(1: Nr, 1: Nt)
    do jt=1, Nt
        do jr=1, Nr
            psin(jr, jt) = (psi(jr, jt) - psiax)/(psib - psiax)
        enddo
        psin( 1, jt) = 0
        psin(Nr, jt) = 1.

        call linterp_feqis(psin(:, jt), lambda2dold(:, jt), Nr, &
            psig( 2: Nr-1), lambda2d( 2: Nr-1, jt), Nr-2)
        call linterp_feqis(psin(:, jt), lambda2dold(:, jt), Nr, &
            psigp(1: Nr-1), lambda2dp(1: Nr-1, jt), Nr-1)
        do jr=2, Nr-1
            lambda2d( jr, jt) = min(1., lambda2d( jr, jt))
            lambda2dp(jr, jt) = min(1., lambda2dp(jr, jt))
        enddo
    enddo
endif 

do jt=1, Nt
    do jr=1, Nr
        r(  jr, jt) = lambda2d( jr, jt)*dXb0(jt)
        rr2(jr, jt) = lambda2dp(jr, jt)*dXb0(jt)
        X( jr, jt) = X0 + r(  jr, jt)*cos(thetap(jt))
        Y( jr, jt) = Y0 + r(  jr, jt)*sin(thetap(jt))
        X2(jr, jt) = X0 + rr2(jr, jt)*cos(thetap(jt))
        Y2(jr, jt) = Y0 + rr2(jr, jt)*sin(thetap(jt))
    enddo
enddo

do jt=1, Nt-1
    lambda2di( :, jt) = (lambda2d( :, jt) + lambda2d( :, jt+1))/2.0
    lambda2dpi(:, jt) = (lambda2dp(:, jt) + lambda2dp(:, jt+1))/2.0
enddo
jt = Nt
lambda2di( :, jt) = (lambda2d( :, jt) + lambda2d( :, 1))/2.0
lambda2dpi(:, jt) = (lambda2dp(:, jt) + lambda2dp(:, 1))/2.0

do jt=1, Nt
    do jr = 1, Nr
        r_i1(jr, jt) = lambda2di( jr, jt)*dXb0_i(jt)
        r_i( jr, jt) = lambda2dpi(jr, jt)*dXb0_i(jt)
        X_i1(jr, jt) = X0 + r_i1(jr, jt)*cos(thetap_i(jt))
        Y_i1(jr, jt) = Y0 + r_i1(jr, jt)*sin(thetap_i(jt))
        X_i( jr, jt) = X0 + r_i( jr, jt)*cos(thetap_i(jt))
        Y_i( jr, jt) = Y0 + r_i( jr, jt)*sin(thetap_i(jt))
    enddo
enddo
   
! Now computes jacobian
! J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r, theta)

do jr=1, Nr-1
    dpsi = psig(jr+1) - psig(jr)
    do jt=1, Nt
        if (jt == 1) then
            jthe_l = Nt
        else
            jthe_l = jt - 1
        endif
        dthe = thetap_i(jthe_l+1) - thetap_i(jthe_l) 
        dXdr2(jr, jt) = (X(jr+1, jt) - X(jr, jt))/dpsi 
        dYdr2(jr, jt) = (Y(jr+1, jt) - Y(jr, jt))/dpsi
        dXdh2(jr, jt) = (X_i(jr, jt) - X_i(jr, jthe_l))/dthe
        dYdh2(jr, jt) = (Y_i(jr, jt) - Y_i(jr, jthe_l))/dthe
    enddo
enddo

Jcbn2 = dXdr2*dYdh2 - dXdh2*dYdr2   ! i+1/2, j 

do jr=2, Nr-1
    dpsi = psigp(jr) - psigp(jr-1)
    do jt=1, Nt
        jthe_l = jt - 1
        jthe_r = jt + 1
        if (jt == 1) then
            jthe_l = Nt
        elseif (jt == Nt) then
            jthe_r = 1
        endif 
        dXdr(  jr, jt) = (X2(  jr, jt) - X2(  jr-1, jt))/dpsi
        dYdr(  jr, jt) = (Y2(  jr, jt) - Y2(  jr-1, jt))/dpsi
        dXdri1(jr, jt) = (X_i( jr, jt) - X_i( jr-1, jt))/dpsi
        dYdri1(jr, jt) = (Y_i( jr, jt) - Y_i( jr-1, jt))/dpsi
        dXdh(  jr, jt) = (X_i1(jr, jt) - X_i1(jr, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
        dYdh(  jr, jt) = (Y_i1(jr, jt) - Y_i1(jr, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
        dXdhi1(jr, jt) = (X(jr, jthe_r) - X(jr, jt))/(thetap(jt+1) - thetap(jt))  
        dYdhi1(jr, jt) = (Y(jr, jthe_r) - Y(jr, jt))/(thetap(jt+1) - thetap(jt))
    enddo
enddo

jr = 1
do jt=1, Nt
    jthe_l = jt - 1
    jthe_r = jt + 1
    if (jt == 1) then
        jthe_l = Nt
    elseif (jt == Nt) then
        jthe_r = 1
    endif 
    dXdr(  jr, jt) = ( X2(jr, jt) - X(jr, jt))/ psigp(jr)
    dYdr(  jr, jt) = ( Y2(jr, jt) - Y(jr, jt))/ psigp(jr)
    dXdri1(jr, jt) = (X_i(jr, jt) - X(jr, jt))/ psigp(jr)
    dYdri1(jr, jt) = (Y_i(jr, jt) - Y(jr, jt))/ psigp(jr)
    dXdh(  jr, jt) = (X_i(jr, jt) - X_i(jr, jthe_l))/ (thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dYdh(  jr, jt) = (Y_i(jr, jt) - Y_i(jr, jthe_l))/ (thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dXdhi1(jr, jt) = (X2(jr, jthe_r) - X2(jr, jt))/ (thetap(jt+1) - thetap(jt))  
    dYdhi1(jr, jt) = (Y2(jr, jthe_r) - Y2(jr, jt))/ (thetap(jt+1) - thetap(jt))
enddo

jr = Nr
dpsi = psigp(jr) - psigp(jr-1)
do jt=1, Nt
    jthe_l = jt - 1
    jthe_r = jt + 1
    if (jt == 1) then
        jthe_l = Nt
    elseif (jt == Nt) then
        jthe_r = 1
    endif 
    dXdr(  jr, jt) = (X(   jr, jt) - X2( jr-1, jt))/dpsi
    dYdr(  jr, jt) = (Y(   jr, jt) - Y2( jr-1, jt))/dpsi
    dXdri1(jr, jt) = (X_i1(jr, jt) - X_i(jr-1, jt))/dpsi
    dYdri1(jr, jt) = (Y_i1(jr, jt) - Y_i(jr-1, jt))/dpsi
    dXdh(  jr, jt) = (X_i1(jr, jt) - X_i1(jr, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dYdh(  jr, jt) = (Y_i1(jr, jt) - Y_i1(jr, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dXdhi1(jr, jt) = (X(jr, jthe_r) - X(jr, jt))/(thetap(jt+1) - thetap(jt))  
    dYdhi1(jr, jt) = (Y(jr, jthe_r) - Y(jr, jt))/(thetap(jt+1) - thetap(jt))
enddo

Jcbn   = dXdr  *dYdh   - dXdh  *dYdr    ! i, j
Jcbni1 = dXdri1*dYdhi1 - dXdhi1*dYdri1  ! i, j+1/2  

do jt=1, Nt
    do jr = 1, Nr-1
        Mdet =  dXdr(jr, jt)*dYdh(jr, jt) - dXdh(jr, jt)*dYdr(jr, jt)
        drdX =  dYdh(jr, jt)/Mdet     ! dlambda/dR
        dhdX = -dYdr(jr, jt)/Mdet     ! dtheta/dR
        drdY = -dXdh(jr, jt)/Mdet     ! dlambda/dZ
        dhdY =  dXdr(jr, jt)/Mdet     ! dtheta/dZ 
        dphdR(jr, jt) = drdX
        dphdZ(jr, jt) = drdY
        gradr(jr, jt) = drdX**2.0 + drdY**2.0   ! ]grad(lambda)]^2
        gradh(jr, jt) = dhdX**2.0 + dhdY**2.0   ! grad(theta)^2
        grt(  jr, jt) = drdX*dhdX + drdY*dhdY   !   grad(theta)*grad(lambda)  ! i, j

        Mdet =  dXdr2(jr, jt)*dYdh2(jr, jt) - dXdh2(jr, jt)*dYdr2(jr, jt)
        drdX =  dYdh2(jr, jt)/Mdet
        dhdX = -dYdr2(jr, jt)/Mdet
        drdY = -dXdh2(jr, jt)/Mdet
        dhdY =  dXdr2(jr, jt)/Mdet
        gradr2(jr, jt) = drdX**2.0 + drdY**2.0 
        gradh2(jr, jt) = dhdX**2.0 + dhdY**2.0
        grt2(  jr, jt) = drdX*dhdX + drdY*dhdY   ! i+1/2, j
  
        Mdet = gradr2(jr, jt)*gradh2(jr, jt) - grt2(jr, jt)**2.0
        ghht2(jr, jt) = gradr2(jr, jt)/Mdet

        Mdet = dXdri1(jr, jt)*dYdhi1(jr, jt) - dXdhi1(jr, jt)*dYdri1(jr, jt)
        drdX = dYdhi1(jr, jt)/Mdet
        dhdX = -dYdri1(jr, jt)/Mdet
        drdY = -dXdhi1(jr, jt)/Mdet
        dhdY = dXdri1(jr, jt)/Mdet
        gradhi1(jr, jt) = dhdX**2.0 + dhdY**2.0  
        grti1(  jr, jt) = drdX*dhdX + drdY*dhdY       ! i, j+1/2
    enddo
enddo

do jt=1, Nt    
    gradr(  1, jt) = 0. 
    dphdR(  1, jt) = 0. 
    dphdZ(  1, jt) = 0. 
    gradh(  1, jt) = 0. 
    gradhi1(1, jt) = 0. 
    grt(    1, jt) = 0. 
    grti1(  1, jt) = 0. 
enddo
   
!Define major radius

do jt=1, Nt
    do jr=1, Nr
        Rmaj(  jr, jt) = X(   jr, jt) 
        Rmaj2( jr, jt) = X2(  jr, jt) 
        Rmaji( jr, jt) = X_i( jr, jt)
        Rmaji1(jr, jt) = X_i1(jr, jt)
    enddo
enddo

! Area at grid points

do jr=1, Nr
    if (jr == 1) then
        dpsi = 0.5*psigp(jr)
    else
        dpsi = (psigp(jr) - psigp(jr-1))
    endif
    do jt=1, Nt
        if (jt == 1) then
            dthe = (thetap_i(Nt+1) - thetap_i(Nt))
        else
            dthe = (thetap_i(jt) - thetap_i(jt-1))
        endif
        dArea(jr, jt) = Jcbn(jr, jt)*dpsi*dthe
    enddo
enddo

! Half grid points

do jr=1, Nr-1
    dpsi = (psig(jr+1) - psig(jr))
    do jt=1, Nt
        if (jt == 1) then
            dthe = (thetap_i(Nt+1) - thetap_i(Nt))
        else
            dthe = (thetap_i(jt) - thetap_i(jt-1))
        endif
        dArea2(jr, jt) = Jcbn2(jr, jt)*dpsi*dthe
    enddo
enddo

!Compute total area, note that this area, in the limit of large Nr, Nt, coincides with pi r^2 if it is
!circle, so is ok!

Area = SUM(dArea)

!Now compute arc lengths do fluxes
 
jr = 1
do jt=1, Nt
    if (jt == 1) then
        dthe = thetap_i(Nt+1) - thetap_i(Nt)
    else
        dthe = thetap_i(jt) - thetap_i(jt-1)
    endif
    dArc_rp1 (jr, jt) = Jcbn2(jr, jt)*gradr2(jr, jt)*dthe/Rmaj2(jr, jt)
    dArc_rpt1(jr, jt) = Jcbn2(jr, jt)*grt2(  jr, jt)*dthe/Rmaj2(jr, jt)
    dArc_rm1(jr, jt) = 0.   ! i-1/2, j
    dArc_tp1(jr, jt) = 0.   ! i, j+1/2
    dArc_tm1(jr, jt) = 0.   ! i, j-1/2
enddo

do jr=2, Nr
    dpsi = psigp(jr) - psigp(jr-1)
    do jt=1, Nt
        if (jt == 1) then
            jthe_l = Nt
        else
            jthe_l = jt-1
        endif
        dthe = thetap_i(jthe_l+1) - thetap_i(jthe_l)
        dArc_rp1( jr, jt) = Jcbn2(jr  , jt)*gradr2(jr  , jt)*dthe/Rmaj2(jr  , jt)
        dArc_rm1( jr, jt) = Jcbn2(jr-1, jt)*gradr2(jr-1, jt)*dthe/Rmaj2(jr-1, jt)
        dArc_rpt1(jr, jt) = Jcbn2(jr  , jt)*grt2(  jr  , jt)*dthe/Rmaj2(jr  , jt)
        dArc_rmt1(jr, jt) = Jcbn2(jr-1, jt)*grt2(  jr-1, jt)*dthe/Rmaj2(jr-1, jt)
        dArc_tp1( jr, jt) = Jcbni1(jr, jt    )*gradhi1(jr, jt    )*dpsi/Rmaji1(jr, jt    )
        dArc_tm1( jr, jt) = Jcbni1(jr, jthe_l)*gradhi1(jr, jthe_l)*dpsi/Rmaji1(jr, jthe_l)
        dArc_tpr1(jr, jt) = Jcbni1(jr, jt    )*grti1(  jr, jt    )*dpsi/Rmaji1(jr, jt    )
        dArc_tmr1(jr, jt) = Jcbni1(jr, jthe_l)*grti1(  jr, jthe_l)*dpsi/Rmaji1(jr, jthe_l)
    enddo
enddo

! Compute differentials

do jr=1, Nr-1
    dpsi = psig(jr+1) - psig(jr)
    do jt=1, Nt
        ddr(jr, jt) = dpsi
    enddo
enddo

do jr=2, Nr
    dpsi = psigp(jr) - psigp(jr-1)
    do jt=1, Nt
        ddr_i(jr, jt) = dpsi
    enddo
enddo

do jr=1, Nr
    do jt=1, Nt
        jthe_l = jt - 1
        if (jt == 1) then
            jthe_l = Nt
        endif 
        dtp( jr, jt) = thetap(jt+1) - thetap(jt)
        dtm( jr, jt) = thetap  (jthe_l+1) - thetap  (jthe_l)
        dt_i(jr, jt) = thetap_i(jthe_l+1) - thetap_i(jthe_l)
    enddo
enddo

return
end subroutine build_2dgrid
