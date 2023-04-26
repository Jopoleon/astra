subroutine build_2dgrid(nrho, ntheta, Rb, Zb, X0, Y0, lambda2d, lambda2dp, &
    psig, psiax, psib, j_ok, psigp, PSI, &
! Output
    dArea, X2, dArea2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, X, Y, &
    rmin, thetap, thetap_i, Jcbn, gradr2, Jcbn2)

implicit none

double precision, parameter :: GPI=3.141592653589793, GPI2=2.*GPI

integer, intent(in) :: ntheta, nrho, j_ok
double precision, intent(in) :: psiax, psib, X0, Y0
double precision, intent(in), dimension(ntheta) :: Rb, Zb
double precision, intent(in), dimension(nrho) :: psig, psigp
double precision, intent(in), dimension(nrho, ntheta) :: psi, lambda2d, lambda2dp

double precision, intent(out), dimension(ntheta+1) :: thetap, thetap_i
double precision, intent(out), dimension(nrho, ntheta) :: dArea, dArea2, X2, &
    dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
    dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
    ddr, ddr_i, dtp, dtm, dt_i, X, Y, &
    rmin, Jcbn, Jcbn2, &
    gradr2

integer :: jrho, jthe, jthe_l, jthe_r, j, k
double precision :: drdX, dhdX, drdY, dhdY, Mdet_inv, dpsi, dthe
double precision, dimension(ntheta) :: dXb0, dXb0_i
double precision, dimension(nrho, ntheta) :: lambda2dold, &
    lambda2di, lambda2dpi, psin, Y2, X_i, Y_i, &
    dXdr2, dYdr2, dXdr, dYdr, dXdri1, dYdri1, &
    dXdh2, dYdh2, dXdh, dYdh, dXdhi1, dYdhi1, &
    Jcbni1, gradhi1, grti1, &
    rr2, r_i, r_i1, X_i1, Y_i1, grt2

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
 
do jthe=1, ntheta-1
    dXb0_i(jthe) = 0.5*(dXb0(jthe) + dXb0(jthe+1))
enddo
dXb0_i(ntheta) = dXb0_i(1)

if (j_ok == 1) then !relambda
!relambda
    lambda2dold(1: nrho, 1: ntheta) = lambda2d(1: nrho, 1: ntheta)
    do jthe=1, ntheta
        do jrho=1, nrho
            psin(jrho, jthe) = (psi(jrho, jthe) - psiax)/(psib - psiax)
        enddo
        psin( 1, jthe) = 0
        psin(nrho, jthe) = 1.

        call linterp_feqis(psin(:, jthe), lambda2dold(:, jthe), nrho, &
            psig (2: nrho-1), lambda2d (2: nrho-1, jthe), nrho-2)
        call linterp_feqis(psin(:, jthe), lambda2dold(:, jthe), nrho, &
            psigp(1: nrho-1), lambda2dp(1: nrho-1, jthe), nrho-1)
    enddo
endif 

do jthe=1, ntheta
    do jrho=1, nrho
        rmin(jrho, jthe) = lambda2d (jrho, jthe)*dXb0(jthe)
        rr2 (jrho, jthe) = lambda2dp(jrho, jthe)*dXb0(jthe)
        X   (jrho, jthe) = X0 + rmin(jrho, jthe)*cos(thetap(jthe))
        Y   (jrho, jthe) = Y0 + rmin(jrho, jthe)*sin(thetap(jthe))
        X2  (jrho, jthe) = X0 + rr2 (jrho, jthe)*cos(thetap(jthe))
        Y2  (jrho, jthe) = Y0 + rr2 (jrho, jthe)*sin(thetap(jthe))
    enddo
enddo

do jthe=1, ntheta-1
    lambda2di (:, jthe) = 0.5*(lambda2d (:, jthe) + lambda2d (:, jthe+1))
    lambda2dpi(:, jthe) = 0.5*(lambda2dp(:, jthe) + lambda2dp(:, jthe+1))
enddo
jthe = ntheta
lambda2di (:, jthe) = 0.5*(lambda2d (:, jthe) + lambda2d (:, 1))
lambda2dpi(:, jthe) = 0.5*(lambda2dp(:, jthe) + lambda2dp(:, 1))

do jthe=1, ntheta
    do jrho = 1, nrho
        r_i1(jrho, jthe) = lambda2di (jrho, jthe)*dXb0_i(jthe)
        r_i (jrho, jthe) = lambda2dpi(jrho, jthe)*dXb0_i(jthe)
        X_i1(jrho, jthe) = X0 + r_i1(jrho, jthe)*cos(thetap_i(jthe))
        Y_i1(jrho, jthe) = Y0 + r_i1(jrho, jthe)*sin(thetap_i(jthe))
        X_i (jrho, jthe) = X0 + r_i (jrho, jthe)*cos(thetap_i(jthe))
        Y_i (jrho, jthe) = Y0 + r_i (jrho, jthe)*sin(thetap_i(jthe))
    enddo
enddo
   
! Now computes jacobian
! J = det( {dX/dr dX/dtheta}  {dY/dr dY/dtheta} ) = dX/dr*dY/dtheta - dX/dtheta*dY/dr= J(r, theta)

do jrho=1, nrho-1
    dpsi = psig(jrho+1) - psig(jrho)
    do jthe=1, ntheta
        if (jthe == 1) then
            jthe_l = ntheta
        else
            jthe_l = jthe - 1
        endif
        dthe = thetap_i(jthe_l+1) - thetap_i(jthe_l) 
        dXdr2(jrho, jthe) = (X(jrho+1, jthe) - X(jrho, jthe))/dpsi 
        dYdr2(jrho, jthe) = (Y(jrho+1, jthe) - Y(jrho, jthe))/dpsi
        dXdh2(jrho, jthe) = (X_i(jrho, jthe) - X_i(jrho, jthe_l))/dthe
        dYdh2(jrho, jthe) = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l))/dthe
    enddo
enddo

Jcbn2 = dXdr2*dYdh2 - dXdh2*dYdr2   ! i+1/2, j 

do jrho=2, nrho-1
    dpsi = psigp(jrho) - psigp(jrho-1)
    do jthe=1, ntheta
        jthe_l = jthe - 1
        jthe_r = jthe + 1
        if (jthe == 1) then
            jthe_l = ntheta
        elseif (jthe == ntheta) then
            jthe_r = 1
        endif 
        dXdr  (jrho, jthe) = (X2  (jrho, jthe) - X2 (jrho-1, jthe))/dpsi
        dYdr  (jrho, jthe) = (Y2  (jrho, jthe) - Y2 (jrho-1, jthe))/dpsi
        dXdri1(jrho, jthe) = (X_i (jrho, jthe) - X_i(jrho-1, jthe))/dpsi
        dYdri1(jrho, jthe) = (Y_i (jrho, jthe) - Y_i(jrho-1, jthe))/dpsi
        dXdh  (jrho, jthe) = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
        dYdh  (jrho, jthe) = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
        dXdhi1(jrho, jthe) = (X (jrho, jthe_r) - X(jrho, jthe))/(thetap(jthe+1) - thetap(jthe))  
        dYdhi1(jrho, jthe) = (Y (jrho, jthe_r) - Y(jrho, jthe))/(thetap(jthe+1) - thetap(jthe))
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
    dXdr  (jrho, jthe) = (X2 (jrho, jthe) - X(jrho, jthe))/ psigp(jrho)
    dYdr  (jrho, jthe) = (Y2 (jrho, jthe) - Y(jrho, jthe))/ psigp(jrho)
    dXdri1(jrho, jthe) = (X_i(jrho, jthe) - X(jrho, jthe))/ psigp(jrho)
    dYdri1(jrho, jthe) = (Y_i(jrho, jthe) - Y(jrho, jthe))/ psigp(jrho)
    dXdh  (jrho, jthe) = (X_i(jrho, jthe) - X_i(jrho, jthe_l))/ (thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dYdh  (jrho, jthe) = (Y_i(jrho, jthe) - Y_i(jrho, jthe_l))/ (thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dXdhi1(jrho, jthe) = (X2(jrho, jthe_r) - X2(jrho, jthe))/ (thetap(jthe+1) - thetap(jthe))  
    dYdhi1(jrho, jthe) = (Y2(jrho, jthe_r) - Y2(jrho, jthe))/ (thetap(jthe+1) - thetap(jthe))
enddo

jrho = nrho
dpsi = psigp(jrho) - psigp(jrho-1)
do jthe=1, ntheta
    jthe_l = jthe - 1
    jthe_r = jthe + 1
    if (jthe == 1) then
        jthe_l = ntheta
    elseif (jthe == ntheta) then
        jthe_r = 1
    endif
    dXdr  (jrho, jthe) = (X   (jrho, jthe) - X2 (jrho-1, jthe))/dpsi
    dYdr  (jrho, jthe) = (Y   (jrho, jthe) - Y2 (jrho-1, jthe))/dpsi
    dXdri1(jrho, jthe) = (X_i1(jrho, jthe) - X_i(jrho-1, jthe))/dpsi
    dYdri1(jrho, jthe) = (Y_i1(jrho, jthe) - Y_i(jrho-1, jthe))/dpsi
    dXdh  (jrho, jthe) = (X_i1(jrho, jthe) - X_i1(jrho, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dYdh  (jrho, jthe) = (Y_i1(jrho, jthe) - Y_i1(jrho, jthe_l))/(thetap_i(jthe_l+1) - thetap_i(jthe_l))
    dXdhi1(jrho, jthe) = (X(jrho, jthe_r) - X(jrho, jthe))/(thetap(jthe+1) - thetap(jthe))  
    dYdhi1(jrho, jthe) = (Y(jrho, jthe_r) - Y(jrho, jthe))/(thetap(jthe+1) - thetap(jthe))
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

do jthe=1, ntheta    
    gradhi1(1, jthe) = 0. 
    grti1  (1, jthe) = 0. 
enddo

! Area at grid points

do jrho=1, nrho
    if (jrho == 1) then
        dpsi = 0.5*psigp(jrho)
    else
        dpsi = (psigp(jrho) - psigp(jrho-1))
    endif
    do jthe=1, ntheta
        if (jthe == 1) then
            dthe = (thetap_i(ntheta+1) - thetap_i(ntheta))
        else
            dthe = (thetap_i(jthe) - thetap_i(jthe-1))
        endif
        dArea(jrho, jthe) = Jcbn(jrho, jthe)*dpsi*dthe
    enddo
enddo

! Half grid points

do jrho=1, nrho-1
    dpsi = (psig(jrho+1) - psig(jrho))
    do jthe=1, ntheta
        if (jthe == 1) then
            dthe = (thetap_i(ntheta+1) - thetap_i(ntheta))
        else
            dthe = (thetap_i(jthe) - thetap_i(jthe-1))
        endif
        dArea2(jrho, jthe) = Jcbn2(jrho, jthe)*dpsi*dthe
    enddo
enddo

!Now compute arc lengths do fluxes
 
jrho = 1
do jthe=1, ntheta
    if (jthe == 1) then
        dthe = thetap_i(ntheta+1) - thetap_i(ntheta)
    else
        dthe = thetap_i(jthe) - thetap_i(jthe-1)
    endif
    dArc_rp1 (jrho, jthe) = Jcbn2(jrho, jthe)*gradr2(jrho, jthe)*dthe/X2(jrho, jthe)
    dArc_rpt1(jrho, jthe) = Jcbn2(jrho, jthe)*grt2  (jrho, jthe)*dthe/X2(jrho, jthe)
    dArc_rm1(jrho, jthe) = 0.   ! i-1/2, j
    dArc_tp1(jrho, jthe) = 0.   ! i, j+1/2
    dArc_tm1(jrho, jthe) = 0.   ! i, j-1/2
enddo

do jrho=2, nrho
    dpsi = psigp(jrho) - psigp(jrho-1)
    do jthe=1, ntheta
        if (jthe == 1) then
            jthe_l = ntheta
        else
            jthe_l = jthe-1
        endif
        dthe = thetap_i(jthe_l+1) - thetap_i(jthe_l)
        dArc_rp1 (jrho, jthe) = Jcbn2(jrho  , jthe)*gradr2(jrho  , jthe)*dthe/X2(jrho  , jthe)
        dArc_rm1 (jrho, jthe) = Jcbn2(jrho-1, jthe)*gradr2(jrho-1, jthe)*dthe/X2(jrho-1, jthe)
        dArc_rpt1(jrho, jthe) = Jcbn2(jrho  , jthe)*grt2  (jrho  , jthe)*dthe/X2(jrho  , jthe)
        dArc_rmt1(jrho, jthe) = Jcbn2(jrho-1, jthe)*grt2  (jrho-1, jthe)*dthe/X2(jrho-1, jthe)
        dArc_tp1 (jrho, jthe) = Jcbni1(jrho, jthe  )*gradhi1(jrho, jthe  )*dpsi/X_i1(jrho, jthe  )
        dArc_tm1 (jrho, jthe) = Jcbni1(jrho, jthe_l)*gradhi1(jrho, jthe_l)*dpsi/X_i1(jrho, jthe_l)
        dArc_tpr1(jrho, jthe) = Jcbni1(jrho, jthe  )*grti1  (jrho, jthe  )*dpsi/X_i1(jrho, jthe  )
        dArc_tmr1(jrho, jthe) = Jcbni1(jrho, jthe_l)*grti1  (jrho, jthe_l)*dpsi/X_i1(jrho, jthe_l)
    enddo
enddo

! Compute differentials

do jrho=1, nrho-1
    dpsi = psig(jrho+1) - psig(jrho)
    do jthe=1, ntheta
        ddr(jrho, jthe) = dpsi
    enddo
enddo

do jrho=2, nrho
    dpsi = psigp(jrho) - psigp(jrho-1)
    do jthe=1, ntheta
        ddr_i(jrho, jthe) = dpsi
    enddo
enddo

do jrho=1, nrho
    do jthe=1, ntheta
        jthe_l = jthe - 1
        if (jthe == 1) then
            jthe_l = ntheta
        endif
        dtp (jrho, jthe) = thetap  (jthe  +1) - thetap  (jthe)
        dtm (jrho, jthe) = thetap  (jthe_l+1) - thetap  (jthe_l)
        dt_i(jrho, jthe) = thetap_i(jthe_l+1) - thetap_i(jthe_l)
    enddo
enddo

return
end subroutine build_2dgrid
