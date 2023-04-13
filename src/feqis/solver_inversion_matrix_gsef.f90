subroutine solver_inversion_matrix_gsef(PSIb, Nr, Nt, &
     known_term, Ndims, LDAB, &
     dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
     dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
     ddr, ddr_i, dtp, dtm, dt_i, PSI)

implicit none

integer, intent(in) :: Nt, Nr, Ndims, LDAB
double precision, intent(in) :: Psib
double precision, intent(in), dimension(Nr, Nt) :: known_term, &
   dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, dArc_tp1, dArc_tm1, &
   dArc_tpr1, dArc_tmr1, ddr, ddr_i, dtp, dtm, dt_i
double precision, intent(out), dimension(Nr, Nt) :: PSI

integer :: jt, jr, NRHS, KL, KU, INFO, iilow, iiup, &
   jdim00, jdim0p, jdimp0, jdim0m, jdimm0, jdimmm, &
   jdimmp, jdimpm, jdimpp, jtm, jtp, jrawm
integer, dimension(Ndims) :: IPIV
double precision :: rm_arc, tp_arc, tm_arc, rp_arc, denom
double precision, dimension(Ndims, 1) :: BB
double precision, dimension(LDAB, Ndims) :: AB

NRHS = 1
KL = 2*Nt
KU = 2*Nt

AB = 0.0
BB  = 0.0
PSI = 0.0

jr = 1
do jt=1, Nt
    jtm = jt - 1
    jtp = jt + 1
    if (jt == 1) then
        jtm = Nt
    endif
    if (jt == Nt) then
        jtp = 1
    endif

    jdim00 = 1
    jdimp0 = 1 + jt  + (jr-1)*Nt
    jdimpp = 1 + jtp + (jr-1)*Nt
    jdimpm = 1 + jtm + (jr-1)*Nt
    jrawm  = KL + KU + 1 + jdim00

    BB(1, 1) = BB(1, 1) + known_term(1, jt)

    iilow = max(    1, 1 - KU)
    iiup  = min(Ndims, 1 + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdim00, jdim00)  =  AB(jrawm - jdim00, jdim00) - &
             dArc_rp1(1, jt)/ddr(1, jt)
    endif

    iilow = max(    1, jdimp0 - KU)
    iiup  = min(Ndims, jdimp0 + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimp0, jdimp0) = AB(jrawm - jdimp0, jdimp0) +  &
           dArc_rp1(jr, jt)/ddr(jr, jt)
    endif

    iilow = max(    1, jdimpp - KU)
    iiup  = min(Ndims, jdimpp + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimpp, jdimpp) = AB(jrawm - jdimpp, jdimpp) +  &
           1./4.*dArc_rpt1(jr, jt)/dt_i(jr, jt)
    endif

    iilow = max(1, jdimpm-KU)
    iiup = min(Ndims, jdimpm + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimpm, jdimpm) = AB(jrawm - jdimpm, jdimpm)- &
            1./4.*dArc_rpt1(jr, jt)/dt_i(jr, jt)
    endif

enddo

do jr=2, Nr-1

    do jt=1, Nt

        jtm = jt - 1
        jtp = jt + 1
        if (jt == 1) then
            jtm = Nt
        endif
        if (jt == Nt) then
           jtp = 1
        endif
        jdim00 = 1 + jt  + (jr - 2)*Nt
        jdimp0 = 1 + jt  + (jr - 1)*Nt
        jdimm0 = 1 + jt  + (jr - 3)*Nt
        jdimmm = 1 + jtm + (jr - 3)*Nt
        jdimmp = 1 + jtp + (jr - 3)*Nt
        jdim0p = 1 + jtp + (jr - 2)*Nt
        jdim0m = 1 + jtm + (jr - 2)*Nt
        jdimpp = 1 + jtp + (jr - 1)*Nt
        jdimpm = 1 + jtm + (jr - 1)*Nt			
        jrawm  = KL + KU + 1 + jdim00

        rm_arc = -dArc_tmr1(jr, jt)/ddr_i(jr, jtm)
        rp_arc =  dArc_tpr1(jr, jt)/ddr_i(jr, jt)
        tm_arc = -dArc_rmt1(jr, jt)/dt_i(jr-1, jt)
        tp_arc =  dArc_rpt1(jr, jt)/dt_i( jr, jt)
        denom  = -(dArc_rp1(jr, jt)/ddr(jr, jt) + dArc_rm1(jr, jt)/ddr(jr-1, jt) + &
                   dArc_tp1(jr, jt)/dtp(jr, jt) + dArc_tm1(jr, jt)/dtm(jr , jt))

        if (jr == Nr - 1) then
            BB(jdim00, 1) = known_term(jr, jt) - dArc_rp1(jr, jt)/ddr(jr, jt)*PSIb - &
                1./2.*(rp_arc + rm_arc)*PSIb
        else
            BB(jdim00, 1) = known_term(jr, jt)
        endif

        iilow = max(    1, jdim00 - KU)
        iiup  = min(Ndims, jdim00 + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim00, jdim00) = denom
        endif

        iilow = max(    1, jdim0p - KU)
        iiup  = min(Ndims, jdim0p + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim0p, jdim0p) = dArc_tp1(jr, jt)/dtp(jr, jt) + &
                1./4.*(tp_arc + tm_arc)
        endif

        iilow = max(    1, jdim0m - KU)
        iiup  = min(Ndims, jdim0m + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim0m, jdim0m) = dArc_tm1(jr, jt)/dtm(jr, jt) + &
                1./4.*(-tp_arc - tm_arc)
        endif

        if (jr < Nr - 1) then
            iilow = max(    1, jdimp0 - KU)
            iiup  = min(Ndims, jdimp0 + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimp0, jdimp0) = dArc_rp1(jr, jt)/ddr(jr, jt) + &
                    1./4.*(rp_arc + rm_arc)
            endif

            iilow = max(    1, jdimpp - KU)
            iiup  = min(Ndims, jdimpp + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimpp, jdimpp) = 1./4.*(rp_arc + tp_arc)
            endif

            iilow = max(    1, jdimpm - KU)
            iiup  = min(Ndims, jdimpm + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimpm, jdimpm) = 1./4.*(-tp_arc + rm_arc)
            endif
        endif

        if (jr == 2) then
            jdimm0 = 1
        endif
        iilow = max(    1, jdimm0 - KU)
        iiup  = min(Ndims, jdimm0 + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            if (jr == 2) then
                AB(jrawm - jdimm0, jdimm0) = dArc_rm1(jr, jt)/ddr(jr-1, jt) + &
                    1./2.*(-rp_arc - rm_arc)
            else
                AB(jrawm - jdimm0, jdimm0) = dArc_rm1(jr, jt)/ddr(jr-1, jt) + &
                    1./4.*(-rp_arc - rm_arc)
            endif
        endif

        if (jr > 2) then
            iilow = max(    1, jdimmp - KU)
            iiup  = min(Ndims, jdimmp + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimmp, jdimmp) = 1./4.*(tm_arc - rp_arc)
            endif

            iilow = max(    1, jdimmm - KU)
            iiup  = min(Ndims, jdimmm + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimmm, jdimmm) = 1./4.*(-rm_arc - tm_arc)
            endif
        endif
    enddo
enddo

call DGBSV(Ndims, KL, KU, NRHS, AB, LDAB, IPIV, BB, Ndims, INFO)
PSI(1, :) = BB(1, 1)
do jt=1, Nt
    do jr=2, Nr-1
        jdim00 = 1 + jt + (jr - 2)*Nt
        PSI(jr, jt) = BB(jdim00, 1)
    enddo
enddo
PSI(Nr, :) = PSIb

return
end subroutine solver_inversion_matrix_gsef
