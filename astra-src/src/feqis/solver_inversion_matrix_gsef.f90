subroutine solver_inversion_matrix_gsef(PSIb, nrho, ntheta, &
     known_term, Ndims, LDAB, &
     dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
     dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1, &
     ddr, ddr_i, dtp, dtm, dt_i, PSI)

implicit none

integer, intent(in) :: ntheta, nrho, Ndims, LDAB
double precision, intent(in) :: Psib
double precision, intent(in), dimension(nrho) :: ddr, ddr_i
double precision, intent(in), dimension(ntheta) :: dtp, dtm, dt_i
double precision, intent(in), dimension(nrho, ntheta) :: known_term, &
   dArc_rp1, dArc_rm1, dArc_rpt1, dArc_rmt1, &
   dArc_tp1, dArc_tm1, dArc_tpr1, dArc_tmr1
double precision, intent(out), dimension(nrho, ntheta) :: PSI

integer :: jthe, jrho, NRHS, KL, KU, INFO, iilow, iiup, &
   jdim00, jdim0p, jdimp0, jdim0m, jdimm0, jdimmm, &
   jdimmp, jdimpm, jdimpp, jtm, jtp, jrawm
integer, dimension(Ndims) :: IPIV
double precision :: rm_arc, tp_arc, tm_arc, rp_arc, denom
double precision, dimension(Ndims, 1) :: BB
double precision, dimension(LDAB, Ndims) :: AB

NRHS = 1
KL   = 2*ntheta
KU   = 2*ntheta
AB   = 0.0
BB   = 0.0
PSI  = 0.0

jrho = 1
do jthe=1, ntheta
    jtm = jthe - 1
    jtp = jthe + 1
    if (jthe == 1) then
        jtm = ntheta
    else if (jthe == ntheta) then
        jtp = 1
    endif
    jdim00 = 1
    jdimp0 = 1 + jthe
    jdimpp = 1 + jtp
    jdimpm = 1 + jtm
    jrawm  = KL + KU + 1 + jdim00

    BB(1, 1) = BB(1, 1) + known_term(1, jthe)

    iilow = max(    1, 1 - KU)
    iiup  = min(Ndims, 1 + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdim00, jdim00) = AB(jrawm - jdim00, jdim00) - &
             dArc_rp1(jrho, jthe)/ddr(jrho)
    endif
    iilow = max(    1, jdimp0 - KU)
    iiup  = min(Ndims, jdimp0 + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimp0, jdimp0) = AB(jrawm - jdimp0, jdimp0) +  &
           dArc_rp1(jrho, jthe)/ddr(jrho)
    endif
    iilow = max(    1, jdimpp - KU)
    iiup  = min(Ndims, jdimpp + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimpp, jdimpp) = AB(jrawm - jdimpp, jdimpp) +  &
           0.25*dArc_rpt1(jrho, jthe)/dt_i(jthe)
    endif
    iilow = max(1, jdimpm-KU)
    iiup = min(Ndims, jdimpm + KL)
    if (jdim00 >= iilow .and. jdim00 <= iiup) then
        AB(jrawm - jdimpm, jdimpm) = AB(jrawm - jdimpm, jdimpm) - &
            0.25*dArc_rpt1(jrho, jthe)/dt_i(jthe)
    endif
enddo

do jrho=2, nrho-1

    do jthe=1, ntheta

        jtm = jthe - 1
        jtp = jthe + 1
        if (jthe == 1) then
            jtm = ntheta
        else if (jthe == ntheta) then
           jtp = 1
        endif
        jdimm0 = 1 + jthe + (jrho - 3)*ntheta
        jdim00 = 1 + jthe + (jrho - 2)*ntheta
        jdimp0 = 1 + jthe + (jrho - 1)*ntheta
        jdimmm = 1 + jtm  + (jrho - 3)*ntheta
        jdim0m = 1 + jtm  + (jrho - 2)*ntheta
        jdimpm = 1 + jtm  + (jrho - 1)*ntheta
        jdimmp = 1 + jtp  + (jrho - 3)*ntheta
        jdim0p = 1 + jtp  + (jrho - 2)*ntheta
        jdimpp = 1 + jtp  + (jrho - 1)*ntheta
        jrawm  = KL + KU + 1 + jdim00

        rm_arc = -dArc_tmr1(jrho, jthe)/ddr_i(jrho)
        rp_arc =  dArc_tpr1(jrho, jthe)/ddr_i(jrho)
        tm_arc = -dArc_rmt1(jrho, jthe)/dt_i(jthe)
        tp_arc =  dArc_rpt1(jrho, jthe)/dt_i(jthe)
        denom  = -(dArc_rp1(jrho, jthe)/ddr(jrho) + dArc_rm1(jrho, jthe)/ddr(jrho-1) + &
                   dArc_tp1(jrho, jthe)/dtp(jthe) + dArc_tm1(jrho, jthe)/dtm(jthe))

        iilow = max(    1, jdim00 - KU)
        iiup  = min(Ndims, jdim00 + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim00, jdim00) = denom
        endif
        iilow = max(    1, jdim0p - KU)
        iiup  = min(Ndims, jdim0p + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim0p, jdim0p) = dArc_tp1(jrho, jthe)/dtp(jthe) + &
                0.25*(tp_arc + tm_arc)
        endif
        iilow = max(    1, jdim0m - KU)
        iiup  = min(Ndims, jdim0m + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            AB(jrawm - jdim0m, jdim0m) = dArc_tm1(jrho, jthe)/dtm(jthe) + &
                0.25*(-tp_arc - tm_arc)
        endif

        if (jrho == nrho - 1) then
            BB(jdim00, 1) = known_term(jrho, jthe) - dArc_rp1(jrho, jthe)/ddr(jrho)*PSIb - &
                0.5*(rp_arc + rm_arc)*PSIb
        else
            BB(jdim00, 1) = known_term(jrho, jthe)

            iilow = max(    1, jdimp0 - KU)
            iiup  = min(Ndims, jdimp0 + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimp0, jdimp0) = dArc_rp1(jrho, jthe)/ddr(jrho) + &
                    0.25*(rp_arc + rm_arc)
            endif
            iilow = max(    1, jdimpp - KU)
            iiup  = min(Ndims, jdimpp + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimpp, jdimpp) = 0.25*(rp_arc + tp_arc)
            endif
            iilow = max(    1, jdimpm - KU)
            iiup  = min(Ndims, jdimpm + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimpm, jdimpm) = 0.25*(-tp_arc + rm_arc)
            endif
        endif

        if (jrho == 2) then
            jdimm0 = 1
        else
            iilow = max(    1, jdimmp - KU)
            iiup  = min(Ndims, jdimmp + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimmp, jdimmp) = 0.25*(tm_arc - rp_arc)
            endif

            iilow = max(    1, jdimmm - KU)
            iiup  = min(Ndims, jdimmm + KL)
            if (jdim00 >= iilow .and. jdim00 <= iiup) then
                AB(jrawm - jdimmm, jdimmm) = 0.25*(-rm_arc - tm_arc)
            endif
        endif
        iilow = max(    1, jdimm0 - KU)
        iiup  = min(Ndims, jdimm0 + KL)
        if (jdim00 >= iilow .and. jdim00 <= iiup) then
            if (jrho == 2) then
                AB(jrawm - jdimm0, jdimm0) = dArc_rm1(jrho, jthe)/ddr(jrho-1) + &
                    0.5 *(-rp_arc - rm_arc)
            else
                AB(jrawm - jdimm0, jdimm0) = dArc_rm1(jrho, jthe)/ddr(jrho-1) + &
                    0.25*(-rp_arc - rm_arc)
            endif
        endif

    enddo
enddo

call DGBSV(Ndims, KL, KU, NRHS, AB, LDAB, IPIV, BB, Ndims, INFO)

PSI(1, :) = BB(1, 1)
do jthe=1, ntheta
    do jrho=2, nrho-1
        jdim00 = 1 + jthe + (jrho - 2)*ntheta
        PSI(jrho, jthe) = BB(jdim00, 1)
    enddo
enddo
PSI(nrho, :) = PSIb

end subroutine solver_inversion_matrix_gsef
