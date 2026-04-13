module nbibce

use nbstatus, only: n_rho, n_energy
  
implicit none

integer, parameter :: IV1=321, n_ch=9, n_fun=50

integer :: IV, IT
double precision :: HV, HM, DT, EZ
double precision, dimension(n_energy) :: YFI, CNSFI
double precision, dimension(n_ch) :: VB, RNB, ZB, RMB, EB
double precision, dimension(IV1) :: DV2, A, B, A1, B1, D, AL, BT, &
   AE, BE, AI, BI
double precision, dimension(n_fun) :: YM1, YM2
double precision, dimension(n_rho) :: PBCX, YFCUR, YLNI, YLNE, YLNZ
double precision, dimension(IV1, n_fun) :: FSRS, FVM, RMN

contains

! Subroutines for time dependent Fokker-Planck calculations

!---------------------------------------------------------------------
    double precision function YERF(x_in)
    
    double precision, intent(in) :: x_in
    double precision :: EX2, T, T2, T3

    EX2 = dEXP(-x_in**2)
    T  = 1.d0 / (1.d0 + 0.3275911d0 * x_in)
    T2 = T**2
    T3 = T2 * T

    YERF = (1.d0 - (0.254829592d0 * T - &
        0.284496736d0 * T2 + 1.421413741d0 * T3 - &
        1.453152027d0 * T2 * T2 + 1.061405429d0 * T2 * T3) * EX2)

    end function YERF
    
!---------------------------------------------------------------------
    double precision function YCERF(x_in)
!----------------------------------------------------- 17.11.89
! YcERF =( exp(-x2) / x / sqrt(pi) + erf(x)*(1 - .5/x2) ) / (2 x)
!-------------------------------------------------------------

    double precision, intent(in) :: x_in
    double precision :: X2, EX2

    if (x_in < 0.06d0) then
        YCERF = 0.3761264d0
        return
    endif

    X2 = x_in**2
    EX2 = dEXP(-X2)
    YCERF = (0.56418959d0 * EX2 + YERF(x_in) * (x_in - 0.5d0 / x_in)) / (2.d0 * X2)

    end function YCERF

!---------------------------------------------------------------------
    double precision function YAERF(x_in)
!---------------------------------------------------------------------
! YAERF = erf(x) - 2 * x * exp(-x2) / sqrt(pi)
!
! ERF = integral exp(-s2) ds * 2 / sqrt(pi)
! Abramowitz, Stegun p.122
!---------------------------------------------------------------------

    double precision, intent(in) :: x_in
    double precision :: EX2

    if (x_in < 0.055d0) then
        YAERF = 0.7522528d0 * x_in**3
        return
    endif

    EX2 = dexp(-x_in**2)
    YAERF = YERF(x_in) - 1.1283792d0 * x_in * EX2

    end function YAERF

!---------------------------------------------------------------------
    subroutine NBCOEF
! Fokker-Planck coefficients

    use nbstatus, only: ISPEND

    integer :: I, J, I1

    double precision :: VB2, YE, YA, YB, YD, YI, YDEL, YEXPI, YHDVB, YIP, YAB

    do I=1, IV1
        A1(I) = 0.
        B(I)  = 0.
        A(I)  = ((I + 1) * I / (HV * (I + 0.5))**3)
        B1(I) = 0.
        D(I)  = 0.
    enddo

! Coefficients for plasma species Aelectron, Ap, Ad

    do J=1, ISPEND

        VB2  = VB(J) * VB(J)
        YE   = HV * HV / EB(J)
        YDEL = ZB(J)**2 * RNB(J) / RNB(1)
        YB = .5 * VB2 * YDEL * dexp(0.75 * YE)
        YA = dexp(YE)
        YD = YDEL / (HM * HM) / VB(J)
        YHDVB = HV / VB(J)

        do I=1, IV1
            I1 = I + 1
            YI  = I * YHDVB
            YIP = (I + .5d0) * YHDVB
            YEXPI = dexp(I * YE)
            YAB = YB * A(I) * YAERF(YIP)
            B(I)  = B(I)  + YAB * YEXPI
            A1(I) = A1(I) + YAB / (YEXPI * YA)
            D(I)  = D(I) + YD * YCERF(YI)
        enddo

        if (J == 1) then
            do I = 1, IV1
                AE(I) = A1(I)
                BE(I) = B(I)
            enddo
        endif
    enddo

    A(1) = 0.

    do I=1, IV1
        AI(I) = A1(I) - AE(I)
        BI(I) = B(I)  - BE(I)
    enddo

    do I=1, IV
        I1 = I + 1
        D(I)   = D(I)  * DV2(I)
        A(I1)  = A1(I) * DV2(I1)
        A1(I)  = A1(I) * DV2(I)
        B1(I1) = B(I)  * DV2(I1)
        B(I)   = B(I)  * DV2(I)

    enddo

    B(IV1)  = B(IV1)  * DV2(IV1)
    D(IV1)  = D(IV1)  * DV2(IV1)
    A1(IV1) = A1(IV1) * DV2(IV1)

    end subroutine NBCOEF

!---------------------------------------------------------------------
    subroutine NBMESH
! V, MU mesh

    integer :: IT1, JV, JT

! ... Vi = i * HV  i = 1..IV1
! ... MUj = -1 + (j - 1/2) * Hmu

    IT1 = IT + 1
    HM = 2. / IT
    HV = 4. / 3. / IV

! ... YM2j = 1 - (Mj + 1/2)**2
    do JT = 1, IT
        YM1(JT) = -1. + HM * (JT - .5)
        YM2(JT) = JT * HM * (2. - JT * HM)
    enddo

! ... DV2i = 1 / Vi**2
    do JV = 1, IV1
        DV2(JV) = 1. / (HV * JV)**2
    enddo

    end subroutine NBMESH

!---------------------------------------------------------------------
    subroutine NBPOMU
! MU sweep

    integer :: I, ITM1, J, J1, JM
    double precision :: AJ, BJ, CJ, DJ, Y

    ITM1 = IT - 1

! ... Calculations of ALFA, BETA

    do I=1, IV
        BJ = D(I) * YM2(1)
        CJ = BJ + DT
        DJ = FSRS(I, 1) + DT * FVM(I, 1)
        AL(2) = BJ / CJ
        BT(2) = DJ / CJ

        do J=2, ITM1
            J1 = J + 1
            DJ = FSRS(I, J) + DT * FVM(I, J)
            AJ = BJ
            BJ = D(I) * YM2(J)
            CJ = AJ + BJ + DT
            Y = (CJ - AL(J) * AJ)
            AL(J1) = BJ / Y
            BT(J1) = (DJ + AJ * BT(J)) / Y
        enddo

! ... Sweeping
        DJ = FSRS(I, IT) + DT * FVM(I, IT)
        AJ = BJ
        CJ = AJ + DT
        FVM(I, IT) = (DJ + AJ * BT(IT)) / (CJ - AL(IT) * AJ)

        do JM=1, ITM1
            J  = IT - JM
            J1 = J + 1
            FVM(I, J) = AL(J1) * FVM(I, J1) + BT(J1)
        enddo

    enddo

    end subroutine NBPOMU

!---------------------------------------------------------------------
    subroutine NBPOVE
! V sweep

    integer :: J, JP, I, I1, IM
    double precision :: Y

    do J=1, IT ! ..... Calculations of alfa, beta
        JP = J + 1
        Y = DT + A1(1) + RMN(1, J)
        AL(2) = B(1)/Y
        BT(2) = (FSRS(1, J) + FVM(1, J) * DT)/Y

! ... Sweeping
        do I=2, IV
            I1 = I + 1
            Y = DT + A1(I) + B1(I) - AL(I) * A(I) + RMN(I, J)
            AL(I1) = B(I)/Y
            BT(I1) = (A(I) * BT(I) + FSRS(I, J) + FVM(I, J) * DT)/Y
        enddo

        FVM(IV1, J) = BT(IV1) * A1(IV) / (B(IV) - A1(IV) * AL(IV1))

        do IM=1, IV
            I  = IV - IM + 1
            I1 = I + 1
            FVM(I, J) = AL(I1) * FVM(I1, J) + BT(I1)
        enddo

    enddo

    end subroutine NBPOVE

!---------------------------------------------------------------------
    subroutine NBIONR(EBEAM, ABEAM, RTOR, NA1, TAU, NNCL, NNWM, &
      CBM1, CBM3, CBM4, CBMI2, CBMI3, JSRREC, YEXTARR)

!====================================================== 22-MAR-99
! Update 05-NOV-2012, 14-JAN-13
! linearization 29-OCT-2015
! 2D (v, cosTET) Fokker-Planck solver for different
! magnetic surfaces, ripple losses
!======================================================== Polevoy

    use nbstatus, only: PEBM, PIBM, CUBM, CUFI, PBPER, PBLON, &
        NIBM, TE, NE, TI, AMAIN, AMETR, SHIF, ZEF, ISPE, ISPEND, &
        NN, TN, NNBM1, NNBM2, NNBM3
    use cross_sections, only: SPEX, FNBF
    use pi_const, only: GP
    use nbicom, only: YASBA, YASBA1

    integer, intent(in) :: NA1, JSRREC
    double precision, intent(in) :: EBEAM, ABEAM, RTOR, &
        TAU, NNCL, NNWM, CBM1, CBM3, CBM4, CBMI2, CBMI3, YEXTARR(n_rho, 9)

    integer :: JFPBEG, JN1OLD, N, N1, ITC, IT1, NTET1, IEB, J1BEG, J1END, &
        I, J, JT, JV, JN, JE, JBMS4, JNA, JNAC, J2, JSP, ISP, &
        ITRAP, JTDTS, ITIME, NTET, JN22, IE, IVE, J1, JTIME, &
        I1, I2, JLREC, JSRNUM, JDBL

    double precision :: YET, CBMI33, YM2F, &
        SQPI, DTION, YEV21, YEV22, YEV23, YJ2, YEPS, T0, TSNBI, &
        CNSTN, CNSTQ, CNSTP, CNSTC, CNSTE, CNSTT, CNSTQT, &
        DTAU, EBDTI, EBDTI0, CNSTE0, CNSNN, CNSNN0, YNN0, YV2, &
        YRMN, YSRSE, YEBEAM, YPB, YIP, YV4, FVMMIN, F0J, YEXARG, &
        YE, YI, YPEBM, YPIBM, YPBPER, YPBLON, YNB, YCUFI, YMF, &
        PTHBM, PTHERM, YDELPE, YDELPI, YDELMI, YDELME, YMAXMU

    save JFPBEG, JN1OLD

    data JFPBEG /0/
    data JN1OLD /1/

    JDBL = 2
    SQPI = sqrt(GP)

!----- time step
    if (TAU <= 0.) return
    DTION = CBMI2 * TAU

!----- flux surface
    N1 = (NA1 - 1) / CBMI3 + 1
    N  = N1 - 1

!----- pitch angle
    ITC   = 50 / 2
    NTET  = ITC
    NTET1 = ITC + 1
    IT  = 2 * ITC
    IT1 = IT + 1

!----- velocity
    IV = IV1 - 1

!----- beam energy components
    YFI(3) = (IV * 3) / 4.d0
    YFI(2) = YFI(3) / 1.414214d0
    YFI(1) = YFI(3) / 1.732051d0
    IEB = 3

!----- clear previous distributions
    do JN=1, NA1
        PEBM(JN) = 0.
        PIBM(JN) = 0.
        CUBM(JN) = 0.
        CUFI(JN) = 0.
        PBPER(JN) = 0.
        NIBM(JN) = 0.
        PBLON(JN) = 0.
    enddo

!----- file for Fij
    JLREC = JDBL * 4 * IV1 * IT

    open(31, file = 'dat/fij.dat', form='unformatted', access='direct', recl=JLREC)

!----- zero initial distribution
    if (N1 > JN1OLD) then
        J1BEG = JN1OLD
        J1END = N1
    else
        J1END = JN1OLD
        J1BEG = N1
    endif

    if (J1BEG /= J1END) then
        do JN=J1BEG, J1END
            do I=1, IV1
                do J=1, IT
                    FVM(I, J) = 0.
                enddo
            enddo
            write(31, rec=JN, err=999) ((FVM(JV, JT), JV = 1, IV1), JT = 1, IT)
        enddo
    endif

    JFPBEG = 1
    JN1OLD = N1

!----- source arrays
    do JT=1, NTET1
        do JN=1, N1
            do JE=1, 3
                YASBA1(JE, JN, JT) = 0.d0
                YASBA (JE, JN, JT) = 0.d0
            enddo
        enddo
    enddo

    JBMS4 = CBM1

    call NBMESH

    EZ = 0.d0
    CNSTN = HM*HV
    CNSTQ = 3.2d-3*EBEAM*HV**2 * CNSTN
    CNSTP = EBEAM*CNSTN*HV*HV
    CNSTC = 0.7*SQRT(EBEAM/ABEAM)*HV*CNSTN
    CNSTE = 1.6d-3*CNSTP/DTION
    CNSTT = 3.5d-2*EBEAM*SQRT(EBEAM*ABEAM)
    YEV21 = EBEAM/ABEAM
    YEV22 = YEV21/2.
    YEV23 = YEV21/3.
! Change to program units
    YJ2 = 0.
! Radial distribution cycle
    open(35, file='dat/srsfi.dat', form='unformatted', access='DIRECT', recl=JSRREC)

    do JN=1, N ! Output for distribution function        
        JNA = 1 + CBMI3*(JN - 1)
        JNAC = JNA - 1 + CBMI3
        J2 = JNA - CBMI3*YJ2
        if (JN == N) JNAC=NA1-1
        if (CBMI3 > 1.d0) YJ2 = 0.5d0
        YLNE(JN) = 15.85d0 + log(TE(J2)/sqrt(NE(J2)))
        if (EBEAM > 100.*ABEAM) then
            YLNI(JN) = 23.7d0 + log(AMAIN(J2)/(AMAIN(J2) + ABEAM)* &
              sqrt(1.d-3*ABEAM*EBEAM*TE(J2)/NE(J2)))
        else
            YLNI(JN) = 25.4d0 + log(1.d-3*EBEAM*AMAIN(J2)/ &
                  (AMAIN(J2) + ABEAM)*sqrt(TE(J2)/NE(J2)))
        endif
        YLNZ(JN) = YLNI(JN)
        YEPS = AMETR(J2)/(RTOR + SHIF(J2))
        ITRAP = (1.d0 - sqrt(2.*YEPS/(1. + YEPS)))/HM - 1
        YFCUR(JN) = (1. - FNBF(ZEF(J2), YEPS)/ZEF(J2))
        PBCX(JN) = 0.
        RNB(1) = NE(J2)*YLNE(JN)
        EB(1) = TE(J2)/EBEAM
        VB(1) = sqrt(EB(1)*ABEAM/RMB(1))

        do JSP=2, ISPEND
            ISP = ISPE(JSP)
            if (EBEAM > 100.*ABEAM) then
                YLNZ(JN) = 23.7d0 + log(RMB(JSP)/(RMB(JSP) + ABEAM)* &
                 sqrt(1.d-3*ABEAM*EBEAM*TE(J2)/NE(J2)))
            else
                YLNZ(JN) = 25.4d0 + log(1.d-3*EBEAM*RMB(JSP)/ &
             (RMB(JSP)+ABEAM)*SQRT(TE(J2)/NE(J2)))
            endif
            RNB(JSP) = YEXTARR(J2, ISP)*YLNZ(JN)
            EB(JSP) = TI(J2)/EBEAM
            VB(JSP) = sqrt(EB(JSP)*ABEAM/RMB(JSP))
        enddo
  
! Time step DTAU[s]
        T0 = CNSTT/(NE(J2)*YLNE(JN))

! Beam prtcls. slowing down time
        TSNBI = 2.d0*ABEAM*sqrt(TE(J2))*TE(J2)/YLNE(JN)/NE(J2)
        CBMI33 = 2.
        JTDTS = CBMI33*DTION/TSNBI
        if (JTDTS > 1) then
            ITIME = JTDTS
        else
            ITIME = 1
        endif
        DTAU = DTION/ITIME
        DT   = T0/DTAU
! Constants
        CNSTQT = CNSTQ/T0
        EBDTI  = 1.d0/EB(2)
        EBDTI0 = EBDTI/DV2(1)
        CNSTE0 = EB(2)*sqrt(EB(2))*SQPI/2.d0

! Sources and losses distributions
        CNSNN = 4.373E7*T0*CBM3*0.5d0
        CNSNN0 = (NNCL + NNWM)*4.373d7*T0*CBM4*0.5d0
        do JV=1, IV1
            do JT=1, IT
                RMN(JV, JT) = 0.d0
                FSRS(JV, JT) = 0.d0
            enddo
        enddo

! Fast ions CX due to cold neutrals
! Fast ions CX due to NB neutrals

        if (CBM4 > 0.d0 .or. CBM3 > 0.d0) then
            YNN0 = CNSNN0*NN(J2) + CNSNN * &
               (NNBM1(J2) + NNBM2(J2) + NNBM3(J2))
            do JV=1, IV1
                YV2 = YEV21/DV2(JV)
                YET = max(TN(J2), YV2)
                YRMN = YNN0*SPEX(YET)*sqrt(YV2)
                do JT=1, IT
                    RMN(JV, JT) = YRMN
                enddo
            enddo
        endif
! End of CX losses

        YSRSE = 0.
        do JSRNUM=1, JBMS4
            read(35, rec=JSRNUM, err=211) YEBEAM, &
                (((YASBA(JE, JN22, JT), JE=1, 3), JN22=1, N), JT=1, IT)
 211        continue

            do JT=1, IT
                do IE=1, 3
                    if (YASBA(IE, JN, JT) > 0.) then
                        IVE = YFI(IE)*sqrt(YEBEAM/EBEAM)
                        CNSFI(IE) = DV2(IVE)/CNSTN * 0.5 * &
                                 (YFI(IE)/IVE)**2*YEBEAM/EBEAM
! Correction of power balance
                        FSRS(IVE, JT) = YASBA(IE, JN, JT)* &
                         CNSFI(IE)*T0 + FSRS(IVE,JT)
                        YSRSE = YSRSE + YASBA(IE, JN, JT)*T0/CNSTN * &
                             0.5*IVE**2
                    endif
                enddo
            enddo
        enddo 

        YSRSE = YSRSE*CNSTQT
! Coefficients
        call NBCOEF

        read(31, rec=JN, err=996) ((FVM(JV, JT), JV=1, IV1), JT=1, IT)
! For dPb/dt
        YPB = 0.
        do I=1, IV
            YIP = 0.
            YV4 = I**2/DV2(I)
            do J=1, ITC
                J1 = IT - J + 1
                YIP = YIP + FVM(I, J) + FVM(I, J1)
            enddo
            YPB = YPB + YV4*YIP
        enddo

! Fij sweeping
        do JTIME=1, ITIME
            call NBPOMU
            call NBPOVE
        enddo

! Fij linearization Fij b =Fij - Foj exp(-Ei/T)
        FVMMIN = FVM(1, 1)
        do J=2, IT
            FVMMIN = min(FVMMIN, FVM(1, J))
        enddo
        F0J    = FVMMIN
        PTHBM  = F0J*TI(J2)*CNSTE0*1.5
        PTHERM = 1.6d-3*PTHBM/DTION
        YE     = 0.
        YMAXMU = 0.d0

        do JT =1, IT
            YMAXMU = YMAXMU + FVM(1, JT) - F0J
        enddo

        do JV=1, IV1
            YE = 0.
            YEXARG = EBDTI*(1.d0/DV2(1) - 1.d0/DV2(JV))
            YE = dEXP(YEXARG)*F0J
            do JT=1, IT
                FVM(JV, JT) = FVM(JV, JT) - YE
            enddo
        enddo

        write(31, rec=JN, err=996) &
         ((FVM(JV, JT), JV=1, IV1), JT=1, IT)

! Power to plasma, beam pressure, density, current
        YIP    = 0.
        YI     = 0.
        YPEBM  = 0.
        YPIBM  = 0.
        YPBPER = 0.
        YPBLON = 0.
        YNB    = 0.
        YCUFI  = 0.
        YMF    = 0.
        YM2F   = 0.
        do J=1, ITC
            J1 = IT - J + 1
            YM2F = YM2F + YM2(J)*(FVM(1, J) + FVM(1, J1))
            YIP  = YIP + FVM(1, J) + FVM(1, J1)
        enddo
        do J=1, ITRAP
            J1 = IT - J + 1
            YMF = YMF + YM1(J)*(FVM(1, J) - FVM(1, J1))
        enddo
        YNB = YNB + YIP/DV2(1)
        YV4 = 1./DV2(1)
        YPBPER = YPBPER + YV4*YM2F
        YPBLON = YPBLON + YV4*(YIP - YM2F)
        YCUFI = YCUFI + YMF/DV2(1)

        YDELPE = 0.
        YDELPI = 0.
        do I=1, IV
            I1  = I+1
            I2  = I**2
            YI  = YIP
            YIP = 0.
            YMF = 0.
            YM2F = 0.
            YV4 = (I1**2/DV2(I1))
            do J=1, ITC
                J1 = IT - J + 1
                YM2F = YM2F + YM2(J)*(FVM(I1, J) + FVM(I1, J1))
                YIP  = YIP + FVM(I1, J) + FVM(I1, J1)
            enddo

            do J=1, ITRAP
                J1 = IT - J + 1
                YMF = YMF + YM1(J)*(FVM(I1, J) - FVM(I1, J1))
            enddo
            YDELMI = YDELPI
            YDELPI = BI(I)*YIP - AI(I)*YI
            YDELME = YDELPE
            YDELPE = BE(I)*YIP - AE(I)*YI
            YNB    = YNB + YIP/DV2(I1)
            YCUFI  = YCUFI + YMF*I1/DV2(I1)
            YPBPER = YPBPER + YV4*YM2F
            YPBLON = YPBLON + YV4*(YIP - YM2F)
            YPEBM  = YPEBM - I**2 * (YDELPE - YDELME)
            YPIBM  = YPIBM - I**2 * (YDELPI - YDELMI)
        enddo
        do J=JNA, JNAC
! Beam power distribution
            PEBM(J) = YPEBM*CNSTQT/2.	
            PIBM(J) = YPIBM*CNSTQT/2. + PTHERM
! Beam energy Wbeam = Pbper + Pblon/2
! Beam perpendicular pressure <Mb Vort2/2> [10*19 keV/m3]
            PBPER(J) = YPBPER*CNSTP
! Beam parallel pressure   <Mb Vpar2> [10*19 keV/m3]
            PBLON(J) = 2.*YPBLON*CNSTP
! Beam density	     [10*19 /m3]
            NIBM(J) = YNB*CNSTN
! Fast ion current without trapping correction for Eb [MA/m2]
            CUFI(J) = YCUFI*CNSTC
            CUBM(J) = CUFI(J)*YFCUR(JN)
        enddo
    enddo ! radial loop

    close(31)
    close(35)

    return

 999  write(*, *) 'error NBION2'
    return

 996  write(*, *) 'R/W error in NBION2'

    end subroutine NBIONR

end module nbibce
