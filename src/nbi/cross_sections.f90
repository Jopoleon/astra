module cross_sections

  implicit none

contains

!---------------------------------------------------------------------
    double precision function SPEX(E_in)
!---------------------------------------------- [10#-13 cm2]
! Neutral beam charge-exchange cross section by proton impact
! Reviere A.C //Nucl.Fusion.v.11(1971).p.363
! E [kev] = Ebeam*Mp/Mb
!------------------------------------------- Polevoy A.R. 21.03.91

    double precision, intent(in) :: E_in
    double precision :: Y

    Y = 1.d3*E_in
    SPEX = 0.06937d0*(1.d0 - 0.155d0*LOG10(Y))**2 / (1.d0 + 0.1112d-14 * Y**3.3)

    end function SPEX

!---------------------------------------------------------------------
    double precision function FNBF(Z_in, E_in)
! toroidal correction for NB driven current
! Kim Y.B.,Callen J.D.,Hamnen H.//Nucl.Fus.,(1988)
! Neocl.Cur.And Transp.In Aux.Heat.Tokamaks

    double precision, intent(in) :: Z_in, E_in
    double precision :: FT, G, Z2

    FT = sqrt(E_in)*(1.46d0 - 0.46d0*E_in)
    G = FT/(1.d0 - FT)
    Z2 = Z_in**2
    FNBF = ((Z2 + 1.41d0*Z_in) + (Z2 + 0.45d0*Z_in)*G)/((Z2 + 1.41d0*Z_in) + &
         (2.d0*Z2 + 2.66d0*Z_in + 0.75d0)*G + (Z2 + 1.24*Z_in + 0.35) * G**2)

    end function FNBF

!---------------------------------------------------------------------
    double precision function FNBI1(X2, Y)
!---------------------------------------------------------------------
! aproximation for 1/X((1+X3)/X3)  S(U3/(1+U3))    dU
!    0<Y<4                         0
! Mikkelsen D.R., Singer C.E.//Nucl.Tecnol./Fus., V.4, Sept.1983, PP.237-252
! with Ar stopping crossection by Leonov

    double precision, intent(in) :: X2, Y
    double precision :: X

    X = sqrt(X2)
    FNBI1 = X2*X/(4.d0 + 3.d0*Y + X2*(X + 1.39d0 + 0.61d0*Y**0.7))

    end function FNBI1

!---------------------------------------------------------------------
    double precision function FNB2(X2)
! FNB2 = S udu/(1+u3)

    double precision, intent(in) ::X2
    double precision :: X

    X = sqrt(X2)
    FNB2 = (0.166666667d0*LOG((1.d0 - X + X2)/(1.d0 + 2.d0*X + X2)) + &
        0.57735026d0*(DATAN(0.57735026d0*(2.d0*X - 1.d0)) + 0.52359874d0))

    end function FNB2

!---------------------------------------------------------------------
    double precision function FNBP(X, B)
!----------------------------------------------------
!        3b    X             3b+1     
! FNBP  =(1+1/x3)/x2* S z*(z3/(1+z3)) dz
!              0   step 0.05
!----------------------------Polevoy 18.05.89 -------

    double precision, intent(in) :: X, B

    integer :: JEND, J
    double precision :: YA, A, YS, Z, Z1, Z3, Y, YS1, YS2
    
    YA = 3.d0*B
    A = YA + 1.d0
    JEND = 20.000001*X
    if (JEND < 7) then
        if (X <= 0.0d0) then
            FNBP = 0.0d0
            return
        endif
        FNBP = X**(3.*A + 2.)/(3.*A + 2)*(1.d0 + 1./X**3)**YA/X**2
        return
    endif
    YS = 0.0d0
    Z  = 0.0d0
    if (JEND >= 100) then
        do J=1, 100
            Z =Z + 0.050d0
            Y =Z*(1.d0 - 1.d0/(1.d0 + Z**3))**A
            YS =YS+Y
        enddo
        FNBP = ((YS - 0.5d0*Y)*0.050d0 + 0.5d0*(X - 5.0d0)*(X + 5.0d0)) * &
            (1.d0 + 1.d0/X**3)**YA/X**2
        return
    else
        do J=1, JEND
            Z = Z + 0.050d0
            Z3 = Z**3
            Y = Z*(1.d0 - 1.d0/(1.d0 + Z3))**A
            YS =YS+Y
        enddo
        YS1 = (YS - 0.5d0*Y)*(1.d0+1.d0/Z3)**YA/Z**2
        Z1 = Z
        Z = Z + 0.05d0
        YS2 = (YS + 0.5d0*Z*(1.d0 - 1.d0/(1. + Z**3))**A) * (1.d0+ 1./Z**3)**YA/Z**2
        FNBP = 0.05d0*YS1+(YS2 - YS1)*(X - Z1)
    endif

    end function FNBP

!---------------------------------------------------------------------
    double precision function SIMPI(E_keV, Zq)
!----------------------------------------------- Simpi [10#-13 cm2]
! Neutral beam impurity impact ionization cross section for
! H(D, T) beam (3He, 4He, C, O, Fe) impurity
! Janev R.K., Boley C.D., Post D.E., "Penetration of Neutral
! Beams into Fusion Plasma", Nucl.Fusion, Vol.29., No.12, (1989),
!  pp. 2125-2137.
! Ebeam [KeV]=0-10000,
! Simpi =Zq*C1*[1/(1+C2*E)+C3*ln(1+C5*E)/(C4+E)]
! E [keV] =Eb/Mb/Zq
!--------------------------------------------Polevoy A.R. 21.03.91

    double precision, intent(in) :: E_keV, Zq
    double precision :: E

    E = E_keV/ZQ
    SIMPI = 7.457d-3*ZQ*(1./(1. + 0.08095*E)+2.754*LOG(1. + 1.27*E)/(64.58 + E))

    end function SIMPI

!---------------------------------------------------------------------
    double precision function SEIV(Te_keV)
!------------------------------------------ <S*Ve> [10#-13 cm3/s]
! The electron velocity weighted Maxwellian average for
! neutral beam ionization cross section by electrons
! Te [keV] -  electron temperature
!------------------------------------------- Polevoy A.R. 21.03.91

    double precision, intent(in) :: Te_keV
    double precision :: Y

    Y = 13.6d-3/Te_keV
    if (Te_keV < 0.02) then
        SEIV = 0.61d6*EXP(-Y)/(sqrt((1. + Y)/Y)*(Y + 0.73d0))
    else
        Y = LOG10(Te_keV) + 3.d0
        SEIV = 10.**(7.769d0 - 0.5151d0*Y - 2.563d0/Y)
    endif

    end function SEIV

!---------------------------------------------------------------------
    double precision function SPII(E_keV)
!---------------------------------------------- Spii [10#-13 cm2]
! Neutral beam ionization cross section by proton impact
! Reviere A.C //Nucl.Fusion.v.11(1971).p.363
! E [kev] = Ebeam*Mp/Mb
!------------------------------------------- Polevoy A.R. 21.03.91

    double precision, intent(in) :: E_keV
    double precision :: Y1

    if (E_keV >= 150.d0) then
        SPII = 36.d-3*(LOG10(0.1666d0*E_keV) + 3.d0)/E_keV
    else if (E_keV > 3.d0) then
        Y1 = LOG10(E_keV) + 3.d0
        SPII = 10**((-0.8712d0*Y1 + 8.156d0)*Y1 - 21.833d0)
    else
        SPII = 0.d0
    endif

    end function SPII

!---------------------------------------------------------------------
    double precision function STOTQ1(ABEAM, E_keV, NE19, Te_keV, Zq, Aq)
!----------------------------------------------- Stotq1 [10#-13 cm2]
! Partial neutral beam stopping cross section for
! H(D, T) beam H(D, T, 3He, 4He, Li, Be, B, C, N, O, Ne, Fe, Ar) plasma
!       S Suzuki,  et al "Attenuation og high-energy neutral hydrogen
!       beams in hugh-density plasmas"
! Ne [10#19 m-3]=0.1-100,  Ebeam [KeV]=100-10000,  Te [KeV]=1-50
! E [keV]=Ebeam*Mi/Mb
! Stotq =Zq*(1+Sq*(Zq-1))*S1(E, ne, Te),
! n*Stot =nq*Stotq+... ,  q=1, 2, ..., N for all ion species
!--------------------------------------------- Polevoy A.R. 25.01.09
!---------------------------Ne, Ar added by Leonov V.M. in 1999
      
    double precision, intent(in) :: ABEAM, E_keV, NE19, Te_keV, Zq, Aq

    integer :: i, j, k, iSTOTQ1
    double precision :: A(10), B(3, 2, 2), T, ALT, ALN, ALE, AN, AT, AE, S1, SQ

    data iSTOTQ1/0/
    if (E_keV < 10.d0) write(*, *) 'Illegal use of STOT: Eb<10 keV'
    if (abeam < 1 .or. abeam > 3.) then
        write(*, *) 'NBI (stotq1): no data stopping CS for ABEAM = ',  ABEAM
        STOTQ1 = 0.d0
        return
    endif
 
    if (E_keV < 100.d0) then
        if (abeam == 1.d0) then
            A(1) = -5.29d1
            A(2) = -1.36d0
            A(3) = 7.19d-2
            A(4) = 1.37d-2
            A(5) = 4.54d-1
            A(6) = 4.03d-1
            A(7) = -2.2d-1
            A(8) = 6.660d-2
            A(9) = -6.77d-2
            A(10) = -1.48d-3
        else if (abeam == 2.d0) then
            A(1) = -6.79d1
            A(2) = -1.22d0
            A(3) = 8.14d-2
            A(4) = 1.39d-2
            A(5) = 4.54d-1
            A(6) = 4.65d-1
            A(7) = -2.73d-1
            A(8) = 7.51d-2
            A(9) = -6.3d-2
            A(10) = -5.08-4
        else if (abeam == 3.d0) then
            A(1) = -7.42d1
            A(2) = -1.18d0
            A(3) = 8.43d-2
            A(4) = 1.39d-2
            A(5) = 4.53d-1
            A(6) = 4.91d-1
            A(7) = -2.94d-1
            A(8) = 7.88d-2
            A(9) = -6.12d-2
            A(10) = -1.85d-4
         endif
    else
        if (abeam == 1.d0) then
            A(1) = 1.27d1
            A(2) = 1.25d0
            A(3) = 4.52d-1
            A(4) = 1.05d-2
            A(5) = 5.47d-1
            A(6) = -1.02d-1
            A(7) = 3.6d-1
            A(8) = -2.98d-2
            A(9) = -9.59d-2
            A(10) = 4.21d-3
        else if (abeam == 2.d0) then
            A(1) = 1.41d1
            A(2) = 1.11d0
            A(3) = 4.08d-1
            A(4) = 1.05d-2
            A(5) = 5.47d-1
            A(6) = -4.03d-2
            A(7) = 3.45d-1
            A(8) = -2.88d-2
            A(9) = -9.71d-2
            A(10) = 4.74d-3
        else if (abeam == 3.d0) then
            A(1) = 1.27d1
            A(2) = 1.26d0
            A(3) = 4.49d-1
            A(4) = 1.05d-2
            A(5) = 5.47d-1
            A(6) = -5.77d-3
            A(7) = 3.36d-1
            A(8) = -2.82d-2
            A(9) = -9.74d-2
            A(10) = 4.87d-3
        endif
    endif
 
    if (ZQ <= 1.) then
        B = 0.
    else if (AQ < 6.d0) then  ! He
        if (E_keV < 100.d0) then
            B(1, 1, 1) = -7.92d-1
            B(1, 1, 2) = 4.2d-2
            B(1, 2, 1) = 5.3d-2
            B(1, 2, 2) = -1.39d-2
            B(2, 1, 1) = 3.01d-1
            B(2, 1, 2) = -2.64d-2
            B(2, 2, 1) = -2.99d-2
            B(2, 2, 2) = 6.07d-3
            B(3, 1, 1) = 2.72d-4
            B(3, 1, 2) = 6.11d-3
            B(3, 2, 1) = 3.47d-3
            B(3, 2, 2) = -9.19d-4
        else
            B(1, 1, 1) = 2.31d-1
            B(1, 1, 2) = 3.43d-1
            B(1, 2, 1) = -1.85d-1
            B(1, 2, 2) = -1.62d-2
            B(2, 1, 1) = 1.05d-1
            B(2, 1, 2) = -7.03d-2
            B(2, 2, 1) = 5.31d-2
            B(2, 2, 2) = 3.42d-3
            B(3, 1, 1) = -8.38d-3
            B(3, 1, 2) = 4.15d-3
            B(3, 2, 1) = -3.35d-3
            B(3, 2, 2) = -2.21d-4
        endif
    else if (AQ > 6d0 .and. AQ < 8.d0) then ! Li
        if (E_keV < 100.d0) then
            B(1, 1, 1) = -4.97d-1
            B(1, 1, 2) = 3.38d-2
            B(1, 2, 1) = 3.87d-2
            B(1, 2, 2) = -2.74d-3
            B(2, 1, 1) = 1.55d-1
            B(2, 1, 2) = -2.61d-2
            B(2, 2, 1) = -1.95d-2
            B(2, 2, 2) = 5.07d-4
            B(3, 1, 1) = 1.13d-2
            B(3, 1, 2) = 6.09d-3
            B(3, 2, 1) = 1.63d-3
            B(3, 2, 2) = -1.36d-4
        else
            B(1, 1, 1) = -4.41d-1
            B(1, 1, 2) = 1.29d-1
            B(1, 2, 1) = -1.7d-1
            B(1, 2, 2) = -1.62d-2
            B(2, 1, 1) = 2.77d-1
            B(2, 1, 2) = -1.56d-2
            B(2, 2, 1) = 4.66d-2
            B(2, 2, 2) = 3.79d-3
            B(3, 1, 1) = -1.93d-2
            B(3, 1, 2) = 7.53d-4
            B(3, 2, 1) = -2.86d-3
            B(3, 2, 2) = -2.39d-4
        endif
    else if (AQ > 8.d0 .and. AQ < 10.d0) then ! Be
        if (E_keV < 100.d0) then
            B(1, 1, 1) = 1.12d-1
            B(1, 1, 2) = 4.95d-2
            B(1, 2, 1) = 1.16d-2
            B(1, 2, 2) = -2.86d-3
            B(2, 1, 1) = -1.49d-1
            B(2, 1, 2) = -3.31d-2
            B(2, 2, 1) = -4.26d-3
            B(2, 2, 2) = 9.8d-4
            B(3, 1, 1) = 4.47d-2
            B(3, 1, 2) = 6.52d-3
            B(3, 2, 1) = -3.56d-4
            B(3, 2, 2) = -2.03d-4
        else
            B(1, 1, 1) = -6.13d-1
            B(1, 1, 2) = 5.52d-2
            B(1, 2, 1) = -1.67d-1
            B(1, 2, 2) = -1.59d-2
            B(2, 1, 1) = 3.04d-1
            B(2, 1, 2) = 1.54d-3
            B(2, 2, 1) = 4.36d-2
            B(2, 2, 2) = 3.78d-3
            B(3, 1, 1) = -2.01d-2
            B(3, 1, 2) = -2.16d-4
            B(3, 2, 1) = -2.51d-3
            B(3, 2, 2) = -2.27d-4
        endif
    else if (AQ > 10.d0 .and. AQ < 12.d0) then ! B
        if (E_keV < 100.d0) then
            B(1, 1, 1) = 1.22d-1
            B(1, 1, 2) = 5.27d-2
            B(1, 2, 1) = -4.3d-4
            B(1, 2, 2) = -3.18d-3
            B(2, 1, 1) = -1.51d-1
            B(2, 1, 2) = -3.64d-2
            B(2, 2, 1) = 3.43d-3
            B(2, 2, 2) = 1.51d-3
            B(3, 1, 1) = 4.2d-2
            B(3, 1, 2) = 6.92d-3
            B(3, 2, 1) = -1.41d-3
            B(3, 2, 2) = -2.9d-4
        else
            B(1, 1, 1) = -7.32d-1
            B(1, 1, 2) = 1.83d-2
            B(1, 2, 1) = -1.55d-1
            B(1, 2, 2) = -1.72d-2
            B(2, 1, 1) = 3.21d-1
            B(2, 1, 2) = 9.46d-3
            B(2, 2, 1) = 3.97d-2
            B(2, 2, 2) = 4.2d-3
            B(3, 1, 1) = -2.04d-2
            B(3, 1, 2) = -6.19d-4
            B(3, 2, 1) = -2.24d-3
            B(3, 2, 2) = -2.54d-4
        endif
    else if (AQ > 11.d0 .and. AQ < 13.d0) then ! C
        if (E_keV < 100.d0) then
            B(1, 1, 1) = 1.61d-1
            B(1, 1, 2) = 5.98d-2
            B(1, 2, 1) = -3.36d-3
            B(1, 2, 2) = -4.26d-3
            B(2, 1, 1) = -1.57d-1
            B(2, 1, 2) = -3.96d-2
            B(2, 2, 1) = 4.6d-3
            B(2, 2, 2) = 2.19d-3
            B(3, 1, 1) = 3.91d-2
            B(3, 1, 2) = 7.11d-3
            B(3, 2, 1) = -1.44d-3
            B(3, 2, 2) = -3.85d-4
        else
            B(1, 1, 1) = -1.01d0
            B(1, 1, 2) = -8.65d-3
            B(1, 2, 1) = -1.24d-1
            B(1, 2, 2) = -1.45d-2
            B(2, 1, 1) = 3.91d-1
            B(2, 1, 2) = 1.61d-2
            B(2, 2, 1) = 2.98d-2
            B(2, 2, 2) = 3.32d-3
            B(3, 1, 1) = -2.48d-2
            B(3, 1, 2) = -1.04d-3
            B(3, 2, 1) = -1.52d-3
            B(3, 2, 2) = -1.89d-4
        endif
    else if (AQ > 12.d0 .and. AQ < 15.d0) then ! N
        if (E_keV < 100.d0) then
            B(1, 1, 1) = 1.34d-1
            B(1, 1, 2) = 5.24d-2
            B(1, 2, 1) = -4.69d-3
            B(1, 2, 2) = -3.06d-3
            B(2, 1, 1) = -1.3d-1
            B(2, 1, 2) = -3.52d-2
            B(2, 2, 1) = 5.31d-3
            B(2, 2, 2) = 1.64d-3
            B(3, 1, 1) = 3.31d-2
            B(3, 1, 2) = 6.35d-3
            B(3, 2, 1) = -1.51d-3
            B(3, 2, 2) = -3.04d-4
        else
            B(1, 1, 1) = -1.d0
            B(1, 1, 2) = -4.15d-2
            B(1, 2, 1) = -9.76d-2
            B(1, 2, 2) = -1.1d-2
            B(2, 1, 1) = 3.7d-1
            B(2, 1, 2) = 2.28d-2
            B(2, 2, 1) = 2.15d-2
            B(2, 2, 2) = 2.33d-3
            B(3, 1, 1) = -2.22d-2
            B(3, 1, 2) = -1.33d-3
            B(3, 2, 1) = -9.25d-4
            B(3, 2, 2) = -1.15d-4
        endif
    else if (AQ > 15.d0 .and. AQ < 17.d0) then !O
        if (E_keV < 100.) then
            B(1, 1, 1) = 1.07d-1
            B(1, 1, 2) = 4.31d-2
            B(1, 2, 1) = -2.83d-3
            B(1, 2, 2) = -1.82d-3
            B(2, 1, 1) = -1.06d-1
            B(2, 1, 2) = -2.89d-2
            B(2, 2, 1) = 3.87d-3
            B(2, 2, 2) = 9.1d-4
            B(3, 1, 1) = 2.77d-2
            B(3, 1, 2) = 5.27d-3
            B(3, 2, 1) = -1.23d-3
            B(3, 2, 2) = -1.89d-4
        else
            B(1, 1, 1) = -9.89d-1
            B(1, 1, 2) = -4.98d-2
            B(1, 2, 1) = -6.36d-2
            B(1, 2, 2) = -7.75d-3
            B(2, 1, 1) = 3.52d-1
            B(2, 1, 2) = 2.34d-2
            B(2, 2, 1) = 1.17d-2
            B(2, 2, 2) = 1.42d-3
            B(3, 1, 1) = -2.04d-2
            B(3, 1, 2) = -1.29d-3
            B(3, 2, 1) = -2.75d-4
            B(3, 2, 2) = -5.45d-5
        endif
    else if (AQ > 19.d0 .and. AQ < 21.d0) then !Ne
        B(1, 1, 1) = -1.32 !Ne by Leonov Ekev > 100
        B(1, 1, 2) = 1.8d-2
        B(1, 2, 1) = -9.7d-2
        B(1, 2, 2) = -1.13d-2
        B(2, 1, 1) = 4.45d-1
        B(2, 1, 2) = -4.5d-3
        B(2, 2, 1) = 2.3d-2
        B(2, 2, 2) = 2.85d-3
        B(3, 1, 1) = -2.78d-2
        B(3, 1, 2) = 2.05d-3
        B(3, 2, 1) = -1.4d-3
        B(3, 2, 2) = -1.7d-4
    else if (AQ > 39.d0 .and. AQ < 41.d0) then !Ar
        B(1, 1, 1) = -1.12 !Ar by Leonov Ekev > 100
        B(1, 1, 2) = 7.3d-2
        B(1, 2, 1) = -7.0d-3
        B(1, 2, 2) = -6.7d-3
        B(2, 1, 1) = 3.62d-1
        B(2, 1, 2) = -2.4d-2
        B(2, 2, 1) = 1.6d-2
        B(2, 2, 2) = 1.68d-3
        B(3, 1, 1) = -2.18d-2
        B(3, 1, 2) = 2.92d-3
        B(3, 2, 1) = -9.5d-4
        B(3, 2, 2) = -9.5d-5
    else  ! Fe and others
        if (iSTOTQ1 < 1) then
            write(*, *) 'NBI: No stop cross section DATA for Aimp=', AQ
            write(*, *) '    cross section for Aimp=56 is substituted instead'
            iSTOTQ1 = 1
        endif
        if (E_keV < 100.) then
            B(1, 1, 1) = -4.65d-5
            B(1, 1, 2) = -7.29d-4
            B(1, 2, 1) = -3.1d-3
            B(1, 2, 2) = -5.17d-4
            B(2, 1, 1) = -1.34d-2
            B(2, 1, 2) = -5.06d-5
            B(2, 2, 1) = 2.87d-3
            B(2, 2, 2) = 2.74d-4
            B(3, 1, 1) = 5.04d-3
            B(3, 1, 2) = 2.13d-4
            B(3, 2, 1) = -6.65d-4
            B(3, 2, 2) = -5.83d-5
        else
            B(1, 1, 1) = -1.93d-1
            B(1, 1, 2) = 5.03d-3
            B(1, 2, 1) = 1.06d-1
            B(1, 2, 2) = 1.25d-2
            B(2, 1, 1) = 2.36d-3
            B(2, 1, 2) = -7.41d-3
            B(2, 2, 1) = -3.88d-2
            B(2, 2, 2) = -4.6d-3
            B(3, 1, 1) = 9.39d-3
            B(3, 1, 2) = 1.58d-3
            B(3, 2, 1) = 3.32d-3
            B(3, 2, 2) = 3.88d-4
        endif
    endif

    SQ = 0.d0
    T = max(Te_keV, 1.d0)
    ALT = LOG(T)
    ALN = LOG(NE19)
    ALE = LOG(E_keV)
    S1 = 1.d-3*A(1)*(1.d0 + (A(2) + A(3)*ALE)*ALE) * (1.d0 + (1.d0 - exp(-A(4)*NE19))**A(5) * &
        (A(6) + (A(7) + A(8)*ALE)*ALE)) * (1.d0 + (A(9) + A(10)*ALT)*ALT)/E_keV
    AT = 1.0d0
    do K=1, 2
        AN = AT
        do J=1, 2
            AE = AN
            do I=1, 2
                SQ = SQ + B(I, J, K)*AE
                AE = AE*ALE
            enddo
            SQ = SQ + B(3, J, K)*AE
            AN = AN*ALN
        enddo
        AT = AT*ALT
    enddo
         
    STOTQ1 = S1*ZQ*(1.d0 + (ZQ - 1.0d0)*SQ)

    end function STOTQ1

!---------------------------------------------------------------------
    double precision function sv_reac(A_main, E_NBI_keV, A_NBI, n_e, &
        Te_keV, Ti_keV, yAi, calc_fus, Acoeff, Bcoeff)
!---------------------------------------------------------------------
! <SigmaV dt> probability for d with EBEAM keV to burn out on maxwellian d
! with Ti  during the slowing down to Ti in reaction:
!---------------------------------------------------------------------
! Crossection by H-S. Bosch,  G.M. Hale
!  NF,  V 32 ,  N 4,  (1992) p 611-631
! Corrected for finit ion temperature according to:
!  D.R.Mikkelsen,  NF V 29,  N 7,  (1989) p 1113-1115
! + second derivative is addedd S' -> S'+ 2 E S''(16-JUL-13)
!
! Use:
!    Sdt245[10^-19/m^3/s] =
!   =Pbeam[MW]*svddnb2*625/EBEAM[keV]*Ndeut[10^19m-3]
!
! input: E_NBI_keV[energy, keV], A_NBI[mass,  a.u.],
!  n_e[Ne, 10^19m-3], Te_keV[Te, keV], yAi[amain,  mass,  a.u.]
! ....Logarithm e
!  YLE =15.85+LOG(TE(j)/sqrt(NE(j)))
!  TauES = 2.d0*A_NBI*TE(j)^1.5/LnE/NE(j)
!  V = Sqrt(2 T/M)
!---------------------------------------------------------------------

    double precision, parameter :: YTMIN=0.01d0

    integer, intent(in) :: calc_fus
    double precision, intent(in) :: A_main, E_NBI_keV, A_NBI, n_e, &
        Te_keV, Ti_keV, yAi, Acoeff(5), Bcoeff(4)
  
    integer :: jk, jend
    double precision :: YX3, YX2, YX, YE, YECM, YSQ, YASS, YBSS, YSS2, YSIG, &
        YXC3, YECDEB, YDS, YLE, YLI, YBG, &
        YA0, YB0, YAS, YBS, YMU, YGAM, YBET,  &
        YB, Vth2, MVth2, MVth24, YVs, YSS, YR, YD, Y27, VtdVb2, YVB, &
        Y13, YEMIN, YRMD3, YRPD3, YVb0, YCOEF

    if (E_NBI_keV <= 0.d0) then
        sv_reac = 0.d0
        return
    endif

    yMU  = A_NBI*A_main/(A_NBI + A_main)
    YBG  = 31.397d0*sqrt(yMU) ! B =44.4d0
    YECM = yMu/A_NBI*E_NBI_keV
    YSQ  = sqrt(YECM)

    if (Ti_keV >= YTMIN) then
        y13    = 1.d0/3.d0
        Y27    = 1.d0/27.d0
        Vth2   = 2.d0*Ti_keV/A_main
        MVth2  = yMU*Vth2
        MVth24 = MVth2/4.d0
        YVB    = sqrt(2.d0*E_NBI_keV/A_NBI)
        VtdVb2 = Vth2/YVB**2
        YB     = 22.2d0*VtdVb2/YVb ! bet =YB/Vb^3
        YEMIN  = Ti_keV
    else
        YEMIN = YTMIN
    endif

! Coulomb Log Le, Li
    YLE = (15.85d0 + LOG(Te_keV/sqrt(n_e)))
    if (E_NBI_keV > 100.d0*A_NBI) then
        YLI = 23.7d0+log(yAi/(yAi+A_NBI) * sqrt(1.d-3*A_NBI*E_NBI_keV*Te_keV/n_e))
    else
        YLI = 25.4d0+log(1.d-3*E_NBI_keV*yAi/(yAi + A_NBI) * sqrt(Te_keV/n_e))
    endif

    jend = 1000
    YDS = (1.d0/jend)
    sv_reac = 0.d0
    YECDEB = 14.6d0*Te_keV*A_NBI/E_NBI_keV/(YLE*yAi/YLI)**0.667
    YXC3 = YECDEB*sqrt(YECDEB)
    do JK=1, jend
        YX  = (YDS*JK)
        YX2 = YX**2
        YX3 = YX2*YX
        YE  = YECM*YX2
        if (E_NBI_keV*YX2 > YEMIN) then
            if (Ti_keV >= YTMIN .and. calc_fus < 2) then
                yBET  = YB/YX3
                YR    = yBET/2.d0+Y27
                YD    = sqrt(YBET*(yBET/4.d0+Y27))
                YRMD3 = (YR - YD)**Y13
                YRPD3 = (YR + YD)**Y13
                YVb0  = Y13 + YRMD3 + YRPD3
                YVS   = YX*YVb0
                YGAM  = 3.0d0 - 2.0/YVb0
                YE    = YECM*YVS**2
                YA0   = Acoeff(1) + YE*(Acoeff(2) + YE*(Acoeff(3) + YE*(Acoeff(4) + YE*Acoeff(5))))
                YB0   = 1.d0 + YE*(Bcoeff(1) + YE*(Bcoeff(2) + YE*(Bcoeff(3) + YE*Bcoeff(4))))
                YAS   = Acoeff(2) + YE*(2.d0*Acoeff(3)+YE*(3.d0*Acoeff(4) + 4.d0*YE*Acoeff(5)))
                YBS   = Bcoeff(1) + YE*(2.d0*Bcoeff(2)+YE*(3.d0*Bcoeff(3) + 4.d0*YE*Bcoeff(4)))
                YASS  = 2.d0*(Acoeff(3) + 3.d0*YE*(Acoeff(4) + 2.d0*YE*Acoeff(5)))
                YBSS  = 2.d0*(Bcoeff(2) + 3.d0*YE*(Bcoeff(3) + 2.d0*YE*Bcoeff(4)))
                YSS   = YAS/YA0 - YBS/YB0
                YSS2  = YSS*(1.d0 - 4.d0*YE*YBS/YB0) + 2.d0*YE*(YASS/YA0 - YBSS/YB0)

                Ycoef = exp(-YX2*(YVb0 - 1.d0)**2/VtdVb2)*yVb0/sqrt(YGAM)*(1.d0 + MVth24*YSS2/YGAM + &
                    1.5d0*(1.d0-1.d0/yVb0)/YGAM**2*(MVth2*YSS-VtdVb2/YVS**2))
                YSIG = exp(-YBG/YVS/YSQ)/YE*YA0/YB0*yVb0*Ycoef
            else
                YVS = YX
                YE = YECM*YX2
                YA0 = Acoeff(1) + YE*(Acoeff(2) + YE*(Acoeff(3) + YE*(Acoeff(4) + YE*Acoeff(5))))
                YB0 = 1.d0 + YE*(Bcoeff(1) + YE*(Bcoeff(2) + YE*(Bcoeff(3) + YE*Bcoeff(4))))
                YSIG = exp(-YBG/YVS/YSQ)/YE*YA0/YB0
            endif
            sv_reac = sv_reac + YSIG/(1. + YXC3/YX3)
        endif
    enddo

    sv_reac = sv_reac*YDS*1.d-3*4.38d-4*sqrt(E_NBI_keV/A_NBI)*2.d0*A_NBI/n_e*Te_keV*sqrt(Te_keV)/YLE
 
    end function sv_reac
    
!---------------------------------------------------------------------
    double precision function svddnp1(E_NBI_keV, A_NBI, n_e, &
        Te_keV, Ti_keV, yAi, calc_fus)

! d(Ebeam) + d(Ti) -> He3(870 keV) + n (2450 keV)

    double precision, parameter :: A_main=2.d0, &
        Acoeff(5) = (/5.3701d4, 3.3027d2, -0.12706d0, 2.9327d-5, -2.5151d-9/), &
        Bcoeff(4) = (/0.d0, 0.d0, 0.d0, 0.d0/)

    integer, intent(in) :: calc_fus
    double precision, intent(in) :: E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi

    if (A_NBI /= 2.d0) then
        svddnp1 = 0.d0
    else
        svddnp1 = sv_reac(A_main, E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi, calc_fus, Acoeff, Bcoeff)
    endif

    end function svddnp1

!---------------------------------------------------------------------
    double precision function svddnp2(E_NBI_keV, A_NBI, n_e, &
        Te_keV, Ti_keV, yAi, calc_fus)

! d(Ebeam) + d(Ti) -> t(1008 keV) + p(3025 keV)

    double precision, parameter :: A_main=2.d0, &
        Acoeff(5) = (/5.5576d4, 2.1054d2, -3.2638d-2, 1.4987d-6, 1.1881d-10/), &
        Bcoeff(4) = (/0.d0, 0.d0, 0.d0, 0.d0/)

    integer, intent(in) :: calc_fus
    double precision, intent(in) :: E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi

    if (A_NBI /= 2.d0) then
        svddnp2 = 0.d0
    else
        svddnp2 = sv_reac(A_main, E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi, calc_fus, Acoeff, Bcoeff)
    endif

    end function svddnp2
  
!---------------------------------------------------------------------
    double precision function svdtbp(E_NBI_keV, A_NBI, n_e, &
        Te_keV, Ti_keV, yAi, calc_fus)

! d(t)(Ebeam) + t(d)(Ti) -> He4(3524 keV) + n(14072 keV)

    double precision, parameter :: Acoeff(5) = (/6.927d4, 7.454d8, 2.05d6, 5.2002d4, 0.d0/), &
        Bcoeff(4) = (/63.8d0, -0.995d0, 6.981d-5, 1.728d-4/)

    integer, intent(in) :: calc_fus
    double precision, intent(in) :: E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi
    double precision :: A_main

    if (A_NBI /= 2.d0 .and. A_NBI /= 3.d0) then
        svdtbp = 0.d0
    else
        if (A_NBI == 3.d0) then
            A_main = 2.d0 ! t NBI in d bulk plasma
        else
            A_main = 3.d0 ! d NBI in t bulk plasma
        endif
        svdtbp = sv_reac(A_main, E_NBI_keV, A_NBI, n_e, Te_keV, Ti_keV, yAi, calc_fus, Acoeff, Bcoeff)
    endif

    end function svdtbp

end module cross_sections
