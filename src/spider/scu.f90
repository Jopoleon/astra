subroutine lsq_sur6(r, z, u, n, c, rm, zm, um, dp)

!...    list square surf. fitting
!       c(1)+c(2)*(r-r1)+c(3)*(z-z1)+c(4)*(r-r1)**2+c(5)*(z-z1)**2+c(6)*(r-r1)*(z-z1)

implicit none

real*8, intent(in), dimension(*) :: r, z, u
real*8, intent(out) :: rm, zm
real*8, intent(out), dimension(*) :: dp

integer :: k, n
integer, dimension(6) :: IP
real*8 :: s_r, s_r2, s_r3, s_r4, s_z, s_z2, s_z3, s_z4, s_rz, s_rz2, s_r2z, &
    s_r2z2, s_r3z, s_rz3, s_u, s_ur, s_ur2, s_uz, s_uz2, s_urz, &
    det, det_r, det_z, um, dudr, dudr2, dudz, dudz2, dudrdz
real*8, dimension(6) :: B, C
real*8, dimension(6, 6) :: A

s_r   = 0.d0
s_z   = 0.d0
s_r2  = 0.d0
s_z2  = 0.d0
s_rz  = 0.d0
s_r3  = 0.d0
s_rz2 = 0.d0
s_r2z = 0.d0
s_z3  = 0.d0
s_r4  = 0.d0
s_r2z2= 0.d0
s_r3z = 0.d0
s_z4  = 0.d0
s_rz3 = 0.d0
s_u   = 0.d0
s_ur  = 0.d0
s_uz  = 0.d0
s_ur2 = 0.d0
s_uz2 = 0.d0
s_urz = 0.d0

do k=1, n
   s_r   = s_r + (r(k) - r(1))
   s_z   = s_z + (z(k) - z(1))
   s_r2  = s_r2 + (r(k) - r(1))**2
   s_z2  = s_z2 + (z(k) - z(1))**2
   s_rz  = s_rz + (r(k) - r(1))*(z(k) - z(1))
   s_u   = s_u + u(k)
   s_r3  = s_r3  + (r(k) - r(1))**3
   s_rz2 = s_rz2 + (r(k) - r(1))*(z(k) - z(1))**2
   s_r2z = s_r2z + (r(k) - r(1))**2*(z(k) - z(1))
   s_ur  = s_ur  + u(k)*(r(k) - r(1))
   s_z3  = s_z3 + (z(k) - z(1))**3
   s_uz  = s_uz + u(k)*(z(k) - z(1))
   s_r4  = s_r4 + (r(k) - r(1))**4
   s_r2z2= s_r2z2 + (r(k) - r(1))**2*(z(k) - z(1))**2
   s_r3z = s_r3z  + (r(k) - r(1))**3*(z(k) - z(1))
   s_ur2 = s_ur2  +  u(k)*(r(k)-r(1))**2
   s_z4  = s_z4  + (z(k) - z(1))**4
   s_rz3 = s_rz3 + (r(k) - r(1))*(z(k) - z(1))**3
   s_uz2 = s_uz2 + u(k)*(z(k) - z(1))**2
   s_urz = s_urz + u(k)*(r(k) - r(1))*(z(k) - z(1))
enddo

! THE CREATION THE MATRIX A and RIGHT HAND B:

A(1,1) = FLOAT(n)
A(1,2) = s_r
A(1,3) = s_z
A(1,4) = s_r2
A(1,5) = s_z2
A(1,6) = s_rz

B(1) = s_u

A(2,1) = s_r
A(2,2) = s_r2
A(2,3) = s_rz
A(2,4) = s_r3
A(2,5) = s_rz2
A(2,6) = s_r2z

B(2) = s_ur

A(3,1) = s_z
A(3,2) = s_rz
A(3,3) = s_z2
A(3,4) = s_r2z
A(3,5) = s_z3
A(3,6) = s_rz2

B(3) = s_uz

A(4,1) = s_r2
A(4,2) = s_r3
A(4,3) = s_r2z
A(4,4) = s_r4
A(4,5) = s_r2z2
A(4,6) = s_r3z

B(4) = s_ur2

A(5,1) = s_z2
A(5,2) = s_rz2
A(5,3) = s_z3
A(5,4) = s_r2z2
A(5,5) = s_z4
A(5,6) = s_rz3

B(5) = s_uz2

A(6,1) = s_rz
A(6,2) = s_r2z
A(6,3) = s_rz2
A(6,4) = s_r3z
A(6,5) = s_rz3
A(6,6) = s_r2z2

B(6) = s_urz

!  THE SOLUTION THE MATRIX EQUATION Ac=B.

CALL GE(6, 6, A, B, c, IP)

!   magnetic axis

det = 4.d0*c(4)*c(5) - c(6)**2

det_r = -2.d0*c(2)*c(5) + c(6)*c(3)
det_z = -2.d0*c(3)*c(4) + c(6)*c(2)
rm = det_r/det
zm = det_z/det
um = c(1) + c(2)*rm + c(3)*zm + c(4)*rm**2 + c(5)*zm**2 + c(6)*rm*zm
rm = rm + r(1)
zm = zm + z(1)

dudr   = c(2)
dudz   = c(3)
dudr2  = 2.d0*c(4)
dudz2  = 2.d0*c(5)
dudrdz = c(6)

dp(1) = dudr
dp(2) = dudz
dp(3) = dudr2
dp(4) = dudrdz
dp(5) = dudz2

return
end subroutine lsq_sur6

!---------------------------------------------------------------------
subroutine DERIV5(X, Y, F, M, N, U)

!...    DEFINITON OF THE FIRST AND SECOND DERIVATIONS.

implicit none

integer, intent(in) :: M, N
real*8, intent(in), dimension(*) :: X, Y, F
real*8, intent(out), dimension(*) :: U

integer :: i, j, ifail
integer, dimension(5) :: IP
real*8 :: sdx2, sdx3, sdx4, sdy2, sdy3, sdy4, sdxdy, SDX2DY, SDXDY2, SDX3DY, &
    SDX2Y2, SDXDY3, FDX, FDX2, FDY, FDY2, FDXDY, DX, DX2, DY, DY2, DF, &
    DAVER, DAVER2
real*8, dimension(5) :: B, S
real*8, dimension(5, 5) :: A

SDX2   = 0.d0
SDY2   = 0.d0
SDXDY  = 0.d0
SDX3   = 0.d0
SDX2DY = 0.d0
SDXDY2 = 0.d0
SDY3   = 0.d0
SDX4   = 0.d0
SDX3DY = 0.d0
SDX2Y2 = 0.d0
SDXDY3 = 0.d0
SDY4   = 0.d0
DAVER2 = 0.d0
FDX    = 0.d0
FDY    = 0.d0
FDX2   = 0.d0
FDXDY  = 0.d0
FDY2   = 0.d0

do I=2, M
    DX  = X(I) - X(1)
    DY  = Y(I) - Y(1)
    DF  = F(I) - F(1)
    DX2 = DX**2
    DY2 = DY**2
    SDX2   = SDX2   + DX2
    SDY2   = SDY2   + DY2
    SDXDY  = SDXDY  + DX *DY
    SDX3   = SDX3   + DX2*DX
    SDX2DY = SDX2DY + DX2*DY
    SDXDY2 = SDXDY2 + DX *DY2
    SDY3   = SDY3   + DY2*DY
    SDX4   = SDX4   + DX2*DX2
    SDX3DY = SDX3DY + DX2*DX  *DY
    SDX2Y2 = SDX2Y2 + DX2*DY2
    SDXDY3 = SDXDY3 + DX *DY2 *DY
    SDY4   = SDY4   + DY2*DY2
    DAVER2 = DAVER2 + DX2 + DY2
    FDX   = FDX   + DF*DX
    FDY   = FDY   + DF*DY
    FDX2  = FDX2  + DF*DX2
    FDXDY = FDXDY + DF*DX*DY
    FDY2  = FDY2  + DF*DY2
enddo

! Creating matrix A:

A(1, 1) = SDX2
A(1, 2) = SDXDY
A(1, 3) = SDX3  *0.5d0
A(1, 4) = SDX2DY
A(1, 5) = SDXDY2*0.5d0
A(2, 1) = SDXDY
A(2, 2) = SDY2
A(2, 3) = SDX2DY*0.5d0
A(2, 4) = SDXDY2
A(2, 5) = SDY3  *0.5d0
A(3, 1) = SDX3
A(3 ,2) = SDX2DY
A(3, 3) = SDX4  *0.5d0
A(3, 4) = SDX3DY
A(3, 5) = SDX2Y2*0.5d0
A(4, 1) = SDX2DY
A(4, 2) = SDXDY2
A(4, 3) = SDX3DY*0.5d0
A(4, 4) = SDX2Y2
A(4, 5) = SDXDY3*0.5d0
A(5, 1) = SDXDY2
A(5, 2) = SDY3
A(5, 3) = SDX2Y2*0.5d0
A(5, 4) = SDXDY3
A(5, 5) = SDY4  *0.5d0

! Creating right hand B:

B(1) = FDX
B(2) = FDY
B(3) = FDX2
B(4) = FDXDY
B(5) = FDY2

!  THE WEIGHTS S(I):

DAVER = SQRT(DAVER2)

S(1) = 1.d0/SQRT(SDX2)/DAVER
S(2) = 1.d0/SQRT(SDY2)/DAVER
S(3) = 1.d0/SDX2/DAVER
S(4) = 1.d0/SQRT(SDX2*SDY2)/DAVER
S(5) = 1.d0/SDY2/DAVER

do I=1, N
    do J=1, N
        A(I, J) = S(I)*A(I, J)
    enddo
    B(I) = S(I)*B(I)
enddo

!  THE SOLUTION THE MATRIX EQUATION AU=B.

IFAIL = 0

CALL GE(5, 5, A, B, U, IP)

return
end subroutine DERIV5

!---------------------------------------------------------------------
SUBROUTINE DERIV9(X, Y, F, M, N, U)
! Cubic fit

implicit none

integer, intent(in) :: M, N
real*8, intent(in), dimension(M) :: X, Y, F
real*8, intent(out), dimension(N) :: U

integer :: i, k, l, ifail
integer, dimension(9) :: IP
real*8 :: dx, dy, df, dxx, dxy, dyy, dxxx, dxxy, dxyy, dyyy, akl, bk
real*8, dimension(9) :: B
real*8, dimension(9, 9) :: A
real*8, dimension(M) :: BB
real*8, dimension(M, 9) :: AA

do I=1, M-1
    DX = X(I+1) - X(1)
    DY = Y(I+1) - Y(1)
    DF = F(I+1) - F(1)
    DXX = 0.5d0*DX**2
    DXY = DX*DY
    DYY = 0.5d0*DY**2
    AA(I, 1) = DX
    AA(I, 2) = DY
    AA(I, 3) = DXX
    AA(I, 4) = DXY
    AA(I, 5) = DYY
! cubic
    DXXX = DX**3
    DXXY = DX**2 * DY
    DXYY = DX * DY**2
    DYYY = DY**3
    AA(I, 6) = DXXX
    AA(I, 7) = DXXY
    AA(I, 8) = DXYY
    AA(I, 9) = DYYY
    BB(I) = DF
enddo

! Creating matrix A:

do K=1, 9
    do L=1, 9
        AKL = 0.d0
        do I=1, M-1
            AKL = AKL + AA(I, K)*AA(I, L)
        enddo
        A(K, L) = AKL
   enddo
enddo

!...     THE CREATION THE RIGHT HAND B:

do K=1, 9
    BK = 0.d0
    do I=1, M-1
        BK = BK + AA(I, K)*BB(I)
    enddo
    B(K) = BK
enddo

!...     THE SOLUTION THE MATRIX EQUATION AU=B.

IFAIL = 0

CALL GE(N, 9, A, B, U, IP)

IF (IFAIL /= 0) WRITE (*, *) ' DERIV9:IFAIL= ', IFAIL

return
end subroutine DERIV9

!---------------------------------------------------------------------
subroutine GE(N, NZ, A, X, Y, IP)
!...    GAUSS ELIMINATION

implicit none

integer, intent(in) :: N, NZ
real*8, intent(inout), dimension(N) :: X
integer, intent(out), dimension(N) :: IP
real*8, intent(out), dimension(NZ, NZ) :: A
real*8, intent(out), dimension(N) :: Y

integer :: i, j, k, jm, ipe
real*8 :: rm, am, aba, ad, ape

do I=1, N
    IP(I) = I
enddo

do I=1, N-1
!...      MAX ROW ELEMENT SEARCH
    RM = ABS(A(I, I))
    JM = I
    do J=I, N
        ABA = ABS(A(I, J))
        IF (ABA > RM) THEN
            RM = ABA
            JM = J
        ENDIF
    enddo
!...      PERMUTATIONS
    IPE = IP(I)
    IP(I) = IP(JM)
    IP(JM) = IPE
    do K=1, N
        APE = A(K, I)
        A(K, I) = A(K, JM)
        A(K, JM) = APE
    enddo
    AD = 1.d0/A(I, I)
    do K=I+1, N
        AM = A(K, I)*AD
        X(K) = X(K) - AM*X(I)
        do J=I, N
            A(K, J) = A(K, J) - AM*A(I, J)
        enddo
    enddo
enddo

Y(N) = X(N)/A(N, N)
do I=N-1, 1, -1
    AD = 1./A(I, I)
    Y(I) = X(I)
    do J=N, I+1, -1
        Y(I) = Y(I) - A(I, J)*Y(J)
    enddo
    Y(I) = Y(I)*AD
enddo

! BACK PERMUTATION
do I=1, N
    X(IP(I)) = Y(I)
enddo
do I=1, N
    Y(I) = X(I)
enddo

return
end subroutine GE

!---------------------------------------------------------------------
REAL*8 FUNCTION GREENI(R, Z, RP, ZP)
!  GREEN'S FUNCTION FOR TOROIDAL CURRENT LOOP IN INFINITY AREA.

implicit none

real*8, intent(in) :: R, Z, RP, ZP

integer :: IFAILK, IFAILE
real*8 :: T, TT, ELCK, ELCE, S21BBF, S21BCF

T = SQRT(4.d0*R*RP/((R + RP)**2 + (Z - ZP)**2))

! CALCULATION OF COMPLETE ELLIPTIC INTEGRALS OF 1TH AND 2TH KIND.

TT = (1.d0 - T**2)

! THE FUNCTIONS MMDELK, MMDELE FROM IMSL LIBRARY.
!
!        ELCK=MMDELK(2,T,IFAILK)
!        ELCE=MMDELE(2,T,IFAILE)
!
! THE FUNCTIONS S21BBF, S21BCF FROM NAG LIBRARY.

ELCK =                  S21BBF(0.D0, TT, 1.D0, IFAILK)
ELCE = ELCK - T**2/3.D0*S21BCF(0.D0, TT, 1.D0, IFAILE)

IF (IFAILK /= 0.OR.IFAILE /= 0) WRITE(*,*)' ATTENTION! GREEN F.:', &
          '  IFAILK, IFAILE= ', IFAILK, IFAILE

GREENI = ( (1.D0 - T**2 * 0.5D0)*ELCK - ELCE )*( SQRT(R*RP)/T )

return
end function GREENI

!---------------------------------------------------------------------
subroutine grGREN(R0, Z0, R, Z, dGdr, dGdz)

! green's function derivatives at the point (R,Z)
! the source position is (R0,Z0)
! to obtain the real values of derivatives dGdr,dGdz must be devided by pi

implicit none

real*8, intent(in) :: R0, Z0, R, Z
real*8, intent(out) :: dGdr, dGdz

integer :: IFAILK, IFAILE
real*8 :: t, tt, fK, fE, Q, dtdr, dtdz, dQdt, S21BBF, S21BCF

t = SQRT( 4.D0*R*R0/( (R + R0)**2 + (Z - Z0)**2 ) )
tt = (1.D0 - t**2)

! THE FUNCTIONS S21BBF, S21BCF FROM NAG LIBRARY.

fK =               S21BBF(0.D0, tt, 1.D0, IFAILK)
fE = fK -t**2/3.D0*S21BCF(0.D0, tt, 1.D0, IFAILE)

IF (IFAILK /= 0.OR.IFAILE /= 0) WRITE(*, *) ' ATTENTION! GREEN F.:', &
         '  IFAILK, IFAILE= ', IFAILK, IFAILE

Q = ((1.D0 - t**2 * 0.5D0)*fK - fE)/t
dtdr = 0.5d0*(t/R - t**3*(R+R0)/(2.d0*R*R0))
dtdz = -t**3 * (z - z0)/(4.d0*R*R0)
dQdt = -Q/t + 0.5d0*(fE/tt - fK)
dGdr = 0.5d0*SQRT(R0/R)*Q + SQRT(R0*R)*dQdt*dtdr
dGdz = SQRT(R0*R)*dQdt*dtdz

return
end subroutine GRGREN

!---------------------------------------------------------------------
real*8 function hypoten(x1, y1, x2, y2)

implicit none

real*8, intent(in) :: x1, x2, y1, y2

hypoten = SQRT((x1 - x2)**2 + (y1 - y2)**2)

return
end function hypoten

!---------------------------------------------------------------------
real*8 function tint(t, pscal, pvec, eps)

implicit none

real*8, intent(in) :: t, pscal, pvec, eps

tint = (T - PSCAL)*(0.5d0*DLOG(PVEC**2 + (T - PSCAL)**2 + EPS**2) - 1.d0) + &
    PVEC*ATAN((T - PSCAL)/(PVEC + EPS))

return
end function tint

!---------------------------------------------------------------------
subroutine bint(X, Y, R0, Z0, r1, z1, F, I)

use sp_parameters, only: pi

implicit none

integer, intent(in) :: i
real*8, intent(in) :: X, Y, R0, Z0, r1, z1
real*8, intent(out) :: F

real*8 :: epsh, W, H, eps, rm, zm, dtm, dfm, pscal, pvec, dfj, dfj1, dtj, dtj1, E
real*8 :: tint, hypoten, greeni

EPSH = 1.d-9
F = 0.d0
W = 0.d0
H = hypoten(R0,Z0,R1,Z1)
EPS = EPSH*H

!  INTEGRAL ALONG THE EDGE OF REGION

RM = 0.5d0*(R0 + R1)
ZM = 0.5d0*(Z0 + Z1)
DTM = hypoten( X, Y, RM, ZM)
DFM = hypoten(-X, Y, RM, ZM)

IF (I /= 0) THEN
    PSCAL = ((-X - R0)*(R1 - R0) + (Y - Z0)*(Z1 - Z0))/H
    PVEC  = ((-X - R0)*(Z1 - Z0) - (Y - Z0)*(R1 - R0))/H
    W = TINT(H, pscal, pvec, eps) - TINT(0.D0, pscal, pvec, eps)
    PSCAL = ((X - R0)*(R1 - R0) + (Y - Z0)*(Z1 - Z0))/H
    PVEC  = ((X - R0)*(Z1 - Z0) - (Y - Z0)*(R1 - R0))/H
    W = W - (TINT(H, pscal, pvec, eps) - TINT(0.D0, pscal, pvec, eps))
    IF (DTM < 0.25d0*H) THEN
        DFJ  = hypoten(-X, Y, R0, Z0)
        DFJ1 = hypoten(-X, Y, R1, Z1)
        DTJ  = hypoten( X, Y, R0, Z0)
        DTJ1 = hypoten( X, Y, R1, Z1)
        W = W*DFM - 0.5d0*H*(DLOG(DFJ/DTJ)*DFJ + DLOG(DFJ1/DTJ1)*DFJ1)
    ELSE
        W = (W - H * DLOG(DFM/DTM))*DFM
    ENDIF
ENDIF

IF (DTM < 0.25d0*H) THEN
    E = 0.5d0*(GREENI(X, Y, R0, Z0) + GREENI(X, Y, R1, Z1))
ELSE
    E = GREENI(X, Y, RM, ZM)
ENDIF

F = E*H + 0.25D0*W
F = F/pi

return
end subroutine bint

!---------------------------------------------------------------------
subroutine d2GREN(R0, Z0, R, Z, d2Gdrz, d2Gdzz, d2Gdrr)

! green's function second derivatives at the point (R,Z)
! the source position is (R0,Z0)
! to obtain the real values of derivatives dGdr,dGdz must be devided by pi

implicit none

real*8, intent(in) :: R0, Z0, R, Z
real*8, intent(out) :: d2Gdrz, d2Gdzz, d2Gdrr

integer :: IFAILK, IFAILE
real*8 ::t, tt, fK, fE, Q, fgreen, dtdr, dtdz, dQdt, dGdr, &
    d2tdrr, d2tdrz, d2tdzz, d2qdt2, S21BBF, S21BCF

t = SQRT(4.D0*R*R0/((R + R0)**2 + (Z - Z0)**2))
tt = (1.D0 - t**2)

! THE FUNCTIONS S21BBF, S21BCF FROM NAG LIBRARY.

fK=                 S21BBF(0.D0, tt, 1.D0, IFAILK)
fE = fK - t**2/3.D0*S21BCF(0.D0, tt, 1.D0, IFAILE)

IF (IFAILK /= 0.OR.IFAILE /= 0) WRITE(*, *) ' ATTENTION! GREEN F.:', &
         '  IFAILK, IFAILE= ', IFAILK, IFAILE

Q = ((1.D0 - t**2 * 0.5D0)*fK - fE )/t
fgreen = Q*SQRT(R*R0)
dtdr = 0.5d0*(t/R - t**3 *(R+R0)/(2.d0*R*R0))
dtdz = -t**3*(z - z0)/(4.d0*R*R0)
dQdt = -Q/t + 0.5d0*( fE/tt-fK)
d2tdrz = -(1/r/2 - 3.D0/4.D0*t**2 * (r + r0)/r/r0) * t**3 * (z - z0)/r/r0/4.D0
d2tdzz = 3.D0/16.D0*t**5 * (z - z0)**2/r**2/r0**2 - t**3/r/r0/4.D0
d2tdrr = (1/r/2 - 3.D0/4.D0*t**2 * (r + r0)/r/r0)*(t/r/2 - t**3*(r + r0)/r/r0/4)
d2qdt2 = q/t**2 + fe/(1 - t**2)**2*t + 1/(1 - t**2)*(fe - fk)/t/2 - (fe/(1 - t**2) - fk)
dGdr = 0.5d0*SQRT(R0/R)*Q + SQRT(R0*R)*dQdt*dtdr
d2Gdzz = SQRT(R*R0)*(d2qdt2*dtdz**2 + dQdt*d2tdzz)
d2Gdrz = SQRT(R*R0)*(dQdt*dtdz/r/2.d0 + d2qdt2*dtdr*dtdz + dQdt*d2tdrz)
d2Gdrr = SQRT(R*R0)*(dQdt*dtdr/r/2.d0 + d2qdt2*dtdr*dtdr + dQdt*d2tdrr - &
         q/r**2/2.d0) + dGdr/r/2.d0

return
end subroutine d2gren

!---------------------------------------------------------------------
real*8 function FUNSQ(R1, R2, R3, R4, Z1, Z2, Z3, Z4)

implicit none

real*8, intent(in) :: R1, R2, R3, R4, Z1, Z2, Z3, Z4
real*8 :: R13, R24, Z13, Z24

R13 = R3 - R1
R24 = R4 - R2
Z13 = Z3 - Z1
Z24 = Z4 - Z2
FUNSQ = (R13*Z24 - R24*Z13)*0.5D0

return
end function FUNSQ

!---------------------------------------------------------------------
SUBROUTINE PROG1D(M,U, A,B,C,F, ALF,BET)
!
!     A(J)*U(J+1)-C(J)*U(J)+B(J)*U(J-1)=-F(J) , J=1,M
!     B(1)=0 , A(M)=0

implicit none

integer, intent(in) :: M
real*8, intent(in), dimension(*) :: A, B, C, F
real*8, intent(out), dimension(*) :: U, ALF, BET

integer :: j, k
real*8 :: DK

ALF(1) = A(1)/C(1)
BET(1) = F(1)/C(1)

do K=2, M
   DK = C(K) - ALF(K-1)*B(K)
   ALF(K) = A(K)/DK
   BET(K) = (F(K) + BET(K-1)*B(K))/DK
enddo

U(M) = BET(M)

do J=2, M
    K = M + 1 - J
    U(K) = ALF(K)*U(K+1) + BET(K)
enddo

return
end subroutine PROG1D

!---------------------------------------------------------------------
real*8 function blin_(r0, z0, r1, r2, z1, z2, u1, u2, u3, u4)

implicit none

real*8, intent(in) :: r0, z0, r1, r2, z1, z2, u1, u2, u3, u4
real*8 :: s1, s2, s3, s4

s1 = (r2 - r0)*(z2 - z0)
s3 = (r0 - r1)*(z0 - z1)
s2 = (r0 - r1)*(z2 - z0)
s4 = (r2 - r0)*(z0 - z1)

blin_ = (s1*u1 + s2*u2 + s3*u3 + s4*u4)/(s1 + s2 + s3 + s4)

return
end function blin_

!---------------------------------------------------------------------	
real*8 function bqlin_(r0,z0,x,u)

implicit none

double precision, intent(in) :: r0, z0, x(9, 2), u(9)

integer :: i, j
double precision, dimension(9) :: AA
double precision, dimension(9, 9) :: MM, NN

! r1z1 r2z1 r1z2 r2z2 r1z3 r2z3 r3z3 r3z2 r3z1	
!r**2 z**2 rz r z r2z rz2 const r2z2
do j=1, 9
    MM(j, 1) = x(j, 1)**2
    MM(j, 2) = x(j, 2)**2
    MM(j, 3) = x(j, 1)*x(j, 2)
    MM(j, 4) = x(j, 1)
    MM(j, 5) = x(j, 2)
    MM(j, 6) = x(j, 1)**2 * x(j, 2)
    MM(j, 7) = x(j, 1)*x(j, 2)**2
    MM(j, 8) = 1.
    MM(j, 9) = x(j, 1)**2 * x(j, 2)**2
enddo

call inverse_spid(MM, NN, 9)
do j=1, 9
    AA(j) = 0.
    do i=1, 9
        AA(j) = AA(j) + NN(j, i)*u(i)			
    enddo
enddo

bqlin_ = AA(1) * r0**2 + AA(2) * z0**2 + AA(3) * r0*z0 + &
         AA(4) * r0    + AA(5) * z0    + AA(6) * r0**2 * z0 + &
         AA(7) * r0 * z0**2 + AA(8)    + AA(9) * r0**2 * z0**2

return
end function bqlin_
