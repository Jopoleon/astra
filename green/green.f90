!---------------------------------------------------------------------
double precision function ellE_green(X,  DL) ! gives back the first kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

double precision,  intent(in) :: X,  DL

ellE_green = (((0.01736506451D0 *X + 0.04757383546D0)*X + 0.06260601220D0)*X + 0.44325141463D0)*X + 1.0D0 - &
             (((0.00526449639D0 *X + 0.04069697526D0)*X + 0.09200180037D0)*X + 0.24998368310D0)*X*DL

return
end function ellE_green

!---------------------------------------------------------------------
double precision function ellK_green(X,  DL) ! gives back the second kind elliptic integral,  from K. Lackner,  T. Lunt,  IPP - Garching 2022

double precision,  intent(in) :: X,  DL

ellK_green = ((( 0.01451196212D0*X + 0.03742563713D0)*X + 0.03590092383D0)*X + 0.09666344259D0)*X + 1.38629436112D0 - &
            (((( 0.00441787012D0*X + 0.03328355346D0)*X + 0.06880248576D0)*X + 0.12498593597D0)*X + 0.5D0)*DL

return 
end function ellK_green

!----------------------------------------------------------------------------------- 
double precision function green_function(r1, z1, r2, z2)

implicit none

double precision, intent(in) :: r1, z1, r2, z2
double precision :: TT, K, ELCK, ELCE, acl, alg
double precision, external :: ellk_green, elle_green

K = sqrt(4.*r1*r2/ ((r2 + r1)**2 + (z2 - z1)**2))

TT = 1. - K**2

acl = tt
alg = dlog(acl)

ELCK = ellK_green(acl, alg) !              S21BBF(0.D0, TT, 1.D0, IFAILK)
ELCE = ellE_green(acl, alg) ! ELCK-K**2/3.D0*S21BCF(0.D0, TT, 1.D0, IFAILE)

green_function = ( (1.D0 - K**2/2.)*ELCK - ELCE )*( SQRT(r1*r2)/K )

return
end function green_function

!-----------------------------------------------------------------------------------
double precision function green_function_includingsamepoint(r1, z1, r2, z2, dl)

implicit none

double precision, intent(in) :: r1, r2, z1, z2, dl
double precision :: TT
double precision, external :: green_function

if (abs(r1 - r2) < 1.e-6 .and. abs(z1 - z2) < 1.e-6) then
    tt = dl/(4*r1)
    green_function_includingsamepoint = -tt*(log(tt**2) - 4*log(2.) + 2.)*r1**2/dl
else
    green_function_includingsamepoint = green_function(r1, z1, r2, z2)
endif

return
end function green_function_includingsamepoint

!-----------------------------------------------------------------------------------
double precision function green_function_identity(r1, dr, dz, ntype)

implicit none

integer, intent(in) :: ntype
double precision, intent(in) :: r1, dr, dz
integer :: i, j
double precision :: ELCK, ELCE, greenf
double precision, external :: green_function

if (ntype == 2) then
    green_function_identity = r1*(log(8.*r1/(0.2236*(dr + dz))) - 2.0) ! From A. Kavin,  used in SPIDER (A. A. Ivanov and S. Yu. Medvedev),  self inductance of a rectangular coil in toroidal direction
else if (ntype == 1) then
    elce = SQRT(dr**2 + dz**2)
    do i=-1, 1
        do j=-1, i
            if (i == j) then
                greenf = greenf + (r1 + dr*i/3.d0) * (log( 8.*(r1 + dr*i/3.d0)/(elce/3.d0) ) - 0.5D0)
            else
                elck = green_function(r1 + dr*i/3., dz*i/3., r1 + dr*j/3., dz*j/3.)
                greenf = greenf + 4.*elck
            endif
        enddo
    enddo
    green_function_identity = greenf/9.d0
endif

return
end function green_function_identity

!-----------------------------------------------------------------------------------
double precision function green_function_non_identity(r1, z1, r2, z2, dr1, dz1, dr2, dz2, ntype1, ntype2)

implicit none

integer, intent(in) :: ntype1, ntype2
double precision, intent(in) :: r1, z1, r2, z2, dr1, dz1, dr2, dz2

integer :: i, j
double precision :: ELCK, greenf
double precision, external :: green_function

greenf = 0.

if (ntype1 == 1 .and. ntype2 == 1) then
    do i=-1, 1
        do j=-1, 1
            elck = green_function(r1+dr1*i/3., z1+dz1*i/3., r2+dr2*j/3., z2+dz2*j/3.)
            greenf = greenf + elck
        enddo
    enddo
    green_function_non_identity = greenf/9.d0
else if (ntype1 == 1 .and. ntype2 == 2) then
    do i=-1, 1
        elck = green_function(r1+dr1*i/3., z1+dz1*i/3., r2, z2)
        greenf = greenf + elck
    enddo
    green_function_non_identity = greenf/3.d0
else if (ntype1 == 2 .and. ntype2 == 1) then
    do j=-1, 1
        elck = green_function(r1, z1, r2+dr2*j/3., z2+dz2*j/3.)
        greenf = greenf + elck
    enddo
    green_function_non_identity = greenf/3.d0
else if (ntype1 == 2 .and. ntype2 == 2) then
    green_function_non_identity = green_function(r1, z1, r2, z2)
endif

return
end function green_function_non_identity
