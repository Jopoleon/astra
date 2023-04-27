double precision function green_function(r1, z1, r2, z2)

implicit none

double precision, intent(in) :: r1, z1, r2, z2
double precision :: TT, K, acl, alg
double precision, external :: ellk_green, elle_green

K = sqrt(4.*r1*r2/((r2 + r1)**2 + (z2 - z1)**2))
TT = 1. - K**2

acl = tt
alg = dlog(acl)

green_function = ( (1.D0 - K**2/2.)*ellK_green(acl, alg) - ellE_green(acl, alg))*( SQRT(r1*r2)/K )

return
end function green_function

!---------------------------------------------------------------------
double precision function ellE_green(X, DL) ! first kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

double precision, intent(in) :: X, DL

ellE_green = (((0.01736506451D0 *X + 0.04757383546D0)*X + 0.06260601220D0)*X + 0.44325141463D0)*X + 1.0D0 - &
             (((0.00526449639D0 *X + 0.04069697526D0)*X + 0.09200180037D0)*X + 0.24998368310D0)*X*DL

return
end function ellE_green

!-----------------------------------------------------------------------------------
double precision function ellK_green(X, DL) ! second kind elliptic integral, from K. Lackner, T. Lunt, IPP - Garching 2022

double precision, intent(in) :: X, DL

ellK_green= ((( 0.01451196212D0*X + 0.03742563713D0)*X + 0.03590092383D0)*X + 0.09666344259D0)*X + 1.38629436112D0 - &
            ((((0.00441787012D0*X + 0.03328355346D0)*X + 0.06880248576D0)*X + 0.12498593597D0)*X + 0.5D0)*DL

return
end function ellK_green
