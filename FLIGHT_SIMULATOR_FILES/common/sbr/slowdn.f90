      subroutine SLOWDN

!EFable 2016, compute t_fast and p_fast assuming slowing down. R. Bilato, http://pubman.mpdl.mpg.de/pubman/item/escidoc:2039793/component/escidoc:2183257/Bilato_On.pdf, page 4-5
			
	use parameter_inc
	use status_inc
	use const_inc
	use outcmn_inc

			implicit none

	integer J,i,ishot
				
!      include	'for/parameter.inc'
!      include	'for/const.inc'
!      include	'for/status.inc'
!      include	'for/outcmn.inc'
!      include	'tmp/declar.usr'

	double precision zwe(NRD)
	double precision zishot
	
	double precision e_crit(NRD),t_fast(NRD)
	double precision e_birth,gg(NRD),xx(NRD)
	
!alphas
		e_birth=3500.

		e_crit = 4.*TE*( &
     & 3./4.*sqrt(GP)*42.8561* &
     & ZMAIN**2.*(F1+F2)/NE/AMAIN+ &
     & ZIM1**2.*NIZ1/NE/AIM1+ &
     & ZIM2**2.*NIZ2/NE/AIM2+ &
     & ZIM3**2.*NIZ3/NE/AIM3 &
     & )**(2./3.)

	xx=sqrt(e_birth/e_crit)
	gg=1./sqrt(3.)*atan((2.*xx-1.)/sqrt(3.))+ &
     & 1./6.*log((xx**2.-xx+1)/(xx+1)**2.)+ &
     & sqrt(3.*GP)/18.
	
	t_fast=e_crit* &
     & (log(1+(e_birth/e_crit)**(1.5)))**(-1.)* &
     & (e_birth/e_crit-2*gg)

!	write(*,*)	e_crit(1),t_fast(1)

	PFAST=t_fast*NALF !output











						
						end




