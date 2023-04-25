subroutine find_new_X0Y0(Nr, Nt, PSI, X, Y, X0, Y0, & 
& psiax, iax, jax,derivs)

implicit none
integer, intent(in) :: Nr, Nt
double precision, intent(in), dimension(Nr, Nt) :: PSI, X, Y
integer, intent(out) :: iax, jax
double precision, intent(out) :: X0, Y0, psiax
integer :: jt, info, jmin(2),jj,j,ii(nt)
double precision :: g3(2*Nt+1)
double precision :: work(2*(2*Nt+1)*6), matrix(2*Nt+1, 6)
double precision :: derivs(8),c(6),yy(nt)


!jmin = minloc(PSI)
!iax = jmin(1)
!jax = jmin(2)

!if (iax /= 1) then
!    x0 = x(iax, jax)
!    y0 = y(iax, jax)
!    psiax = PSI(iax, jax)
!    return
!endif

	do j=1,nt
		ii(j)=minloc(psi(1:nr,j),1)
		yy(j)=psi(ii(j),j)
	enddo
	jj=minloc(yy(1:nt),1)
	psiax=psi(ii(jj),jj)
	x0=x(ii(jj),jj)
	y0=y(ii(jj),jj)
	iax=ii(jj)
	jax=jj

	if (iax.ne.1) return
	if (jax.ne.1) return


!refine axis 
	g3(2:nt+1)=psi(2,1:nt)
	g3(nt+2:2*nt+1)=psi(3,1:nt)
	g3(1)=psi(1,1)
	matrix(1,1)=x(1,1)		 
	matrix(1,2)=y(1,1)		 
	do jj=1,nt
	matrix(jj+1,1)=x(2,jj)		 
	matrix(jj+1,2)=y(2,jj)		 
	enddo 
	do jj=1,nt
	matrix(nt+jj+1,1)=x(3,jj)		 
	matrix(nt+jj+1,2)=y(3,jj)		 
	enddo 
	call least_square_biquad_ef( &
     & matrix(1:1*nt+1,1),matrix(1:1*nt+1,2), &
     & g3(1:1*nt+1),1*nt+1,c,x0, &
     & y0,psiax,derivs)
	

return
end subroutine find_new_X0Y0
