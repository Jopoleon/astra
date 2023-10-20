subroutine full_system_advance_ef(j_init)

use errors_params
use ef_circuit
use astra2fbe
use parameters_a2equil
use feqis_tools, only: psi_external_calc, plasma_psi_to_coils

implicit none

integer j_init,j_iter
double precision error_temp	
double precision cur_temp(300)	

! time stepping
! at iteration 0, dpsidt = 0
if (j_init.eq.0) then
write(*,*) 'init full system'
!first do full equilibrium solution at time t=0
call psi_external_calc
call solve_gse2d_fbe_full_ef(0)	
call plasma_psi_to_coils
psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)

write(*,*) 'init done'
call circuit_eq_advance_ef(0)	
j_init=-1
return
endif

! full iterations
write(*,*) 'iter full system'

if (fast_mode.eq.1.and.execute_plasma.eq.1) then
psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)
call psi_external_calc
call solve_gse2d_fbe_full_ef_1turn(1,0,0.d0,0.d0)	
call plasma_psi_to_coils
endif

do j_iter=1,2*max_iter

cur_temp(1:nconduc)=curconduc(1:nconduc)
call circuit_eq_advance_ef(1)	

if (fast_mode.eq.0) then
call psi_external_calc
call solve_gse2d_fbe_full_ef_1turn(1,0,0.d0,0.d0)	
call plasma_psi_to_coils
endif

error_temp=sum(abs(cur_temp(1:nconduc)-curconduc(1:nconduc)))/(nconduc+err_epsilon)/iplasma

if (j_iter.eq.1.and.fast_mode.eq.0) then
error_temp=100.
endif

write(*,*) 'error circuit equation',error_temp,j_iter
if (error_temp.le.err_circ_plasma_iter) then
if (execute_plasma.eq.1.or.fast_mode.eq.0) then
j_init=-1
endif

write(*,*) 'circuit eq done'

return	
endif

if (j_iter.gt.max_iter) then
write(*,*) 'circuit equations not converging, max number of iterations override... stopping',error_temp
stop
endif

enddo

return
end subroutine full_system_advance_ef

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_ef(j_init)	

use errors_params, only: err_find_psistab
use ef_circuit
use astra2fbe
use feqis_tools, only: get_closest_index

implicit none

integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
integer j_iter,j_iter2,iax,jax,j_cyclo
double precision g(300,300),temp_err,raxold,zaxold
double precision g2(300,300)
double precision temp_err2,raxoldo,zaxoldo
double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
double precision dum1,dum2,raxtmp,zaxtmp,dum3,det
double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
double precision curpastmp(npassive),correction
double precision curpastmpo(npassive),psicdtmp(300)
double precision voltdz(300),curref(300),rdtl,psro,pszo,dist1,dist2
double precision przsav(6),psistabrr,psistabzz,cibapr,cibazr
double precision rleft,rright,zup,zdown,dcrdr,dcrdz,dczdr,dczdz
integer jeppa
!first, initialized initial guess coming from prescribed boundary current density: jrhoteta

if (refit_mode.eq.-1) then ! 1 turn only
call solve_gse2d_fbe_full_ef_1turn(j_init,0,0.d0,0.d0)	
return
endif

if (refit_mode.eq.0) then
!start iterations to find self-consistent solution
g000(1:nr2,1:nz2)=psiextrz(1:nr2,1:nz2)
call get_closest_index(raxp,zaxp,iaxis,jaxis)
rax=r(iaxis)
zax=z(jaxis)
raxold=rax
zaxold=zax
raxoldo=rax
zaxoldo=zax
psistabR=0.
psistabZ=0.
psistab1o=0.
psistab2o=0.
temp_err=100.
temp_err2=100.
rleft=raxold
rright=raxold
zup=zaxold
zdown=zaxold
! outer cycle, calculate new axis
j_cyclo=0

jeppa=1
dist1=dr*dr_factor_init
dist2=dz*dz_factor_init

do j_iter2=1,10000000

if (j_cyclo.eq.1) then
raxold=raxoldo+dist1
zaxold=zaxoldo
psistab1o=psistabr
psistab2o=psistabz
endif
if (j_cyclo.eq.2) then
raxold=raxoldo
zaxold=zaxoldo+dist2
dcrdr=(psistabr-psistab1o)/dist1
dczdr=(psistabz-psistab2o)/dist1
endif

if (j_cyclo.eq.3) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
dcrdz=(psistabr-psistab1o)/dist2
dczdz=(psistabz-psistab2o)/dist2
det=(dcrdr*dczdz-dcrdz*dczdr)
cibapr=(dczdz*psistab1o-dcrdz*psistab2o)/det
cibazr=(-dczdr*psistab1o+dcrdr*psistab2o)/det
raxold=raxoldo-cibapr
zaxold=zaxoldo-cibazr
raxoldo=raxold
zaxoldo=zaxold
endif

!inner cycle, calculate psi1 and psi2

psro=1000.
pszo=1000.
redo_bnd=1

do j_iter=1,10000

raxtmp=trax
zaxtmp=tzax

call solve_gse2d_fbe_full_ef_1turn(j_iter-1+j_iter2-1,1,raxold,zaxold)	

temp_err=(abs(psro-psistabr)+abs(pszo-psistabz))

if (temp_err.le.err_find_psistab)	goto 941

psro=psistabr
pszo=psistabz

enddo
write(*,*) 'total iterations passed, stopping!'
stop
return

941 continue

if (j_cyclo.eq.3) then	
temp_err2=abs(psistabr)+abs(psistabZ)
j_cyclo=1
write(*,*) n_of_newton_iterations,nint((0.+j_iter2)/3.)
if (nint((0.+j_iter2)/3.).ge.n_of_newton_iterations) goto 609
else
j_cyclo=j_cyclo+1
endif

if (j_iter2.ge.400000)	goto 609
if (temp_err2.le.err_find_psistab)		goto 609

enddo

609 continue

return

endif	

if (refit_mode.eq.101) then !only vertical stab
!start iterations to find self-consistent solution
g000(1:nr2,1:nz2)=psiextrz(1:nr2,1:nz2)
call get_closest_index(raxp,zaxp,iaxis,jaxis)
rax=r(iaxis)
zax=z(jaxis)
raxold=rax
zaxold=zax
raxoldo=rax
zaxoldo=zax
psistabR=0.
psistabZ=0.
psistab1o=0.
psistab2o=0.
temp_err=100.
temp_err2=100.
rleft=raxold
rright=raxold
zup=zaxold
zdown=zaxold
! outer cycle, calculate new axis
j_cyclo=0

jeppa=1
dist1=dr
dist2=dz

do j_iter2=1,10000000
if (j_iter2.gt.250) stop

if (j_cyclo.eq.1) then
zaxold=zaxoldo+dist2
psistab2o=psistabz
endif
if (j_cyclo.eq.2) then      ! inverse of (dcrdr dcrdz  ; dczdr dczdz ) = ( dczdz -dcrdz ; -dczdr dcrdr)
dczdz=(psistabz-psistab2o)/dist2
cibazr=psistab2o/dczdz
zaxold=zaxoldo-cibazr
zaxoldo=zaxold
endif

!inner cycle, calculate psi1 and psi2

psro=1000.
pszo=1000.
redo_bnd=1

do j_iter=1,10000

raxtmp=trax
zaxtmp=tzax

call solve_gse2d_fbe_full_ef_1turn(j_iter-1+j_iter2-1,1,-1.d6,zaxold)	

temp_err=(abs(pszo-psistabz))

if (temp_err.le.err_find_psistab)	goto 9411

psro=psistabr
pszo=psistabz

enddo
write(*,*) 'total iterations passed!'
stop
return

9411 continue

if (j_cyclo.eq.2) then	
temp_err2=abs(psistabZ)
j_cyclo=1
write(*,*) n_of_newton_iterations,nint((0.+j_iter2)/2.)
if (nint((0.+j_iter2)/2.).ge.n_of_newton_iterations) goto 6091
else
j_cyclo=j_cyclo+1
endif

if (j_iter2.ge.400000)	goto 6091
if (temp_err2.le.err_find_psistab)		goto 6091

enddo

6091 continue

write(*,*) 'total iterations passed!'

return

endif	

if (refit_mode.eq.1) then
call restab_axis_with_furier_wall
return	
endif	

if (refit_mode.eq.2) then
call restab_boundary_with_furier_wall !doesnt work well
return	
endif	

if (refit_mode.eq.3) then
call restab_F_function_full_fonfit
return	
endif	

return
end subroutine solve_gse2d_fbe_full_ef

!--------------------------------------------------------------------	
subroutine restab_F_function_full_fonfit

!refits all currents
use errors_params, only: err_find_psistab
use ef_circuit
use astra2fbe
use green_matrix
use feqis_tools, only: interp_j_fromrhotorz, get_closest_index, &
    find_fields_interp_green, inverse_matrix, boundary, &
    psi_external_calc, find_fields_interp_psionly

implicit none

integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
integer j_iter,j_iter2,iax,jax
double precision g(300,300),temp_err,raxold,zaxold
double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
double precision psistab1,psistab2,dum1,dum2,zum1,zum2
double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
double precision voltdz(300),rdtl,tin,tup,cum1,cum2
double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
double precision G_00(nconduc,nteta)
double precision G_00c(nconduc) !average
double precision G_00r(nconduc)
double precision G_00z(nconduc)
double precision f_correction,x1,x2,x3,x4
double precision psibt0,psibt1,psicorr(nteta)
double precision ggep(nteta)
double precision matrix(nconduc,nconduc)
double precision invmatrix(nconduc,nconduc)
double precision Fderiv(nconduc),Ffunc,Ffunc_old
integer info,whichcoil
double precision curref(nconduc),curnow(nconduc),curdiff(nconduc)
double precision raxref,zaxref,rbref(500),zbref(500)
data cum1/0./
data cum2/0./
save cum1,cum2,dum8,dum9

psicorr=0.

Ffunc_old=1.e6

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
write(*,*) 'reinterp curr, restab'
call interp_j_fromrhotorz
!rescale plasma
dum1=0.
do j=1,nz2
do i=1,nr2
dum1=dum1+jrz(i,j)*dr*dz
enddo
enddo
jrz=jrz/dum1*iplasma

rax=raxp
zax=zaxp
write(*,*) raxp,zaxp,dum1

call get_closest_index(rax,zax,iaxis,jaxis)
iax=iaxis
jax=jaxis

!axis is given by raxp, zaxp, boundary by rbndp, zbndp, coilref by curconduc(passive)

sigma_coils(nactive+1:nconduc)=sigma_coils(nactive+1)

curref(1:nconduc)=curconduc(1:nconduc)
curnow(1:nconduc)=curconduc(1:nconduc)
curdiff=0.

raxref=raxp
zaxref=zaxp
rbref(1:nteta)=rbndp(1:nteta)
zbref(1:nteta)=zbndp(1:nteta)

! calculate the matrix F_li of the F function, including the green function terms
matrix=0.
invmatrix=0.
do j=1,nconduc
do k=1,nteta
call find_fields_interp_green(rbref(k),zbref(k),G_00(j,k),j) !give back psi,br,bz at r0,z0
enddo
G_00c(j)=sum(G_00(j,1:nteta))/(0.+nteta)
call find_fields_interp_green(raxref-dr/2.,zaxref,bub(1),j) !give back psi,br,bz at r0,z0
call find_fields_interp_green(raxref+dr/2.,zaxref,bub(2),j) !give back psi,br,bz at r0,z0
call find_fields_interp_green(raxref,zaxref-dz/2.,bub(3),j) !give back psi,br,bz at r0,z0
call find_fields_interp_green(raxref,zaxref+dz/2.,bub(4),j) !give back psi,br,bz at r0,z0
G_00r(j)=(bub(2)-bub(1))/dr
G_00z(j)=(bub(4)-bub(3))/dz
enddo

do j=1,nconduc
do i=1,nconduc
if (i.eq.j) matrix(i,j)=matrix(i,j)+2.*sigma_coils(i)
matrix(i,j)=matrix(i,j)+ & 
& 2.*sigma_B*sum((G_00(i,1:nteta)-G_00c(i))*(G_00(j,1:nteta)-G_00c(j)))+ &
& 2.*sigma_axis*(G_00r(i)*G_00r(j)+G_00z(i)*G_00z(j))
enddo
enddo

!calculate inverse
call inverse_matrix(matrix,invmatrix,nconduc)

do j_iter=1,300000 !iterations to find currents

if (j_iter.gt.50) stop

CALL CPU_TIME(tin)	
g=0.
call solve_gs2d(g) !jrz as right hand side
CALL CPU_TIME(tup)

write(*,*) 'stop here2'
!find boundary condition using g
CALL CPU_TIME(tin)	
call boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary
write(*,*) 'stop here3'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
write(*,*) 'stop here4'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	

!construct correction
do j=1,nz2
do i=1,nr2
f_correction=0.
do k=1,nconduc
f_correction=f_correction+curdiff(k)*greeni(i,j,k) 
enddo
psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
enddo
enddo

CALL CPU_TIME(tup)
write(*,*) tup-tin	

do j=1,nteta
xub(1)=rbref(j)
yub(1)=zbref(j)
call find_fields_interp_psionly(xub(1),yub(1),psicorr(j))  !psi on the boundary
enddo

x1=sum(psicorr)/(nteta+0.) !average psi on the boundary

!derivative at ref axis
call find_fields_interp_psionly(raxref-dr/2.,zaxref,bub(1)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxref+dr/2.,zaxref,bub(2)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxref,zaxref-dz/2.,bub(3)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxref,zaxref+dz/2.,bub(4)) !give back psi,br,bz at r0,z0
x2 = (bub(2)-bub(1))/dr ! psir
x3 = (bub(4)-bub(3))/dz ! psiz

Ffunc=sigma_B*sum((psicorr-x1)**2.)+sum(sigma_coils(1:nconduc)*curdiff**2.)+sigma_axis*(x2**2.+x3**2.)

!calculate F derivative
do i=1,nconduc
Fderiv(i)=2.*(sigma_coils(i)*curdiff(i)+ & 
& sigma_B*sum((psicorr-x1)*(G_00(i,1:nteta)-G_00c(i)))+ &
& sigma_axis*(x2*G_00r(i)+x3*G_00z(i)) & 
& )
enddo

!calculate new currents

do i=1,nconduc
curnow(i)=curnow(i)-sum(invmatrix(i,1:nconduc)*Fderiv)
enddo	
curdiff=curnow-curref

!construct correction
do j=1,nz2
do i=1,nr2
f_correction=0.
do k=1,nconduc
f_correction=f_correction+curdiff(k)*greeni(i,j,k) 
enddo
psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
enddo
enddo

CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) 'stop here41'
CALL CPU_TIME(tin)	
call find_new_axis_part1	
CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

CALL CPU_TIME(tin)	
write(*,*) 'stop here55'
call find_psi_boundary
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call new_jrz_ef  ! calculate new right hand side
write(*,*) 'stop here7'
CALL CPU_TIME(tup)
write(*,*) 'current',tup-tin	

temp_err=abs(Ffunc-Ffunc_old)

Ffunc_old=Ffunc

write(*,*) 'temp err',temp_err,Ffunc

if (temp_err.le.err_find_psistab) goto 1032

enddo

1032 continue

write(*,*) raxp,zaxp,rax,zax

!check
do i=1,nconduc
curconduc(i)= curnow(i)
enddo

call psi_external_calc
psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
call find_new_axis_part1	
call find_psi_boundary
call new_jrz_ef  ! calculate new right hand side

write(*,*) curconduc(1:nconduc),rax,zax	
write(*,*) 'full fonfit eddy currents converged'

return
end subroutine restab_F_function_full_fonfit

!--------------------------------------------------------------------
subroutine restab_boundary_with_furier_wall !not working well

use errors_params, only: err_find_psistab
use ef_circuit
use astra2fbe
use green_matrix
use feqis_tools, only: interp_j_fromrhotorz, get_closest_index, &
    find_angle, find_fields_interp_green, boundary, &
    find_fields_interp_psionly, psi_external_calc, &
    least_square_biquad

implicit none

integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
integer j_iter,j_iter2,iax,jax
double precision g(300,300),temp_err,raxold,zaxold
double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
double precision psistab1,psistab2,dum1,dum2,zum1,zum2
double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
double precision curpastmp(npassive),correction,delr,delz
double precision curpastmpo(npassive),psicdtmp(300),anglr(npassive)
double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
double precision g_0(npassive),g0_r(npassive),g0_z(npassive)
double precision S_00(n_fourier_restab_boundary),C_00(n_fourier_restab_boundary)
double precision G_00c(nteta,n_fourier_restab_boundary)
double precision G_00s(nteta,n_fourier_restab_boundary)
double precision f_correction,x1,x2,x3,x4
double precision psibt0,psibt1,psicorr(nteta)
double precision matrix(nteta,2*n_fourier_restab_boundary)
integer info
double precision work(4*nteta*n_fourier_restab_boundary)

data cum1/0./
data cum2/0./
save cum1,cum2,dum8,dum9

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
write(*,*) 'reinterp curr, restab'
call interp_j_fromrhotorz
!rescale plasma
dum1=0.
do j=1,nz2
do i=1,nr2
dum1=dum1+jrz(i,j)*dr*dz
enddo
enddo
jrz=jrz/dum1*iplasma

rax=raxp
zax=zaxp
write(*,*) raxp,zaxp,dum1


call get_closest_index(rax,zax,iaxis,jaxis)
iax=iaxis
jax=jaxis

G_00c=0.
G_00s=0.

! evaluate coils things
do i=1,npassive
call find_angle(rax,zax,r_cond(nactive+i),z_cond(nactive+i),anglr(i))
!find true axis
do j=1,nteta
xub(1)=rbndp(j)
yub(1)=zbndp(j)
call find_fields_interp_green(xub(1),yub(1),bub(2),nactive+i) !give back psi,br,bz at r0,z0
do k=1,n_fourier_restab_boundary
G_00c(j,k)=G_00c(j,k)+cos(k*anglr(i))*bub(2)
G_00s(j,k)=G_00s(j,k)+sin(k*anglr(i))*bub(2)
enddo
enddo
enddo

do j=1,nteta
	do k=1,n_fourier_restab_boundary
matrix(j,k)=G_00c(j,k)
matrix(j,n_fourier_restab_boundary+k)=G_00s(j,k)
enddo
enddo

psibt0=1000.
psibt1=1000.

do j_iter=1,30

CALL CPU_TIME(tin)	
g=0.
call solve_gs2d(g) !jrz as right hand side
CALL CPU_TIME(tup)
write(*,*) 'solve g0',tup-tin	

write(*,*) 'stop here2'
!find boundary condition using g
CALL CPU_TIME(tin)	
call boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary
write(*,*) 'stop here3'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
write(*,*) 'stop here4'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	

psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux

do j=1,nteta
xub(1)=rbndp(j)
yub(1)=zbndp(j)
call find_fields_interp_psionly(xub(1),yub(1),psicorr(j)) 
enddo

x1=sum(psicorr)/(nteta+0.)
psibt1=x1
psicorr=x1-psicorr

!least square fit solution

call dgels('N', &
     &		nteta,&
     &  	n_fourier_restab_boundary*2,&
     &  	1,&
     &  		matrix,&
     &  	nteta,&
     &  	psicorr,&
     &  	nteta,&
     &  	WORK,&
     &  	2*(nteta)*n_fourier_restab_boundary*2,&
     &  	INFO ) 		 

!construct correction
do j=1,nz2
do i=1,nr2
f_correction=0.
do k=1,n_fourier_restab_boundary
f_correction=f_correction+ & 
& psicorr(k)*sum(greeni(i,j,nactive+1:nactive+npassive)*cos(k*anglr(1:npassive)))	+ & 
& psicorr(n_fourier_restab_boundary+k)*sum(greeni(i,j,nactive+1:nactive+npassive)*sin(k*anglr(1:npassive)))
enddo
psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+f_correction !total flux
enddo
enddo

CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) 'stop here41'
CALL CPU_TIME(tin)	
call find_new_axis_part1	
CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

CALL CPU_TIME(tin)	
write(*,*) 'stop here55'
call find_psi_boundary
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call new_jrz_ef  ! calculate new right hand side
write(*,*) 'stop here7'
CALL CPU_TIME(tup)
write(*,*) 'current',tup-tin	

temp_err=abs(psibt0-psibt1)

write(*,*) 'temp err',temp_err,psibt0,psibt1

psibt0=psibt1

if (temp_err.le.err_find_psistab) goto 1032

enddo

1032 continue

write(*,*) raxp,zaxp,rax,zax,psistabr,psistabz	

!check
do i=1,npassive
do j=1,n_fourier_restab_boundary
curconduc(nactive+i)= curconduc(nactive+i)+& 
 & C_00(j)*cos(j*anglr(i))+S_00(j)*sin(j*anglr(i))
enddo
enddo

call psi_external_calc
psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
call find_new_axis_part1	
call find_psi_boundary
call new_jrz_ef  ! calculate new right hand side

write(*,*) curconduc(1:nconduc),rax,zax	

stop
return
end subroutine restab_boundary_with_furier_wall

!--------------------------------------------------------------------
subroutine restab_axis_with_furier_wall

use errors_params, only: err_find_psistab
use ef_circuit
use astra2fbe
use green_matrix       ! declaration of minimal CPOs
use feqis_tools, only: interp_j_fromrhotorz, get_closest_index, &
    find_angle, least_square_biquad, boundary, &
    find_fields_interp_psionly, psi_external_calc

implicit none

integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
integer j_iter,j_iter2,iax,jax
double precision g(300,300),temp_err,raxold,zaxold
double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
double precision psistab1,psistab2,dum1,dum2,zum1,zum2
double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
double precision curpastmp(npassive),correction,delr,delz
double precision curpastmpo(npassive),psicdtmp(300),anglr(npassive)
double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
double precision bub(9),xub(9),yub(9),ccc(6),ddipsi(8)
double precision g_0(npassive),g0_r(npassive),g0_z(npassive)
double precision S_00r,C_00r,S_00z,C_00z,C_00(258,258),S_00(258,258)
data cum1/0./
data cum2/0./
save cum1,cum2,dum8,dum9

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
write(*,*) 'reinterp curr, restab'
call interp_j_fromrhotorz
!rescale plasma
dum1=0.
do j=1,nz2
do i=1,nr2
dum1=dum1+jrz(i,j)*dr*dz
enddo
enddo
jrz=jrz/dum1*iplasma

rax=raxp
zax=zaxp
write(*,*) raxp,zaxp,dum1

call get_closest_index(rax,zax,iaxis,jaxis)
iax=iaxis
jax=jaxis

! evaluate coils things
do i=1,npassive
call find_angle(rax,zax,r_cond(nactive+i),z_cond(nactive+i),anglr(i))
!find true axis
xub(1)=r(iax-1)
xub(2)=r(iax)
xub(3)=r(iax+1)
xub(4)=r(iax)
xub(5)=r(iax)
xub(6)=r(iax-1)
xub(7)=r(iax-1)
xub(8)=r(iax+1)
xub(9)=r(iax+1)
yub(1)=z(jax)
yub(2)=z(jax)
yub(3)=z(jax)
yub(4)=z(jax-1)
yub(5)=z(jax+1)
yub(6)=z(jax-1)
yub(7)=z(jax+1)
yub(8)=z(jax-1)
yub(9)=z(jax+1)
bub(1)=greeni(iax-1,jax,nactive+i)
bub(2)=greeni(iax,jax,nactive+i)
bub(3)=greeni(iax+1,jax,nactive+i)
bub(4)=greeni(iax,jax-1,nactive+i)
bub(5)=greeni(iax,jax+1,nactive+i)
bub(6)=greeni(iax-1,jax-1,nactive+i)
bub(7)=greeni(iax-1,jax+1,nactive+i)
bub(8)=greeni(iax+1,jax-1,nactive+i)
bub(9)=greeni(iax+1,jax+1,nactive+i)
call least_square_biquad(xub,yub,bub,9,ccc,dum1,dum2,zum1,ddipsi)
g0_r(i)=ddipsi(1)
g0_z(i)=ddipsi(2)
enddo

do j=1,nz2
do i=1,nr2
C_00(i,j)=sum(greeni(i,j,nactive+1:nactive+npassive)*cos(anglr))
S_00(i,j)=sum(greeni(i,j,nactive+1:nactive+npassive)*sin(anglr))
enddo
enddo
C_00r=sum(g0_r*cos(anglr))
S_00r=sum(g0_r*sin(anglr))
C_00z=sum(g0_z*cos(anglr))
S_00z=sum(g0_z*sin(anglr))

write(*,*) 'coils structure',c_00(iax,jax),s_00(iax,jax),c_00r,s_00r,c_00z,s_00z,iplasma

psistab1o=1000.
psistab2o=1000.

do j_iter=1,30000

CALL CPU_TIME(tin)	
g=0.
call solve_gs2d(g) !jrz as right hand side
CALL CPU_TIME(tup)
write(*,*) 'solve g0',tup-tin	

write(*,*) 'stop here2'
!find boundary condition using g
CALL CPU_TIME(tin)	
call boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary
write(*,*) 'stop here3'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
write(*,*) 'stop here4'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
psistabr=0.
psistabz=0.
delr=0.
delz=0.
psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux
call find_new_axis_part1	

dum1=C_00r*S_00z-C_00z*S_00r	
call find_fields_interp_psionly(raxp+dr,zaxp,bub(1)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxp-dr,zaxp,bub(2)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxp,zaxp+dz,bub(3)) !give back psi,br,bz at r0,z0
call find_fields_interp_psionly(raxp,zaxp-dz,bub(4)) !give back psi,br,bz at r0,z0

xub(1)=(bub(1)-bub(2))/(2.*dr)
yub(1)=(bub(3)-bub(4))/(2.*dz)

delr=(S_00z*(-xub(1))+(-C_00z)*(-yub(1)))/dum1
delz=(-S_00r*(-xub(1))+C_00r*(-yub(1)))/dum1

psistabr=delr
psistabz=delz

write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)+ &
       & psistabr*C_00(1:nr2,1:nz2)+psistabz*S_00(1:nr2,1:nz2) !total flux

CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) 'stop here41'
CALL CPU_TIME(tin)	
call find_new_axis_part1	
CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) j_iter,rax,zax,raxp,zaxp,psistabr,psistabz

CALL CPU_TIME(tin)	
write(*,*) 'stop here55'
call find_psi_boundary
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call new_jrz_ef  ! calculate new right hand side
write(*,*) 'stop here7'
CALL CPU_TIME(tup)
write(*,*) 'current',tup-tin	

temp_err=(abs(psistab1o-psistabr)+abs(psistab2o-psistabz))

write(*,*) 'temp err',temp_err,psistab1o,psistab2o,psistabr,psistabz

psistab1o=psistabr
psistab2o=psistabz

if (temp_err.le.err_find_psistab) goto 1032

enddo

1032 continue

write(*,*) raxp,zaxp,rax,zax,psistabr,psistabz	

!check
do i=1,npassive
curconduc(nactive+i)= curconduc(nactive+i)+& 
 & psistabr*cos(anglr(i))+psistabz*sin(anglr(i))
enddo

call psi_external_calc
psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2)
call find_new_axis_part1	
call find_psi_boundary
call new_jrz_ef  ! calculate new right hand side

write(*,*) curconduc(1:nconduc),rax,zax	

write(*,*) 'code restab axis with fourier wall converged'

return
end subroutine restab_axis_with_furier_wall

!--------------------------------------------------------------------
subroutine solve_gse2d_fbe_full_ef_1turn(j_init,j_stab,raxold,zaxold)	

use ef_circuit
use astra2fbe
use feqis_tools, only: interp_j_fromrhotorz, get_closest_index, &
    boundary, nine_point_coeffs_only, find_angle, &
    find_fields_interp_psionly

implicit none

integer j_init,i,j,k,ii,jj,kk,iii,jjj,kkk
integer j_iter,j_iter2,iax,jax,j_stab,j_count
double precision g(300,300),temp_err,raxold,zaxold
double precision g000(300,300),br1,bz1,dum6,dum7,dum8,dum9,dum10
double precision psistab1,psistab2,dum1,dum2,zum1,zum2
double precision psistab1o,psistab2o,dbrr,dbrz,dbzr,dbzz
double precision curpastmp(npassive),correction,delr,delz
double precision curpastmpo(npassive),psicdtmp(300)
double precision voltdz(300),curref(300),rdtl,tin,tup,cum1,cum2
double precision c(9)

data cum1/0./
data cum2/0./
data j_count/0/
save cum1,cum2,dum8,dum9,j_count

!first, initialized initial guess coming from prescribed boundary current density: jrhoteta
if (j_init.eq.0) then
write(*,*) 'reinterp curr'
call interp_j_fromrhotorz
!rescale plasma
dum1=0.
do j=1,nz2
do i=1,nr2
dum1=dum1+jrz(i,j)*dr*dz
enddo
enddo
jrz=jrz/dum1*iplasma

rax=raxp
zax=zaxp
call get_closest_index(rax,zax,iaxis,jaxis)
rax=r(iaxis)
zax=z(jaxis)
write(*,*) raxp,zaxp,rax,zax
endif

CALL CPU_TIME(tin)	
g=0.
call solve_gs2d(g) !jrz as right hand side
CALL CPU_TIME(tup)
write(*,*) 'solve g0',tup-tin	

write(*,*) 'stop here2'
!find boundary condition using g
CALL CPU_TIME(tin)	
call boundary(g)  ! gbound = integral (Green*dg/dn) over the boundary
write(*,*) 'stop here3'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call solve_gs2d(g) ! again jrz as right hand side
psiplasrz(1:nr2,1:nz2)=g(1:nr2,1:nz2)
write(*,*) 'stop here4'
CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) 'stop here5'

CALL CPU_TIME(tin)	
if (j_stab.eq.1) then
psistabr=0.
psistabz=0.
delr=0.
delz=0.
psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux

call find_new_axis_part1	
write(*,*) 'natural ax',rax,zax

 call nine_point_coeffs_only(raxold, zaxold, c, zum1,zum2)
!dpsidr
zum1=(raxold-zum1)/dr
zum2=(zaxold-zum2)/dz

dum1=2.*c(2)*(zum1*zum2**2.+2.*c(2)*zum1*zum2)+ &
   c(3)*zum2**2.+c(4)*zum2+2*c(5)*zum1+c(7)
 
!dpsidz
dum2=2.*c(1)*zum1**2.*zum2+c(2)*zum1**2.+ &
   2.*c(3)*zum1*zum2+c(4)*zum1+2.*c(6)*zum2+c(8)

psistabr=-1./(2.*raxold)*dum1/dr
psistabz=-dum2/dz

do i=1,nr2
do j=1,nz2
psirz(i,j)=psiplasrz(i,j)+psiextrz(i,j)+psistabr*r(i)**2.+psistabz*z(j) !total flux
enddo
enddo
call find_new_axis_part1	

write(*,*) 'axis done' ! not converged, stop'
write(*,*) 'raxx',rax,zax,raxold,zaxold,psistabr,psistabz

else

psirz(1:nr2,1:nz2)=psiplasrz(1:nr2,1:nz2)+psiextrz(1:nr2,1:nz2) !total flux

endif
CALL CPU_TIME(tup)
write(*,*) tup-tin	

write(*,*) 'stop here41'
CALL CPU_TIME(tin)	
call find_new_axis_part1	
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
write(*,*) 'stop here55'
call find_psi_boundary
CALL CPU_TIME(tup)
write(*,*) tup-tin	

CALL CPU_TIME(tin)	
call new_jrz_ef  ! calculate new right hand side
write(*,*) 'stop here7'
CALL CPU_TIME(tup)
write(*,*) 'current',tup-tin,nr2,nz2	

return
end subroutine solve_gse2d_fbe_full_ef_1turn

!--------------------------------------------------------------------
subroutine FEQISUPDATE(machine,coilzzz,time_nowz,nccc)
       
use pi_vars, only: GPI2
use ef_circuit
use astra2fbe

implicit none

integer nccc
character*4 machine
double precision coilzzz(nccc),time_nowz

coilzzz(1:nccc)=curconduc(1:nccc)*1.e3
cur_con_old(1:nconduc)=curconduc(1:nconduc)

if (fast_mode.eq.0) then
psi_cur_old(1:nconduc)=psiplasmatoconduc(1:nconduc)
endif

return
end subroutine FEQISUPDATE

!--------------------------------------------------------------------
subroutine circuit_eq_advance_ef(j_init)

use pi_vars, only: GPI, GPI2
use ef_circuit
use astra2fbe       ! declaration of minimal CPOs

implicit none
integer j_init,i,ic,j,k,iii,jjj,invertcommand,i_equivalence	
double precision restemp(200,200),indtemp(200,200)
double precision dpctemp(i_dim1),vtemp(i_dim1)
double precision curotemp(i_dim1),curtemp(i_dim1)
integer rem_coils(200),i_cnew,firstcall

data firstcall/0/
save restemp,indtemp,rem_coils,i_cnew,firstcall

tau_new=tau_circuit_ef
invertcommand=0

if (tau_new.ne.tau_old) then
invertcommand=1
tau_old=tau_new
endif
if (j_init.eq.0.or.firstcall.eq.0) then
invertcommand=1
tau_old=tau_new
rem_coils=0
endif

ic=nconduc

if (j_init.eq.0) cur_con_old=curconduc
if (j_init.eq.0) i_cnew=ic
if (j_init.eq.0) return

firstcall=1

dpc(1:ic)=GPI2*(psiplasmatoconduc(1:ic)-psi_cur_old(1:ic))/tau_gseq_ef !plasma contribution

do j=1,nconduc
if (activate_coil_ef(j).eq.0) cur_con_old(j)=0.
if (activate_coil_ef(j).eq.0) dpc(j)=0.
enddo

if (reconnect_circuits.eq.1) then
invertcommand=1
i_cnew=ic
indtemp(1:ic,1:ic)=indconduc(1:ic,1:ic)
restemp(1:ic,1:ic)=resconduc(1:ic,1:ic)
do k=1,n_equivalence
i_equivalence=1
do i=1,nconduc
if (new_equivalence(i,k).gt.0.and.i_equivalence.eq.1) then
	i_equivalence=0

	do j=i+1,nconduc
		if (new_equivalence(j,k).eq.new_equivalence(i,k)) then
			indtemp(i,i)=indconduc(i,i)+indconduc(j,j)+indconduc(i,j)+indconduc(j,i)
			restemp(i,i)=resconduc(i,i)+resconduc(j,j)+resconduc(i,j)+resconduc(j,i)
			do iii=1,nconduc
				if (iii.ne.i.and.iii.ne.j) then
					indtemp(i,iii)=indconduc(i,iii)+indconduc(j,iii)
					restemp(i,iii)=resconduc(i,iii)+resconduc(j,iii)
					indtemp(iii,i)=indconduc(iii,i)+indconduc(iii,j)
					restemp(iii,i)=resconduc(iii,i)+resconduc(iii,j)
				endif
			enddo
			rem_coils(j)=new_equivalence(j,k)
		ic=ic-1
		endif			
	enddo	
!reove obsolete columns and rows
do jjj=1,nconduc
j=nconduc-jjj+1
if (rem_coils(j).gt.0) then
			indtemp(1:nconduc,j:nconduc)=indtemp(1:nconduc,j+1:nconduc+1)
			restemp(1:nconduc,j:nconduc)=restemp(1:nconduc,j+1:nconduc+1)
			indtemp(j:nconduc,1:nconduc)=indtemp(j+1:nconduc+1,1:nconduc)
			restemp(j:nconduc,1:nconduc)=restemp(j+1:nconduc+1,1:nconduc)
endif
enddo

endif
enddo	
enddo
i_cnew=ic
endif

if (use_reduce_circuit.eq.1) then
ic=nconduc
dpctemp(1:ic)=dpc(1:ic)
vtemp(1:ic)=voltage(1:ic)
curotemp(1:ic)=cur_con_old(1:ic)
do k=1,n_equivalence
i_equivalence=1
do i=1,nconduc
if (new_equivalence(i,k).gt.0.and.i_equivalence.eq.1) then
	i_equivalence=0
	do jjj=i+1,nconduc
		j=nconduc-jjj+i+1
		if (new_equivalence(j,k).eq.new_equivalence(i,k)) then
			dpctemp(i)=dpctemp(i)+dpctemp(j)
		ic=ic-1
		endif			
	enddo	
endif
enddo	
enddo
do jjj=1,nconduc
j=nconduc-jjj+1
if (rem_coils(j).gt.0) then
			dpctemp(j:nconduc)=dpctemp(j+1:nconduc+1)	
			vtemp(j:nconduc)=vtemp(j+1:nconduc+1)	
			curotemp(j:nconduc)=curotemp(j+1:nconduc+1)
endif
enddo
endif

i=i_cnew

if (use_reduce_circuit.eq.0) then
call solve_circuit_equations(i,& 
& indconduc(1:i,1:i),resconduc(1:i,1:i), & 
& cur_con_old(1:i),curconduc(1:i),& 
& 	voltage(1:i),dpc(1:i),tau_new,invertcommand)
else

call solve_circuit_equations(i,& 
& indtemp(1:i,1:i),restemp(1:i,1:i), & 
& curotemp(1:i),curtemp(1:i),& 
& 	vtemp(1:i),dpctemp(1:i),tau_new,invertcommand)	

!readapt currents
do i=1,nconduc
if (rem_coils(i).eq.0) then
		curconduc(i)=curtemp(i)
endif
if (rem_coils(i).gt.0) then
	curconduc(i)=curtemp(rem_coils(i))
	curtemp(i+1:i_cnew+1)=curtemp(i:i_cnew)
endif
enddo
endif

ic=nconduc
i=ic
do j=1,ic
if (activate_coil_ef(j).eq.0) curconduc(j)=0.
curconduc(j)=max(curconduc(j),current_limit_ef(j,2))
curconduc(j)=min(curconduc(j),current_limit_ef(j,1))
if (sum(force_coil(j,1:i)).gt.0.5) then
do k=1,i
	if (force_coil(j,k).eq.1.) then
		curconduc(j)=curconduc(k)
	endif
enddo
endif
enddo

reconnect_circuits=0

return
end subroutine circuit_eq_advance_ef

!--------------------------------------------------------------------
subroutine solve_circuit_equations(nc,im,rm,I0,I1,& 
& 	V,dpc,tau,invertcommand)

use feqis_tools, only: inverse_matrix

implicit none

integer i,j,k,nc,invertcommand
double precision im(nc,nc),rm(nc,nc),i0(nc), &
    & i1(nc),v(nc),dpc(nc),tau,matrix(nc,nc)	
double precision b(nc),invmatrix(200,200)

save invmatrix

!equation is im*(i1-i0)/tau + rm*i1 = v-dpc

do i=1,nc
b(i) = v(i)-dpc(i)+sum(im(i,1:nc)*i0(1:nc))/tau
enddo

if (invertcommand.eq.1) then
matrix(1:nc,1:nc) = im(1:nc,1:nc)/tau+rm(1:nc,1:nc)
call	inverse_matrix(matrix,invmatrix(1:nc,1:nc),nc)
else
endif

do i=1,nc
i1(i) = sum(invmatrix(i,1:nc)*b(1:nc))
enddo

return
end subroutine solve_circuit_equations

!--------------------------------------------------------------------
subroutine definitions_ef_equil(equil_in,params,j_call,ifplasma)

use pi_vars, only: GPI, GPI2
use errors_params
use parameters_a2equil
use imas_ids
use ef_circuit
use astra2fbe
use numerical_tools, only: linterp
use feqis_tools, only: find_angle

implicit none

integer j_call,i,j,k,ifplasma,idum(10)
double precision rgeaz,zgeaz,dum1(10)
double precision rdum(700),zdum(700),tdum(700)

type(type_parameters) params
       type(type_equilibrium) equil_in

if (j_call.eq.0) then
data_dir=params%prename
nteta=equil_in%eqgeometry%boundary%npoints
nrho=params%neql
psistabR=0.
psistabZ=0.
!normalized psi from 0 axis to 1 edge, equispaced
do i=1,nrho
psigrid(i)=(i-1.)/(nrho-1.)
enddo

dr_factor_init=dr_factor_init_astra
dz_factor_init=dz_factor_init_astra

!errors
 err_circ_plasma_iter = err_circ_in
 err_find_oxpoints = err_find_oxpoints_in
 err_find_oxpoints_derivs = err_find_oxpoints_derivs_in
 err_find_psistab = err_find_psistab_in
 err_find_delr = err_find_delr_in
 err_find_biquad = err_find_biquad_in
 err_epsilon =err_epsilon_in
 err_gaptolez = err_gaptolez_in
 err_fix_boundary = err_fix_boundary_in

mu0=0.4*GPI
Rgeom0=equil_in%global_param%toroid_field%r0
tau_circuit_ef=0.001 !default value	
tau_gseq_ef=0.001	 !default value
activate_coil_ef=1. ! when 0., coil is forced to 0 current
current_limit_ef(:,1)=1e6 ! cant be higher than 1e6 MA
current_limit_ef(:,2)=-1e6 ! cant be lower than -1e6 MA
max_iter=10000 !hardwired
! teta for polar grid, goes from 0 to 2*pi-dteta, but point nt+1 is the periodic one
omega_pl=0.
psi0_astra=equil_in%profiles_1d%psi(1)
psib_astra=equil_in%profiles_1d%psi(nrho)
raxp=raxis_astra
zaxp=zaxis_astra	
psibnd=-1.e6
use_limiter_yesno=use_limiter_astra
voltage_old=0.
voltage=0.
endif

if (ifplasma.eq.1) then

dteta=GPI2/(nteta+1+0.)
do i=1,nteta+1
teta(i)=GPI2*(i-1.)/(nteta+0.)
enddo

btor0=equil_in%global_param%toroid_field%b0
iplasma=equil_in%global_param%i_plasma/1.e6
pressure(1:nrho)=equil_in%profiles_1d%pressure(1:nrho)
pprime(1:nrho)=equil_in%profiles_1d%pprime(1:nrho)
ffprime(1:nrho)=equil_in%profiles_1d%ffprime(1:nrho)
psigrida(1:nrho)=equil_in%profiles_1d%psi(1:nrho) !unnormalized
psigrida(1:nrho)=(psigrida(1:nrho)-psigrida(1))/(psigrida(nrho)-psigrida(1)) !normalized 0 axis, 1 boundary

do i=1,nrho2d
psia_2d(i)=(i-1.)/(nrho2d-1.)
enddo
call linterp(psigrida(1:nrho),ffprime(1:nrho),nrho,psia_2d,ffp_2d,nrho2d)
call linterp(psigrida(1:nrho),pprime(1:nrho),nrho,psia_2d,ppp_2d,nrho2d)
ffp_2d=-GPI2/mu0*ffp_2d
ppp_2d=-GPI2*1.e-6*ppp_2d

ipol(1:nrho)=equil_in%profiles_1d%F_dia(1:nrho)

if (params%k_fixfree.eq.0) then !if 1, comes from free boundary
rexp(1:nteta)=equil_in%eqgeometry%boundary%r(1:nteta)
zexp(1:nteta)=equil_in%eqgeometry%boundary%z(1:nteta)

!define angle not based on mag axis, but on geometrical center
if (j_call.eq.0) then
raxp=raxis_astra
zaxp=zaxis_astra
endif

do i=1,nteta
call find_angle(raxp,zaxp,rexp(i),zexp(i),tetaexp(i))
enddo	

!put points in order
rdum(1:nteta)=rexp(1:nteta)
zdum(1:nteta)=zexp(1:nteta)
tdum(1:nteta)=tetaexp(1:nteta)
do i=1,nteta
j=minloc(tdum(1:nteta),1)
rexp(i)=rdum(j)
zexp(i)=zdum(j)
tetaexp(i)=tdum(j)
tdum(j)=1.e6
enddo

!reorder
if (tetaexp(1).gt.tetaexp(nteta)) then !reorder
j=1
do k=1,nteta-1
if (tetaexp(k+1).lt.tetaexp(k)) j=k+1   ! j is the first teta above 0
enddo		
if (j.gt.1)	tetaexp(1:j-1)=tetaexp(1:j-1)-GPI2
endif

do i=1,nteta
rexp(nteta+i)=rexp(i)
zexp(nteta+i)=zexp(i)
tetaexp(nteta+i)=tetaexp(i)+GPI2
enddo
call linterp(tetaexp(1:nteta*2),rexp(1:nteta*2),nteta*2, &
     &   teta(1:nteta),rbndp(1:nteta),nteta)
call linterp(tetaexp(1:nteta*2),zexp(1:nteta*2),nteta*2, &
     &   teta(1:nteta),zbndp(1:nteta),nteta)

rbnd(1:nteta)=rbndp(1:nteta)
zbnd(1:nteta)=zbndp(1:nteta)

endif
endif

return
end subroutine definitions_ef_equil

!--------------------------------------------------------------------
subroutine equil_ef_init_circ

use pi_vars, only: GPI
use fft_mod_eff, only: sintable, costable
use ef_circuit
use green_matrix
use outcmn_inc, only: machine
use astra2fbe, only: cur_init

implicit none

integer j_files,i,j,k,ii,jj,kk,iii,jjj,kkk,ielem
double precision dummy1,k_fourier
double precision tempcoilr(300),tempcoilz(300), &
     & tempcoilangle(300),tempcoildr(300),tempcoildz(300)
 integer identcoil(5200),tempcoilturns(300),tempcoilelem(300)
integer numeqcump(100),jjelem(300,300),tempnnc(300)
double precision r1,r2,z1,z2,r3,z3,r4,z4,area,gtemp
double precision dr1,dr2,dz1,dz2,dr3,dz3,dr4,dz4,area2,gtemp2
double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,areactmp(300)
double precision rcetmp(5200),zcetmp(5200),datmp(5200)
double precision drcetmp(5200),dzcetmp(5200),areazz(300)
double precision matrix_gs2d(500,500),z_fourier,greenf,ssfw
double precision dummatrix(52,52) !Efable test
integer nctype(5200) !Efable test
double precision dhoriz(200),dvert(200) !blanket elements lengths
integer equivtmp(5200)
integer equivforce(5200)
double precision tatmp(5200)
character(80) fname

fname='exp/cnf/machine_description_out.'//trim(machine)
open(32,file=fname)

read(32,*) nr,nr2,nr1
read(32,*) nz,nz2,nz1
read(32,*) rmin
read(32,*) rmax
read(32,*) zmin
read(32,*) zmax
read(32,*) alpsep

 do i=1,nr2
 	r(i)=rmin+(i-1.)*(rmax-rmin)/nr1     ! computational domain is r(2:nr+1), boundaries are r(1) and r(nr+2)
 enddo
 do i=1,nz2
 	z(i)=zmin+(i-1.)*(zmax-zmin)/nz1
 enddo
rcomp(1:nr)=r(2:nr1)
zcomp(1:nz)=z(2:nz1)
dr=r(2)-r(1)
dz=z(2)-z(1)

do i=1,nz
do j=1,nz
sintable(i,j)=sin(i*j*GPI/(nz+1))
costable(i,j)=cos(i*j*GPI/(nz+1))
enddo
enddo

! load everything from file
  read(32, *) nactive,npassive

read(32, *) ncoils
do i=1, ncoils
    read(32, *) rcoil(i),zcoil(i),drcoil(i),dzcoil(i),anglecoil(i),mequivalence(i)
enddo

read(32, *) nlimiter
do i=1, nlimiter
    read(32, *) limiterr(i),limiterz(i)
enddo
read(32,*) ilim_maxR,lim_maxR
read(32,*) ilim_minR,lim_minR
read(32,*) ilim_maxZ,lim_maxZ
read(32,*) ilim_minZ,lim_minZ

do i=nactive+1,npassive
    read(32, *) r_cond(i),z_cond(i)
enddo

read(32,*) nconduc
do i=1,nconduc	
read(32,*) indconduc(i,1:nconduc)	
enddo
read(32,*) nconduc
do i=1,nconduc	
read(32,*) resconduc(i,1:nconduc)	
enddo
do i=1,nconduc	
do j=1,nr2
read(32,*) greeni(j,1:nz2,i)	
enddo	
enddo
read(32,*) nblocks	
do j=1,nblocks	
do i=1,nblocks
read(32,*) dgreenirj(i,j),dgreenizj(i,j)	
enddo	
enddo
do ii=1,nblocks
do j=1,nz2	
do i=1,nr2
read(32,*) dgreenirpl(i,j,ii),dgreenizpl(i,j,ii)	
enddo	
enddo
enddo
do jj=1,nz2
do ii=1,nr2
read(32,*) zlimpotential(ii,jj)	
enddo	
enddo
read(32,*) ngbnd
read(32,*) green_bnd_f(1:ngbnd)

 close(32)

write(*,*) nactive 

! assign initial currents from astra	
curconduc(1:nconduc)=cur_init(1:nconduc)

return
end subroutine equil_ef_init_circ

!--------------------------------------------------------------------
subroutine fix_boundary_ef(j_init)

use ef_circuit
use astra2fbe
use imas_ids
use metric_coefficients_pbe	
use transfer_functions	
 use pi_vars, only: GPI, GPI2, GPI4, muvac

implicit none

integer j_init,i,j,i_init
type(type_equilibrium) equil_out
double precision psiaxis_new,cnorm,rax_new,zax_new
double precision rhoedge,q_new(nrho),psisave(512,512)
double precision thetap_i(nteta),rmaj2(nrho,nteta), &
 jcbn2(nrho,nteta),darea2(nrho,nteta),effprimp(nrho),epprimp(nrho), &
r_min(nrho,nteta),yy2(nrho,nteta),jrho2(nrho,nteta), gradr2(1:nrho,1:nteta)
integer jr,jt
save psisave

! initial guess
if (j_init.eq.0) then
  raxp=raxis_astra
  zaxp=zaxis_astra
  psiaxisp=psi0_astra
  psibndp=psib_astra	
  do jr=1, Nrho
    do jt=1, Nteta
        lambda2d(jr, jt) = (jr - 1)/(Nr - 1.)
    enddo
  enddo 
  do jt=1, Nteta
  	    psirhoteta(1: Nrho, jt) = psigrida(1: Nrho)
  	enddo
else
   psirhoteta(1:nrho,1:nteta)=psisave(1:nrho,1:nteta)
 psiaxisp=psirhoteta(1,1)
 psibndp=psirhoteta(nrho,1)
endif

!boundary from previous time step

call PHI_EQ_2d_PBE(nrho,nteta,psigrida(1:nrho),iplasma, &
         ffprime(1:nrho),pprime(1:nrho),btor0*rgeom0,&
         rbndp(1:nteta),zbndp(1:nteta),Raxp,Zaxp,psiaxisp, &
         psibndp, solve_fix,j_init, & 
  rpbez(1:nrho,1:nteta), & 
	zpbez(1:nrho,1:nteta), & 
 psirhoteta(1:nrho,1:nteta), & 
   psibez(1:nrho), &   ! psinorm new
 lambda2d(1:nrho,1:nteta), t2dbez(1:nteta), &
  psiaxis_new, cnorm, rax_new, zax_new, thetap_i, rmaj2, & 
	jcbn2, q_new, rhoedge,darea2,epprimp, effprimp,r_min, yy2, gradr2)
	
	raxp=rax_new
	zaxp=zax_new
	psiaxisp=psiaxis_new
	psisave(1:nrho,1:nteta)=psirhoteta(1:nrho,1:nteta)

! notice that qedge is the extrapolate of q_new at neql
ffprimebez(1:nrho)=effprimp(1:nrho)/cnorm
pprimebez(1:nrho)=epprimp(1:nrho)/cnorm
pressbez(1:nrho)=pressure(1:nrho)
ipolbez(1:nrho)=ipol(1:nrho)

!additional calculations

do jt=1,nteta
do jr=1,nrho
 jrhoteta(jr,jt)=( &
      effprimp(jr)/rpbez(jr,jt)+rpbez(jr,jt) &
    *epprimp(jr))/GPI2/0.4/GPI/cnorm
enddo
enddo

do jt=1,nteta
do jr=1,nrho-1
 jrho2(jr,jt)=( &
     & 0.5*(effprimp(jr)+effprimp(jr+1))* &
     &  1./Rmaj2(jr,jt) &
     &  +Rmaj2(jr,jt) &
     &  * &
     & 0.5*(epprimp(jr+1)+epprimp(jr)))/GPI2/0.4/GPI/cnorm
enddo
enddo

jr=nrho-1 !nint(1.*nx)
Z_curr_0D=sum(jrho2(1:jr,:)*YY2(1:jr,:)*darea2(1:jr,:))/ &
     & sum(jrho2(1:jr,:)*darea2(1:jr,:))
R_curr_0D=sqrt(sum(jrho2(1:jr,:)*Rmaj2(1:jr,:)**2.*darea2(1:jr,:))/ &
     & sum(jrho2(1:jr,:)*darea2(1:jr,:)))

!regrid
call build_2dgrid(nrho,nteta,psibez(1:nrho), &
     psirhoteta(1:nrho,1:nteta),rgeom0,pressure(1:nrho), &
   btor0, ipol(1:nrho),iplasma,rpbez(1:nrho,1:nteta),zpbez(1:nrho,1:nteta),&
			rmaj2(1:nrho,1:nteta),r_min(1:nrho,1:nteta),jcbn2(1:nrho,1:nteta), & 
			thetap_i(1:nteta),q_new(1:nrho),rhoedge,gradr2(1:nrho,1:nteta), &
			g2bez(1:nrho), &
      gm1bez(1:nrho), &
      areatbez(1:nrho), &
      perimbez(1:nrho), &
      vbez(1:nrho), &
      g1bez(1:nrho), &
      ggrhobez(1:nrho), &
      bmaxbez(1:nrho), &
      bminbez(1:nrho), &
      gm4bez(1:nrho), &
      bdb0bez(1:nrho), &
      gm5bez(1:nrho), &
      fofbbez(1:nrho), &
      surfbez(1:nrho), & ! lateral surface
			li3,betapol,psplex, &
      bpcell2dbez(1:nrho,1:nteta), &
      bcell2dbez(1:nrho,1:nteta), & 
 					routbez(1:nrho), &
      rinbez(1:nrho), &
      kbez(1:nrho), &
      triaubez(1:nrho), &
      shifbez(1:nrho),gm41bez(1:nrho),qbez(1:nrho)) 

phibez(1:nrho)=0.
 rbp2_b2bez(1:nrho)=0.

rmin2dbez(1:nrho,1:nteta)=r_min(1:nrho,1:nteta)
dator(1:nrho,1:nteta)=darea2(1:nrho,1:nteta)

rpol(1:nrho,1:nteta)=rpbez(1:nrho,1:nteta)
zpol(1:nrho,1:nteta)=zpbez(1:nrho,1:nteta)

rpul(1:nrho,1:nteta)=rpbez(1:nrho,1:nteta)
zpul(1:nrho,1:nteta)=zpbez(1:nrho,1:nteta)

gm41bez(1:nrho)=0.

dpsidvbez(1:nrho)= 0. !to fill 
! additional info from rectangular grid
 
do j=1,nteta
do i=1,nrho
rho(i,j)=sqrt((rpol(i,j)-raxp)**2.+(zpol(i,j)-zaxp)**2.)
enddo
enddo
teta(1:nteta)=t2dbez(1:nteta)
rho(1:nrho,nteta+1)=rho(1:nrho,1)
teta(nteta+1)=teta(1)+GPI2

psia_1d(1:nrho)=psigrida(1:nrho)
ffp_1d(1:nrho)=ffprime(1:nrho)
ppp_1d(1:nrho)=pprime(1:nrho)

return
end subroutine fix_boundary_ef

!--------------------------------------------------------------------
subroutine assignment_of_equilout_stuff(equil_out)

use pi_vars, only: GPI, GPI2
use imas_ids, only: type_equilibrium
use parameters_a2equil
use ef_circuit
use transfer_functions

implicit none

type(type_equilibrium) :: equil_out

equil_out%global_param%psplex=psplex
equil_out%global_param%psibound=-GPI2*psibnd
equil_out%global_param%psiaxis=-GPI2*psiaxis
equil_out%global_param%li3=li3
equil_out%global_param%betpol=betapol	
equil_out%global_param%i_plasma=iplasma*1.e6

equil_out%coord_sys%position%psirz(1:nrho,1:nteta)=psirhoteta(1:nrho,1:nteta)/GPI2
equil_out%coord_sys%position%r(1:nrho,1:nteta)= rpbez(1:nrho,1:nteta)
equil_out%coord_sys%position%z(1:nrho,1:nteta)= zpbez(1:nrho,1:nteta)
equil_out%coord_sys%position%teta2d(1:nteta)=t2dbez(1:nteta)
equil_out%coord_sys%position%rmin(1:nrho,1:nteta)=rmin2dbez(1:nrho,1:nteta)
equil_out%coord_sys%bpcell(1:nrho,1:nteta)=bpcell2dbez(1:nrho,1:nteta)
equil_out%coord_sys%bcell(1:nrho,1:nteta)= bcell2dbez(1:nrho,1:nteta)

equil_out%eqgeometry%rectgrid%npointsr=nr2
equil_out%eqgeometry%rectgrid%npointsz=nz2
equil_out%eqgeometry%rectgrid%r2d(1:nr2)=r(1:nr2)
equil_out%eqgeometry%rectgrid%z2d(1:nz2)=z(1:nz2)	
equil_out%eqgeometry%rectgrid%psirz2d(1:nr2,1:nz2)=psirz(1:nr2,1:nz2)
equil_out%eqgeometry%rectgrid%psi_axis=psiaxis
equil_out%eqgeometry%rectgrid%psi_boundary=psibnd

equil_out%profiles_1d%ffprime(1:nrho)=ffprimebez(1:nrho)
equil_out%profiles_1d%pprime(1:nrho)=pprimebez(1:nrho)
equil_out%profiles_1d%pressure(1:nrho)=pressbez(1:nrho)
equil_out%profiles_1d%rho_tor(1:nrho)=0.
equil_out%profiles_1d%F_dia(1:nrho)=ipolbez(1:nrho)

equil_out%profiles_1d%dPSIdV(1:nrho)= dpsidvbez(1:nrho)
equil_out%profiles_1d%psi(1:nrho)=  psibez(1:nrho) ! psi nnormalized from 0 to 1 (rhopol^2)
equil_out%profiles_1d%g2(1:nrho)=g2bez(1:nrho)
equil_out%profiles_1d%g2int(1:nrho)=g2ibez(1:nrho)
equil_out%profiles_1d%gm1(1:nrho)=gm1bez(1:nrho)
equil_out%profiles_1d%r_outboard(1:nrho)=routbez(1:nrho)
equil_out%profiles_1d%r_inboard(1:nrho)=rinbez(1:nrho)
equil_out%profiles_1d%volume(1:nrho)=vbez(1:nrho)
equil_out%profiles_1d%g1(1:nrho)=g1bez(1:nrho)
equil_out%profiles_1d%gm41(1:nrho)=gm41bez(1:nrho)
equil_out%profiles_1d%ggradro(1:nrho)=ggrhobez(1:nrho)
equil_out%profiles_1d%bmaxt(1:nrho)=bmaxbez(1:nrho)
equil_out%profiles_1d%bmint(1:nrho)=bminbez(1:nrho)
equil_out%profiles_1d%gm4(1:nrho)=gm4bez(1:nrho)
equil_out%profiles_1d%bdb0(1:nrho)=bdb0bez(1:nrho)
equil_out%profiles_1d%gm5(1:nrho)=gm5bez(1:nrho)
equil_out%profiles_1d%fofb(1:nrho)=fofbbez(1:nrho)
equil_out%profiles_1d%areat(1:nrho)=areatbez(1:nrho)
equil_out%profiles_1d%perim(1:nrho)=perimbez(1:nrho)
equil_out%profiles_1d%shif(1:nrho)=shifbez(1:nrho)
equil_out%profiles_1d%elongation(1:nrho)=kbez(1:nrho)
equil_out%profiles_1d%surface(1:nrho)= surfbez(1:nrho) ! lateral surface
equil_out%profiles_1d%tria_upper(1:nrho)=triaubez(1:nrho)
equil_out%profiles_1d%tria_lower=equil_out%profiles_1d%tria_upper(1:nrho)
equil_out%profiles_1d%phi(1:nrho)=phibez(1:nrho)
equil_out%profiles_1d%q(1:nrho)=qbez(1:nrho)
equil_out%profiles_1d%rbp_b2(1:nrho)=rbp2_b2bez(1:nrho)

return
end subroutine assignment_of_equilout_stuff

!--------------------------------------------------------------------
subroutine convert_boundary_to_pbe

use pi_vars, only: GPI2
use ef_circuit
use feqis_tools, only: find_angle, find_fields_interp_psionly

implicit none

double precision teta_fbe(i_dim5)
double precision x1,x2,x3,x4,t1,t2,t3,t4,z1,z2,z3,z4
double precision x11,x22,x33,x44,t11,t22,t33,t44,z11,z22,z33,z44
double precision dx
integer i,j,k,i1,i2,i3,i4,j1,j2,j3,j4

do i=1,nteta+1
teta(i)=GPI2*(i-1.)/(nteta+0.)
enddo
dteta=teta(2)-teta(1)

dx =sqrt(dr**2.+dz**2.)

! find boundary
j=jaxis
do i=iaxis,nr2
if (psirz(i,j).ge.psibnd) k=i
if (psirz(i,j).le.psibnd) goto 10
enddo	

10 continue

if (psirz(k,j).eq.psibnd) then
rbnd(1)=r(k)
else
rbnd(1)=r(k)-(psirz(k,j)-psibnd)/(psirz(k,j)-psirz(k-1,j))*dr		
endif	
zbnd(1)=z(j)
call find_angle(rax,zax,rbnd(1),zbnd(1),teta_fbe(1))	
i=1

do i=2,nteta

dx =sqrt(dr**2.+dz**2.)

x1=sqrt((rbnd(i-1)-rax)**2.+(zbnd(i-1)-zax)**2.)
teta_fbe(i)=teta_fbe(i-1)+dteta

t1=rax+x1*cos(teta_fbe(i))
t2=zax+x1*sin(teta_fbe(i))

if (t1.ge.raus) then
x1=(raus-rax)/cos(teta_fbe(i))
endif
if (t1.le.rinner) then
x1=(rinner-rax)/cos(teta_fbe(i))
endif
if (t2.ge.ztop) then
x1=(ztop-zax)/sin(teta_fbe(i))
endif
if (t2.le.zbot) then
x1=(zbot-zax)/sin(teta_fbe(i))
endif
t1=rax+x1*cos(teta_fbe(i))
t2=zax+x1*sin(teta_fbe(i))
call find_fields_interp_psionly(t1,t2,t3) !give back psi,br,bz at r0,z0
if (t3.eq.psibnd) then
rbnd(i)=t1
zbnd(i)=t2
goto 2
endif

if (t3.lt.psibnd) then
4 continue
z1=rax+(x1-dx)*cos(teta_fbe(i))
z2=zax+(x1-dx)*sin(teta_fbe(i))
call find_fields_interp_psionly(z1,z2,z3) 
if (z3.lt.psibnd) then
	dx=1.1*dx
	goto 4
endif
x11=x1-(t3-psibnd)/(t3-z3)*dx		
rbnd(i)=rax+x11*cos(teta_fbe(i))
zbnd(i)=zax+x11*sin(teta_fbe(i))
goto 2
endif

if (t3.gt.psibnd) then
3	continue
j4=0
x2=x1+dx
z1=rax+x2*cos(teta_fbe(i))
z2=zax+x2*sin(teta_fbe(i))
if (z1.ge.raus) then
x2=(raus-rax)/cos(teta_fbe(i))
j4=1
endif
if (z1.le.rinner) then
x2=(rinner-rax)/cos(teta_fbe(i))
j4=1
endif
if (z2.ge.ztop) then
x2=(ztop-zax)/sin(teta_fbe(i))
j4=1
endif
if (z2.le.zbot) then
x2=(zbot-zax)/sin(teta_fbe(i))
j4=1
endif
z1=rax+x2*cos(teta_fbe(i))
z2=zax+x2*sin(teta_fbe(i))
call find_fields_interp_psionly(z1,z2,z3) 
if (z3.gt.psibnd) then
	if (j4.eq.1) then
		rbnd(i)=z1
		zbnd(i)=z2
		goto 2
	endif
	dx=dx*1.1
	goto 3
endif
x11=x1+(t3-psibnd)/(t3-z3)*dx		
rbnd(i)=rax+x11*cos(teta_fbe(i))
zbnd(i)=zax+x11*sin(teta_fbe(i))
goto 2
endif

2	continue

enddo		

nbnd=nteta

teta(1:nteta)=teta_fbe(1:nteta)
teta(nteta+1)=teta(nteta)+dteta
rbndp(1:nteta)=rbnd(1:nteta)
zbndp(1:nteta)=zbnd(1:nteta)
rbndp(nteta+1)=rbnd(nteta) !shouldnt this be rbnd(1) ???
zbndp(nteta+1)=zbnd(nteta)

raxp=rax
zaxp=zax
psibndp=psibnd
psiaxisp=psiaxis

return
end subroutine convert_boundary_to_pbe

!--------------------------------------------------------------------
subroutine solve_gs2d(g)

use ef_circuit
use fft_mod_eff, only: costable
use feqis_tools, only: discrete_sine_transform

implicit none

double precision g(i_dim2,i_dim2)
double precision gt(i_dim2,i_dim2)
double precision rhs(i_dim2,i_dim2)
double precision wrhs(i_dim2,i_dim2)
double precision trhs(i_dim2,i_dim2)
double precision A(258),z_fourier(258)
double precision B(258),x1,x2
double precision C(258),tin,tup,r1m_1,r2m_1
integer i,j,k,k_fourier
integer j_init

data j_init/0/
save A,B,C,j_init,z_fourier

do j=1,nz2
do i=1,nr2
rhs(i,j)=-mu0*r(i)*jrz(i,j)
enddo
enddo

r1m_1=r(2)/dr**2./((r(1)+r(2))/2.)
r2m_1=r(nr1)/dr**2./((r(nr1)+r(nr2))/2.)

rhs(2:nr1,2)=rhs(2:nr1,2)-g(2:nr1,1)/dz**2
rhs(2:nr1,nz1)=rhs(2:nr1,nz1)-g(2:nr1,nz2)/dz**2

rhs(2,2:nz1)=rhs(2,2:nz1)-g(1,2:nz1)*r1m_1
rhs(nr1,2:nz1)=rhs(nr1,2:nz1)-g(nr2,2:nz1)*r2m_1

wrhs=rhs
!	CALL CPU_TIME(tin)
do i=2,nr1
  	call discrete_sine_transform(nz,wrhs(i,2:nz1))
enddo

k_fourier = nz
!create inverse matrix for gs2d
if (j_init.eq.0) then
A=0.
B=0.
C=0.
do i=2,nr-1		
x1=0.5*(rcomp(i)+rcomp(i-1))		
x2=0.5*(rcomp(i+1)+rcomp(i))		
B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
A(i)=rcomp(i)/x2/dr**2. 
C(i)=rcomp(i)/x1/dr**2. 
enddo
i=1
x1=0.5*(rcomp(i)+r(1))		
x2=0.5*(rcomp(i+1)+rcomp(i))		
B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
A(i)=rcomp(i)/x2/dr**2. 
i=nr
x1=0.5*(rcomp(i)+rcomp(i-1))		
x2=0.5*(r(nr2)+rcomp(i))		
B(i)=-rcomp(i)/dr**2.*(1./x2+1./x1) 
C(i)=rcomp(i)/x1/dr**2. 
j_init=1
do k=2,nz1
z_fourier(k)=2./dz**2.*(costable(1,k-1)-1.)
enddo
endif

gt=0.
!solve matrix
do k=2,nz1
call solve_tridiag_fbe_ef(C(1:nr),B(1:nr)+z_fourier(k),A(1:nr),wrhs(2:nr1,k),gt(2:nr1,k),nr)
enddo
!invert fourier from gt(1:nr,1:kfourier) to g(2:nr1,2:nz1)
!			gt(i,k)=sum(invMM_gs2d(i-1,1:nr,k-1)*wrhs(2:nr1,k))

do i=2,nr1
  	call discrete_sine_transform(nz,gt(i,2:nz1))
enddo
g(2:nr1,2:nz1)=2./(nz+1)*gt(2:nr1,2:nz1)

return
end subroutine solve_gs2d

!--------------------------------------------------------------------
subroutine solve_tridiag_fbe_ef(A,B,C,R,f,Ngrid)

!C Provides solution of the system:
!C
!C   Aj fj-1  + Bj fj  + Cj fj+1 = Rj
!C
!C   where j = 1..Ngrid
!C
!C   bcbound = 1  -> given f_NA1
!C
!C  eximp = 2: implicit
!C
!C by means of this method:

implicit none

integer :: i,j,k,Ngrid,NgridS,bcbound
double precision A(Ngrid),B(Ngrid)
double precision C(Ngrid),R(Ngrid)
double precision f(Ngrid),f_bound
double precision alpha(Ngrid),beta(Ngrid),f0
double precision Bstar,Cstar,Rstar
integer eximp

alpha=0.
beta=0.

         j = 1
   alpha(j) = C(j)/B(j)
         beta(j) = R(j)/B(j)
 
   do j=2,Ngrid-1
    alpha(j) = C(j)/(B(j)-A(j)*alpha(j-1))
   enddo
   do j=2,Ngrid
    beta(j) = (R(j)-A(j)*beta(j-1))/(B(j)-A(j)*alpha(j-1))
   enddo

! Note that boundary value is assumed to be on the main last grid point, so 1-dx/2. This has to be
! corrected later on... CEfable
   j = Ngrid

!C			     f(j-1)=(f_bound-beta(j-1))/alpha(j-1)
!CTHis one should be appropriate with extrapolation... but now go back to real b.c.
!C			     f(j-1)=(2./3.*f_bound-beta(j-1))/(alpha(j-1)-1./3.)
     f(j)=beta(j)
		 do k=1,Ngrid-1
		  j=Ngrid-1-k+1
		  f(j)=beta(j)-alpha(j)*f(j+1)
		 enddo

return
end subroutine solve_tridiag_fbe_ef

!--------------------------------------------------------------------
subroutine find_new_axis_part1	

use ef_circuit
use feqis_tools, only: nine_point_regression

implicit none

integer jaold,iaold,i_mode
integer i,j,k,iax,jax,i_ixpoint
double precision errtol,tolerr,raxm,zaxm
double precision br1,bz1,dbrr,dbrz,dbzr,dbzz
double precision br2,bz2,darax,dazax,bub(100),xub(100),yub(100)
double precision br3,bz3,icase,psi0,ccc(6)
double precision br4,bz4,determ,rleft,rright
double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,psibtmp,x22,x33,x10,x11
double precision ppx(2)
integer i1,i2,i3,i4,i5,i6,i7,i8,i9
integer j1,j2,j3,j4,j5,j6,j7,j8,j9,i0,j0
!the used function is psirz

! 1) find new magnetic axis
!old axis
!	call get_closest_index_ef(rax,zax,iax,jax)
iax=iaxis
jax=jaxis
i1=10000
j1=10000
i2=10000
j2=10000
errtol=1.e6
tolerr=10.
j=1
raxm=rax
zaxm=zax
!	write(44125,*) ' '
  do while(errtol.gt.tolerr)
i1=i2
j1=j2
i2=iax
j2=jax
if (psirz(iax+1,jax).gt.psirz(iax,jax)) then
iax=iax+1
jax=jax
endif
if (psirz(iax-1,jax).gt.psirz(iax,jax)) then
iax=iax-1
jax=jax
endif
if (psirz(iax,jax+1).gt.psirz(iax,jax)) then
iax=iax
jax=jax+1
endif
if (psirz(iax,jax-1).gt.psirz(iax,jax)) then
iax=iax
jax=jax-1
endif
if (psirz(iax+1,jax-1).gt.psirz(iax,jax)) then
iax=iax+1
jax=jax-1
endif
if (psirz(iax-1,jax+1).gt.psirz(iax,jax)) then
iax=iax-1
jax=jax+1
endif
if (psirz(iax-1,jax-1).gt.psirz(iax,jax)) then
iax=iax-1
jax=jax-1
endif
if (psirz(iax+1,jax+1).gt.psirz(iax,jax)) then
iax=iax+1
jax=jax+1
endif

j=j+1

if ((iax.eq.i1).and.(jax.eq.j1)) goto 101
if (j.ge.100000) goto 101		
enddo
101 continue	

iaxis=iax
jaxis=jax	
write(*,*) 'ax',rax,zax,iax,jax,r(iax),z(jax),psirz(iax,jax)

call nine_point_regression(r(iax),z(jax),ppx,derivpsi,psiaxis)

rax=ppx(1)  !r(iaxis)
zax=ppx(2) !z(jaxis)
trax=rax
tzax=zax

write(*,*) 'ax2',rax,zax

if (isnan(rax)) then
write(*,*) 'rax is nan in fbe find axis'
stop
endif
if (rax.lt.0.1.or.rax.gt.100) then
write(*,*) 'rax is nan in fbe find axis',rax,zax,ppx,derivpsi,psiaxis
stop
endif

return
end subroutine find_new_axis_part1

!--------------------------------------------------------------------
subroutine find_psi_boundary

use pi_vars, only: GPI
use ef_circuit
use astra2fbe, only: x_point_save, plasma_config
use errors_params, only: err_find_oxpoints_derivs
use feqis_tools, only: find_closest_xpoints, find_fields_interp_psionly, &
    get_closest_index, find_angle, check_xpoint_connection_axis, &
    nine_point_regression, nine_point_regression_follow, t_find_u_n

implicit none

integer jaold,iaold,i_mode,niter
integer i,j,k,iax,jax,i_ixpoint
double precision errtol,tolerr,raxm,zaxm
double precision br1,bz1,dbrr,dbrz,dbzr,dbzz,rstart
double precision br2,bz2,darax,dazax,bub(9),xub(9),yub(9)
double precision br3,bz3,icase,psi0,ccc(6),psiloc,polfield
double precision br4,bz4,determ,rleft,rright,omega_temp(300,300)
double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,psibtmp,x22,x33,x10,x11
double precision xbnd(i_dim5),ybnd(i_dim5),r_temp(i_dim2),z_temp(i_dim2),rminz
double precision pos_xpoint(2),ddipsi(8),tin,tup,frlim_ef
double precision r_limp,z_limp,psi_limp(500),hard_left,hard_right
double precision psi_xpoint(max_xpoints),ddpsi(5)
double precision rx_add(20),zx_add(20)
 	integer ipath(5),jpath(5),oldpointnum
integer i1,i2,i3,i4,i5,i6,i7,i8,i9,n_adding
integer j1,j2,j3,j4,j5,j6,j7,j8,j9,i0,j0,i_county
data i_county/0/
save i_county,rstart,oldpointnum

i_plasmatype=0

if (n_of_xpoints.ge.1) then
x_point_save(1:n_of_xpoints,1)=r_xpoint(1:n_of_xpoints)
x_point_save(1:n_of_xpoints,2)=z_xpoint(1:n_of_xpoints)
endif

call find_closest_xpoints(rx_add,zx_add,i9,n_adding)
if (i9.eq.0) then
r_xpoint(n_of_xpoints+1:n_of_xpoints+n_adding)=rx_add(1:n_adding)
z_xpoint(n_of_xpoints+1:n_of_xpoints+n_adding)=zx_add(1:n_adding)
n_of_xpoints=n_of_xpoints+n_adding
endif

if (i_county.eq.1) then
!	go through old x-points and see where they end up
if (n_of_xpoints.ge.1) then
iaold=n_of_xpoints
do i=1,iaold

	niter=0

289	niter=niter+1				
	call nine_point_regression_follow(r_xpoint(i),z_xpoint(i),pos_xpoint,ddipsi,x1)

	r_xpoint(i)=pos_xpoint(1)
	z_xpoint(i)=pos_xpoint(2)
	if ( (abs(ddipsi(1))+abs(ddipsi(2))).le.err_find_oxpoints_derivs) goto 290
	if (niter.gt.100) goto 291 !no point found
	if ((pos_xpoint(1).gt.r(nr2)-dr).or.(pos_xpoint(1).lt.r(1)+dr).or. & 
	 & (pos_xpoint(2).gt.z(nz2)-dz).or.(pos_xpoint(2).lt.z(1)+dz)) then
! xpoint doesnt exist anymore
		r_xpoint(i)=1.e6
		z_xpoint(i)=0.
		deriv_x(1:5,i)=1.e6
		goto 317
	endif
	goto 289
291 continue
		r_xpoint(i)=1.e6
		z_xpoint(i)=0.
		deriv_x(1:5,i)=1.e6
		goto 317
290	continue

!check if point outside of domain
	if ((pos_xpoint(1).gt.r(nr2)-dr).or.(pos_xpoint(1).lt.r(1)+dr).or. & 
	 & (pos_xpoint(2).gt.z(nz2)-dz).or.(pos_xpoint(2).lt.z(1)+dz)) then
! xpoint doesnt exist anymore
		r_xpoint(i)=1.e6
		z_xpoint(i)=0.
		deriv_x(1:5,i)=1.e6
		goto 317
	endif

!check if point is on axis
	if ((pos_xpoint(1).gt.rax-dr).and.(pos_xpoint(1).lt.rax+dr).and. & 
	 & (pos_xpoint(2).gt.zax-dz).and.(pos_xpoint(2).lt.zax+dz)) then
! xpoint doesnt exist anymore

		r_xpoint(i)=1.e6
		z_xpoint(i)=0.
		deriv_x(1:5,i)=1.e6
		goto 317
	endif

	r_xpoint(i)=pos_xpoint(1)
	z_xpoint(i)=pos_xpoint(2)
	deriv_x(1:5,i)=ddpsi(1:5)

317	continue
enddo			
endif

!scan the boundary to find new x-points
j=2
do i=2,nr1
call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)

				x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))

if ((pos_xpoint(1).ge.r(i)-dr).and. & 
				& (pos_xpoint(1).le.r(i)+dr).and. & 
	& (pos_xpoint(2).ge.z(j)-dz).and. & 
	& (pos_xpoint(2).le.z(j)+dz).and. & 
	& (x5.ge.0.)) then

n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
r_xpoint(n_of_xpoints)=pos_xpoint(1)
z_xpoint(n_of_xpoints)=pos_xpoint(2)
deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)

endif
enddo
j=nz1
do i=2,nr1
call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)

				x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))

if ((pos_xpoint(1).ge.r(i)-dr).and. & 
				& (pos_xpoint(1).le.r(i)+dr).and. & 
	& (pos_xpoint(2).ge.z(j)-dz).and. & 
	& (pos_xpoint(2).le.z(j)+dz).and. & 
	& (x5.ge.0.)) then

n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
r_xpoint(n_of_xpoints)=pos_xpoint(1)
z_xpoint(n_of_xpoints)=pos_xpoint(2)
deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)

endif
enddo
i=2
do j=2,nz1
call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)

				x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))

if ((pos_xpoint(1).ge.r(i)-dr).and. & 
				& (pos_xpoint(1).le.r(i)+dr).and. & 
	& (pos_xpoint(2).ge.z(j)-dz).and. & 
	& (pos_xpoint(2).le.z(j)+dz).and. & 
	& (x5.ge.0.)) then

n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
r_xpoint(n_of_xpoints)=pos_xpoint(1)
z_xpoint(n_of_xpoints)=pos_xpoint(2)
deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)

endif
enddo
i=nr1
do j=2,nz1
call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)

				x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))

if ((pos_xpoint(1).ge.r(i)-dr).and. & 
				& (pos_xpoint(1).le.r(i)+dr).and. & 
	& (pos_xpoint(2).ge.z(j)-dz).and. & 
	& (pos_xpoint(2).le.z(j)+dz).and. & 
	& (x5.ge.0.)) then

n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
r_xpoint(n_of_xpoints)=pos_xpoint(1)
z_xpoint(n_of_xpoints)=pos_xpoint(2)
deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)

endif
enddo

endif

!-------------------
if (i_county.eq.0) then ! do a full pass to find all X-points
n_of_xpoints=0
do j=2,nz1
do i=2,nr1
	call nine_point_regression(r(i),z(j),pos_xpoint,ddipsi,x1)
				x5=(ddipsi(5)**2.-ddipsi(3)*ddipsi(4))

				if ((pos_xpoint(1).ge.r(i)-dr).and. & 
				& (pos_xpoint(1).le.r(i)+dr).and. & 
	& (pos_xpoint(2).ge.z(j)-dz).and. & 
	& (pos_xpoint(2).le.z(j)+dz).and. & 
	& (x5.ge.0.)) then

					if (abs(pos_xpoint(1)-rax).le.2.*dr.and.abs(pos_xpoint(2)-zax).le.2.*dz) then
						else
						n_of_xpoints=min(max_xpoints,n_of_xpoints+1)
						r_xpoint(n_of_xpoints)=pos_xpoint(1)
						z_xpoint(n_of_xpoints)=pos_xpoint(2)
						deriv_x(1:5,n_of_xpoints)=ddpsi(1:5)

					endif

				endif
enddo
enddo
i_county=1
oldpointnum=n_of_xpoints
endif

!remove disappeared x-points and double counts
318	i=0

319	i=i+1	
!	write(*,*) 'i',i,n_of_xpoints,r_xpoint(i)
if (i.gt.n_of_xpoints) goto 320
if (r_xpoint(i).ge.1.e5) then
		if (n_of_xpoints.eq.1) then
			n_of_xpoints=0
			goto 320
		endif
!					write(*,*) 'big',i
		r_xpoint(i:n_of_xpoints-1)=r_xpoint(i+1:n_of_xpoints)
		z_xpoint(i:n_of_xpoints-1)=z_xpoint(i+1:n_of_xpoints)
		n_of_xpoints=n_of_xpoints-1
		i=i-1
		goto 319
endif

do k=1,i-1 ! check if double counted
if ((abs(r_xpoint(i)-r_xpoint(k)).le.2.*dr).and. & 
 & (abs(z_xpoint(i)-z_xpoint(k)).le.2.*dz)) then					
!		 write(*,*) 'dupli',i
r_xpoint(k)=0.5*(r_xpoint(i)+r_xpoint(k))	
z_xpoint(k)=0.5*(z_xpoint(i)+z_xpoint(k))	
deriv_x(1:5,k)=0.5*(deriv_x(1:5,i)+deriv_x(1:5,k))
		r_xpoint(i:n_of_xpoints-1)=r_xpoint(i+1:n_of_xpoints)
		z_xpoint(i:n_of_xpoints-1)=z_xpoint(i+1:n_of_xpoints)
		n_of_xpoints=n_of_xpoints-1
		i=i-1
goto 319		
endif
enddo

goto 319

320	continue

write(*,*) 'new points ',n_of_xpoints,r_xpoint(1:n_of_xpoints),z_xpoint(1:n_of_xpoints)
write(*,*) 'ax',rax,zax

oldpointnum=n_of_xpoints

!ignore limiter if use_limiter_astra is 0, da trasferirsi in init
if (use_limiter_yesno.eq.0) then
limiterR=r(nr1)
limiterZ=z(nz1)	
endif	
!calculate limiter flux 
do i=1,nlimiter
call find_fields_interp_psionly(limiterR(i),limiterZ(i),psi_limp(i)) !give back psi,br,bz at r0,z0
enddo

if (n_of_xpoints.eq.0) then
   !no x-points, take largest limiter flux
   psibnd = maxval(psi_limp(1:nlimiter),1)
   i_plasmatype=0
endif

! now, remove limiters that are in the shadow of xpoints
ztop=1.e6
zbot=-1.e6
raus=1.e6
rinner=0.
if (n_of_xpoints.ge.1) then
   i_plasmatype=1

!first pass, remove X-points behind the limiter area
   do i=1,n_of_xpoints
      call get_closest_index(r_xpoint(i),z_xpoint(i),j,k)
if (zlimpotential(j,k).lt.0.5) then
         psi_xpoint(i)=-1.e6
  else
         call find_fields_interp_psionly(r_xpoint(i),z_xpoint(i),psi_xpoint(i))
endif
   enddo

!second pass, remove limiter points that are in x-points shadow, simple "straight line method" --> to be refined later on
   do i=1,n_of_xpoints
    if (psi_xpoint(i).gt.-1.e5) then
	call find_angle(rax,zax,r_xpoint(i),z_xpoint(i),x1)
	if ((r_xpoint(i).gt.rax).and.(x1.ge.7./4.*GPI.or.x1.le.GPI/4.)) raus=min(raus,r_xpoint(i))
	if ((z_xpoint(i).gt.zax).and.(x1.ge.GPI/4..and.x1.le.3./4.*GPI)) ztop=min(ztop,z_xpoint(i))
	if ((r_xpoint(i).lt.rax).and.(x1.ge.3./4.*GPI.and.x1.le.5./4.*GPI)) rinner=max(rinner,r_xpoint(i))
	if ((z_xpoint(i).lt.zax).and.(x1.ge.5./4.*GPI.and.x1.le.7./4.*GPI)) zbot=max(zbot,z_xpoint(i))
    endif
   enddo
   do j=1,nlimiter
      if (limiterr(j).le.rinner) psi_limp(j)=-1.e6
  if (limiterr(j).ge.raus) psi_limp(j)=-1.e6
if (limiterz(j).le.zbot) psi_limp(j)=-1.e6
if (limiterz(j).ge.ztop) psi_limp(j)=-1.e6
   enddo			

!third pass, remove x-points that are non-monotonically connected to the plasma.
   do i=1,n_of_xpoints
    if (psi_xpoint(i).gt.-1.e5) then
       call check_xpoint_connection_axis(r_xpoint(i),z_xpoint(i),rax,zax,dr,dz,i1)
       if (i1.eq.0) psi_xpoint(i)=-1.e6
    endif
   enddo

!psi boundary is the highest of all the values
x1=maxval(psi_xpoint(1:n_of_xpoints),1)
x2=maxval(psi_limp(1:nlimiter),1)
i4=maxloc(psi_xpoint(1:n_of_xpoints),1)
i5=maxloc(psi_limp(1:nlimiter),1)
psibnd=max(x1,x2)
if (x2.gt.x1) i_plasmatype=0
if (x1.ge.x2) i_plasmatype=1
endif

write(*,*) 'plasma type',psibnd,i_plasmatype,x1,x2,psi_xpoint(i4),psi_limp(i5),limiterr(i5),limiterz(i5),r_xpoint(i4),z_xpoint(i4)

plasma_config = i_plasmatype
x_point_save(20,1)=r_xpoint(i4)
x_point_save(20,2)=z_xpoint(i4)

!	if (i_plasmatype.eq.1) 
psibnd=psiaxis+(psibnd-psiaxis)*alpsep

!normalized flux
u_n(1:nr2,1:nz2)=(psirz(1:nr2,1:nz2)-psiaxis)/(psibnd-psiaxis)

return
end subroutine find_psi_boundary

!--------------------------------------------------------------------
subroutine new_jrz_ef ! calculate new right hand side given new boundary!

use ef_circuit
use astra2fbe
use feqis_tools, only: t_find_u_n

implicit none

integer i,j,k,i1,i2,i3,i4,i5,j1,j2,j3,j4,j5
double precision dum1,dum2,dum3,zeta,dumc(i_dim2,i_dim2)
double precision t1,t2,t3,t4,x,y,alp,bet,gam,det,det0
double precision je1,je2,je3,je4
integer quadrant
double precision rpluz,zpluz
double precision z11,z12,z13,z14
integer ipluz,jpluz,qipluz,qjpluz
integer ilast,jlast,totpoints,istart,jcallaz
integer external_griddo_j(90000,2),j_griddo_j
integer internal_griddo(90000,2),i_griddo_j

data jcallaz/0/
save jcallaz

! in entry: rbnd,zbnd,nbnd,psiaxis,psibnd,u_n

dumc=0.	
area_eff=dr*dz
i_griddo_j=0
j_griddo_j=0

totpoints=0
do quadrant=1,4
! sweep from axis to exterior and fill in the current
!start from axis position
i1=floor((rax-rmin)/dr+1.)	
j1=floor((zax-zmin)/dz+1.)
!write(6712,*)  ' quadrant',quadrant
if (quadrant.eq.1) then
i=i1+1
j=j1+1
ipluz=1
jpluz=1
endif
if (quadrant.eq.2) then
i=i1
j=j1+1
ipluz=-1
jpluz=1
endif
if (quadrant.eq.3) then
i=i1
j=j1
ipluz=-1
jpluz=-1
endif
if (quadrant.eq.4) then
i=i1+1
j=j1
ipluz=1
jpluz=-1
endif
istart=i

558 continue
totpoints=totpoints+1
if (totpoints.gt.nz2*nr2) then
write(*,*) 'error in find new boundary (totpoints.gt.nz2*nr2)'
stop
endif

call fill_in_current(r(i),z(j), & 
		& nrho2d,psia_2d,ppp_2d,ffp_2d,dumc(i,j),u_n(i,j))
i_griddo_j=i_griddo_j+1
internal_griddo(i_griddo_j,1)=i
internal_griddo(i_griddo_j,2)=j
i1=i
i2=j

call t_find_u_n(i1+ipluz,i2,i1,i2,t1)		
if (t1.lt.0.or.t1.gt.1.) t1=1.e6
call t_find_u_n(i1,i2+jpluz,i1,i2,t2)		
if (t2.lt.0.or.t2.gt.1.) t2=1.e6
if (t2.gt.0..and.t2.le.1.) then
j_griddo_j=j_griddo_j+1
external_griddo_j(j_griddo_j,1)=i
external_griddo_j(j_griddo_j,2)=j+jpluz
endif

if (t1.gt.1.e5) then
!move horizontally to the right
i=i+ipluz
goto 558
endif

ilast=i+ipluz

j_griddo_j=j_griddo_j+1
external_griddo_j(j_griddo_j,1)=ilast
external_griddo_j(j_griddo_j,2)=j

!found boundary, go back, check vertically

i=istart
59781	continue
call t_find_u_n(i,j+jpluz,i,j,t2)		
if (t2.lt.0.or.t2.gt.1.) then
goto 5581
endif

goto 5582

5581	continue
j=j+jpluz
goto 558

!found boundary on Z, need to advance 1 more
5582	continue

j_griddo_j=j_griddo_j+1
external_griddo_j(j_griddo_j,1)=i
external_griddo_j(j_griddo_j,2)=j+jpluz

i=i+ipluz
istart=i
if (i.eq.ilast) then
!terminated quadrant
goto 559
endif
goto 59781

559 continue

enddo !quadrant cycle

!finished quadrants, now do trick at boundary for ciurrent

! fill current in external griddo
do i=1,j_griddo_j
t1=0.
t2=0.
t3=0.
t4=0.
i1=external_griddo_j(i,1)
i2=external_griddo_j(i,2)
call t_find_u_n(i1-1,i2,i1,i2,t1)		
if (t1.lt.0.or.t1.gt.1.) t1=0.
call t_find_u_n(i1,i2+1,i1,i2,t2)		
if (t2.lt.0.or.t2.gt.1.) t2=0.
call t_find_u_n(i1+1,i2,i1,i2,t3)		
if (t3.lt.0.or.t3.gt.1.) t3=0.
call t_find_u_n(i1,i2-1,i1,i2,t4)		
if (t4.lt.0.or.t4.gt.1.) t4=0.

je1=0.
je2=0.
je3=0.
je4=0.
z11=r(i1-1)*(1-t1)+r(i1)*t1
z12=r(i1)
z13=r(i1+1)*(1-t3)+r(i1)*t3
z14=r(i1)
if (t1.gt.0.)	call fill_in_current(z11,z(i2), & 
		& nrho2d,psia_2d,ppp_2d,ffp_2d,je1,1.d0)
if (t2.gt.0.)		call fill_in_current(z12,z(i2+1), & 
		& nrho2d,psia_2d,ppp_2d,ffp_2d,je2,1.d0)
if (t3.gt.0.)	call fill_in_current(z13,z(i2), & 
		& nrho2d,psia_2d,ppp_2d,ffp_2d,je3,1.d0)
if (t4.gt.0.)	call fill_in_current(z14,z(i2-1), & 
		& nrho2d,psia_2d,ppp_2d,ffp_2d,je4,1.d0)

! defining S1 = dR - d1, S2 = dZ-d2, S3 = dZ-d3, S4 = dR-d4
! Jvacuum = C*Jb
! if 1 point only: C = 1 - S1/dR
! if 2 points: C = 1 - S1*S2/(dR*dZ)
! if 3 points: C = 1 - (S1*S2 + S1*S3)/(2*dR*dZ)
! if 4 points: C = 1 - (S1*S2 + S1*S3 + S3*S4 + S4*S2)/(4*dR*dZ)
! t_j = d_j /(dR or dZ depending on direction)

dumc(i1,i2)=t1*je1+t2*je2+t3*je3+t4*je4- & 
 (t1*t2*(je1+je2)/2.+ & 
 t1*t3*(je1+je3)/2.+ & 
 t1*t4*(je1+je4)/2.+ & 
 t2*t3*(je3+je2)/2.+ & 
 t2*t4*(je4+je2)/2.+ & 
 t3*t4*(je3+je4)/2.)

enddo

jrz(1:nr2,1:nz2)=dumc(1:nr2,1:nz2)

!rescale current
dum1=sum(jrz*area_eff) 
jrz=jrz/dum1*iplasma

if (isnan(dum1)) then
write(*,*) 'total current is nan in fbe current rescaling'
stop
endif

t1=0.
t2=0.
t3=0.
do j=1,nz2
do i=1,nr1
t1=t1+r(i)**2.*jrz(i,j)*area_eff(i,j)
t2=t2+z(j)*jrz(i,j)*area_eff(i,j)
t3=t3+jrz(i,j)*area_eff(i,j)
enddo
enddo

R_curr_2D= sqrt(t1/t3)
Z_curr_2D=	t2/t3

return
end subroutine new_jrz_ef

!--------------------------------------------------------------------
subroutine fill_in_current(r0,z0,nx,psia_2d,ppp_2d,ffp_2d,dumc,un)

implicit none

double precision r0,z0,dr,dz,t1,t2,t3,t4,area_eff,zeta,un
double precision psibnd,dumc
integer k,nx
double precision ppp_2d(nx),ffp_2d(nx),psia_2d(nx)

if (un.gt.1.) dumc=0.0
if (un.le.0.) dumc=ppp_2d(1)*r0+ffp_2d(1)/r0
if (un.eq.1.) dumc=ppp_2d(nx)*r0+ffp_2d(nx)/r0

if (un.le.1..and.un.ge.0.) then
	zeta=un*(nx-1.)+1.
	k = floor(zeta)
	if (k.ge.nx) k=nx-1
	if (k.lt.1) k=1
	dumc=((ppp_2d(k+1)*(zeta-k)+ppp_2d(k)*(k+1.-zeta))*r0+ &
			& (ffp_2d(k+1)*(zeta-k)+ffp_2d(k)*(k+1.-zeta))/r0)
endif

return
end subroutine fill_in_current

