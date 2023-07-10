!C======================================================================|
      subroutine XEQDSCA
!C----------------------------------------------------------------------|
!c EQDSK creator from annular equilibrium, valid also for SPIDER adaptive grid!
!C
!C E Fable 2012
!C
!C----------------------------------------------------------------------|
	 use parameter_inc
	 use const_inc
	 use status_inc
	 use outcmn_inc
	 use parameters_a2equil, only: equil_now
		
      implicit none

      integer		j,jf,length
      double precision r,z,br,bz,bt

      integer jna,i,jj
 
	integer Nrrect,Nzrect,icnt,ibgout
	integer Nx,Nt
	integer neqdsk,imtehor
	double precision Rrect(300),Zrect(300),aaa
	double precision PSI_rect(300,300)
	double precision ccc(300,300)
	double precision 	pres_rect(300)
	double precision 	fdia_rect(300)
	double precision 	q_rect(300)
	double precision 	fprime_rect(300)
	double precision 	pprime_rect(300),rhopolr(300)
	double precision 	B_R(300,300)
	double precision 	r_a(300,300)
	double precision 	thetap(300)
	double precision 	B_T(300,300)
	double precision 	B_Z(300,300)
	double precision 	XX(300,300)
	double precision 	YY(300,300)
	double precision 	PSI(300,300)
	double precision 	Rb(300)
	double precision 	Zb(300),pressure(300),qprof(300)
	double precision 	PSIn_grid(300)
	double precision 	dum1(300),dum2(300)
	double precision 	ipol_grid(300)
	double precision 	PSIb,amp_fac
	character*83 file2d_equilr
	character*10 casez(6)
							
	file2d_equilr="exp/equ/"//trim(MACHINE)//"/eqdska.gsef"

nx=nequil
nt=mequil

thetap(1:nt)=equil_now%coord_sys%position%teta2d(1:nt)
rb(1:nt)=equil_now%coord_sys%position%r(nx,1:nt)
zb(1:nt)=equil_now%coord_sys%position%z(nx,1:nt)
xx(1:nx,1:nt)=equil_now%coord_sys%position%r(1:nx,1:nt)
yy(1:nx,1:nt)=equil_now%coord_sys%position%z(1:nx,1:nt)
psi(1:nx,1:nt)=equil_now%coord_sys%position%psirz(1:nx,1:nt)
r_a(1:nx,1:nt)=equil_now%coord_sys%position%rmin(1:nx,1:nt)
psin_grid(1:nx)=equil_now%profiles_1d%psi(1:nx)
pressure(1:nx)=equil_now%profiles_1d%pressure(1:nx)
ipol_grid(1:nx)=equil_now%profiles_1d%F_dia(1:nx)
dum1(1:nx)=equil_now%profiles_1d%pprime(1:nx)
dum2(1:nx)=equil_now%profiles_1d%ffprime(1:nx)
qprof(1:nx)=equil_now%profiles_1d%q(1:nx)

PSIb=PSI(Nx,1)

!divide in prescribed and free boundary
!	write(*,*) 'ddd1'

if (TIME.le.ITFBE) then
 amp_fac=5.
 aaa=abc
 imtehor=4
 nrrect=256
 nzrect=257
!	write(*,*) 'ddd2',nx,nt,nrrect,nzrect,imtehor,psib,psi(1:nx,1),psin_grid(1:nx)
!!	write(*,*) 'ddd2',thetap(1:nt),r_a(1:nx,1)
! From polar to rectangluar grid
	call	fsg2recg_gsef(Nx,Nt,PSI(1:Nx,1:Nt),PSIb, &
     &  aaa,ipol_grid(1:nx),imtehor, &
     &  XX(1:Nx,1:Nt),YY(1:Nx,1:Nt),&
     &  r_a(1:Nx,1:Nt),thetap(1:Nt),&
     &  Rrect(1:nrrect),Zrect(1:nzrect),nrrect,nzrect,&
     &  psi_rect(1:nrrect,1:nzrect),&
     &  ccc(1:nrrect,1:nzrect),&
     &  ccc(1:nrrect,1:nzrect),&
     &  ccc(1:nrrect,1:nzrect),amp_fac)
!	write(*,*) 'ddd3'

else

 nrrect=equil_now%eqgeometry%rectgrid%npointsr
 nzrect=equil_now%eqgeometry%rectgrid%npointsz
 rrect(1:nrrect)=equil_now%eqgeometry%rectgrid%r2d(1:nrrect)
 zrect(1:nzrect)=equil_now%eqgeometry%rectgrid%z2d(1:nzrect)
 psi_rect(1:nrrect,1:nzrect)=equil_now%eqgeometry%rectgrid%psirz2d(1:nrrect,1:nzrect)

endif



	do i=1,Nrrect
	rhopolr(i)=((i-1.)/(Nrrect-1.))**0.5
	enddo


	call qinterp(PSIn_grid(1:Nx),&
     &  pressure(1:Nx),Nx,&
     &  rhopolr(1:Nrrect),pres_rect(1:Nrrect),Nrrect)
	call qinterp(PSIn_grid(1:Nx),&
     &  ipol_grid(1:Nx),Nx,&
     &  rhopolr(1:Nrrect),fdia_rect(1:Nrrect),Nrrect)
	call qinterp(PSIn_grid(1:Nx),&
     &  dum2(1:Nx),Nx,&
     &  rhopolr(1:Nrrect),fprime_rect(1:Nrrect),Nrrect)
	call qinterp(PSIn_grid(1:Nx),&
     &  dum1(1:Nx),Nx,&
     &  rhopolr(1:Nrrect),pprime_rect(1:Nrrect),Nrrect)
	call qinterp(PSIn_grid(1:Nx),&
     &  qprof(1:Nx),Nx,&
     &  rhopolr(1:Nrrect),q_rect(1:Nrrect),Nrrect)


!write output file

	neqdsk=32
	do i=1,6
	casez(i)='XXXXX     '
	enddo
	open(neqdsk,file=file2d_equilr)
	write (neqdsk,2000) (casez(i),i=1,6),3,Nrrect,Nzrect
	write (neqdsk,2020) Rrect(Nrrect)-Rrect(1),&
     & Zrect(Nzrect)-Zrect(1),rtor,&
     & Rrect(1),(Zrect(1)+Zrect(Nzrect))/2.
	write (neqdsk,2020) XX(1,1),YY(1,1),PSI(1,1),PSIb,btor
	write (neqdsk,2020) ipl*1.d6,PSI(1,1),1.,XX(1,1),1.
	write (neqdsk,2020) YY(1,1),1.,Psib,1.,1.
	write (neqdsk,2020) (fdia_rect(i),i=1,Nrrect)
	write (neqdsk,2020) (pres_rect(i),i=1,Nrrect)
	write (neqdsk,2020) (fprime_rect(i),i=1,Nrrect)
	write (neqdsk,2020) (pprime_rect(i),i=1,Nrrect)
	write (neqdsk,2020) ((psi_rect(i,j),i=1,Nrrect),j=1,nzrect)
	write (neqdsk,2020) (q_rect(i),i=1,Nrrect)
	write (neqdsk,2022) nt,nt
	write (neqdsk,2020) (XX(nx,i),YY(nx,i),i=1,nt)
	write (neqdsk,2020) (XX(nx,i),YY(nx,i),i=1,nt)
	close(neqdsk)
2000	format (6a8,3i4)
2020	format (5e16.9)
2022	format (2i5)

return
end

