      subroutine wrb

      use sp_parameters
      use iopath, only: kname, path
      use keys, only: kpr

      implicit none

      integer :: i, j, ib, nrr
      real*8, dimension(nrp) :: psf, sqtor, bj_av, curfi_av, b2_av
      character(len=8) :: etitl(5)
      character(len=80) :: fname

      include 'compol.inc'

      common /compsf/ psf, sqtor
      common /com_jb/ BJ_av, curfi_av
      common /com_b2/ B2_av

      if(kpr.lt.0) return

      write(fname,'(a,a)') path(1:kname),'outp.wr'
      open(1,file=fname)
         write(1,*) nr,nt,nr1,nt1,nr2,nt2,iplas
         write(1,*) ((r(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((z(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((cur(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((psi(i,j),i=1,iplas),j=1,nt)
         write(1,*)  (q(i),i=1,iplas)
         write(1,*)  (f(i),i=1,iplas)
      close(1)

      write(fname,'(a,a)') path(1:kname),'ddp.wr'
      open(1,file=fname)
         write(1,*) iplas
         write(1,*) (q(i),i=1,iplas)
         write(1,*) (f(i),i=1,iplas)
         write(1,*) (dfdpsi(i),i=1,iplas)
         write(1,*) (psia(i),i=1,iplas)
         write(1,*) (sqtor(i),i=1,iplas)
         write(1,*) (dpdpsi(i),i=1,iplas)
         write(1,*) (BJ_av(i),i=1,iplas)
         write(1,*) (b2_av(i),i=1,iplas)
      close(1)
      write(fname,'(a,a)') path(1:kname),'tabppf.wr'
      open(1,file=fname)
         write(1,*) iplas
         do i=1,iplas
            write(1,*) 1.d0-psia(i),dpdpsi(i),dfdpsi(i)
         enddo
      close(1)

      write(fname,'(a,a)') path(1:kname),'q.wr'
      open(1,file=fname)
         do i=1,iplas
            if(i.ne.iplas) then
               write(1,*) 1.d0 - 0.5d0 * 
     &            (psia(i) + psia(i+1)), 0.5d0*q(i)/pi, i
            else
               write(1,*) 1.d0-psia(i),0.5d0*q(i)/pi,i
            endif
         enddo
      close(1)

      nrr=iplas

      write(fname,'(a,a)') path(1:kname),'efit_comp.wr'
      open(1,file=fname)
         write(1,2022) nrr,nt
         write(1,2020) rm,zm,psim*0.4d0*pi,psip*0.4d0*pi,tok*1.d3
         write(1,2020) (f(i)*0.4d0*pi,i=1,nrr-1)
         write(1,2020) (dpdpsi(i)*1.d7/4.d0/pi,i=1,nrr)
         write(1,2020) (dfdpsi(i)*0.4d0*pi,i=1,nrr)
         write(1,2020) ((r(i,j),i=1,nrr),j=1,nt)
         write(1,2020) ((z(i,j),i=1,nrr),j=1,nt)
         write(1,2020) ((psi(i,j)*0.4d0*pi,i=1,nrr),j=1,nt)
         write(1,2020) (q(i),i=1,nrr-1)
         write(1,2020) (r(nrr,j),z(nrr,j),j=1,nt)
      close(1)

      write(fname,'(a,a)') path(1:kname),'tab_bnd.wr'
      open(1,file=fname)
	 write(1,*) nt1 
	 do ib=1,nt1
	    write(1,*) r(iplas,ib),z(iplas,ib) 
	 enddo
      close(1) 

 2020 format(5e16.9)
 2022 format(2i5)

      return
      end subroutine wrb

!----------------------------------------------------------------
      subroutine tab_efit(tokf, psax, eqdfn, rax, zax, b0, r0)
         
      use ppf_modul
      use bnd_modul
      use sp_parameters, only: pi, amu0, twopi
      use iopath, only: kname, path
      use keys, only: kpr

      implicit none

      integer, parameter :: np=1000, nbp=np*4

      character(len=40), intent(in) :: eqdfn
      real*8, intent(out) :: tokf, psax, rax, zax, b0, r0
 
      integer :: i, j, iefit, idum, nw, nh, nbbbs, limitr, i_sign
      real*8 :: simag, sibry, rdim, zdim, rcentr, rleft, zmid, rmaxis,
     &   zmaxis, bcentr, current, xdum, p_intgr, f0c, fvefit, fcefit
      real*8, dimension(np) :: ps, p, f, q, fpol, pres, qpsi, ffprim,
     &   pprime, rlimtr, zlimtr, x, y
      real*8, dimension(np,np) :: psirz, u
      real*8, dimension(nbp) :: rbbbs, zbbbs

      character(len=8) :: ccase(6)
      character(len=80) :: fname

      common/com_eqd/ fpol, pres, qpsi, 
     &                ffprim, pprime, 
     &                rlimtr, zlimtr, ccase, simag, sibry,
     &                idum, nw, nh, limitr

      common/efites/ fcefit,rcentr,iefit
      common/comefi/ x, y, u

      iefit=1                  

      write(*,*) '---------------------------'
      write(*,*) ' Entry of subr."tab_efit":'
      write(*,*) '--------------------------'

      write(fname,'(a,a40)') path(1:kname),eqdfn
      open(1,file=fname,form='formatted')
	 read(1,2000) (ccase(i),i=1,6),idum,nw,nh
         write(*,*) idum,nw,nh
         read(1,2020) rdim,zdim,rcentr,rleft,zmid
         read(1,2020) rmaxis,zmaxis,simag,sibry,bcentr
         read(1,2020) current,simag,xdum,rmaxis,xdum
         read(1,2020) zmaxis,xdum,sibry,xdum,xdum
         read(1,2020) (fpol(i),i=1,nw)
         read(1,2020) (pres(i),i=1,nw)
         read(1,2020) (ffprim(i),i=1,nw)
         read(1,2020) (pprime(i),i=1,nw)
         read(1,2020) ((psirz(i,j),i=1,nw),j=1,nh)
         read(1,2020) (qpsi(i),i=1,nw)
         read(1,2022) nbbbs,limitr
         read(1,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
         read(1,2020) (rlimtr(i),zlimtr(i),i=1,limitr)
      close(1)
      psirz=psirz/twopi
      ffprim=ffprim*twopi
      pprime=pprime*twopi

      p_intgr=0.0d0

      do i=1,nw
         p_intgr=p_intgr+pprime(i)                           
      enddo             

      if(p_intgr.lt.0.d0) then
         i_sign=-1
      else 
         i_sign= 1
      endif               

      simag=simag*i_sign
      sibry=sibry*i_sign
      bcentr=dabs(bcentr)
      current=dabs(current)
      do i=1,nw
         fpol(i)=fpol(i)*i_sign
         ffprim(i)=ffprim(i)*i_sign
         pprime(i)=pprime(i)*i_sign
         do j=1,nh
            psirz(i,j)=psirz(i,j)*i_sign
         enddo             
      enddo             
          
      rax = rmaxis
      zax = zmaxis
      b0  = bcentr
      r0  = rcentr
      f0c = b0*rcentr

      fvefit=fpol(nw)

 2000 format(6a8,3i4)
 2020 format(5e16.9)
 2022 format(2i5)

      if(allocated(pstab)) deallocate( pstab, pptab, fptab )
      allocate( pstab(nw), pptab(nw), fptab(nw) )
     
      nutab=nw

      do i=1,nw
         pstab(i)= dfloat(i-1)/dfloat(nw-1)
         pptab(i)= pprime(i)*amu0*1.d-6
         fptab(i)= ffprim(i)         
      enddo

      nbtab = nbbbs
      if(allocated(rbtab)) deallocate( rbtab, zbtab )
      allocate( rbtab(nbbbs), zbtab(nbbbs) )

      j=0
      do i=nbbbs,1,-1
         j=j+1
         rbtab(j)=rbbbs(i)                 
         zbtab(j)=zbbbs(i)                 
      enddo

      fcefit = fpol(nw)
      tokf   = current*1.d-6
      psax   = -(sibry-simag)

      return
      end subroutine tab_efit

!----------------------------------------------------------------
      subroutine wr_spik

      use sp_parameters
      use iopath, only: kname, path
      use keys, only: kpr

      implicit none

      integer :: i, j, nm
      character(len=80) :: fname

      include 'compol.inc'

      if(kpr.lt.0) return

      write(fname,'(a,a)') path(1:kname),'spik.wr'
      open(1,file=fname,form='formatted')
         nm=iplas*nt1
         write(1,*) iplas,nt1,nm,psim,psibon,1
         write(1,*)(dsqrt(1.d0-psia(i)),i=1,iplas),
     &      (dpdpsi(i),i=1,iplas),
     &      (dfdpsi(i),i=1,iplas),
     &      (r(1,j),j=1,nt1),
     &      (z(1,j),j=1,nt1),
     &      (r(iplas,j),j=1,nt1),
     &      (z(iplas,j),j=1,nt1),
     &      ((ro(i,j)/ro(iplas,j),j=1,nt1),i=1,iplas),
     &      (q(i)/2.d0/pi,i=1,iplas),fvac
      close(1)

      return
      end subroutine wr_spik
