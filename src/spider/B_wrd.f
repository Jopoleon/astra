      subroutine wrb

      use sp_parameters
      use iopath, only: path
      use keys, only: kpr
      use compol, only: iplas

      implicit none

      integer :: i, j, ib, nrr
      real*8, dimension(nrp) :: psf, sqtor, bj_av, curfi_av, b2_av
      character(len=8) :: etitl(5)
      character(len=120) :: fname

      common /compsf/ psf, sqtor
      common /com_b2/ B2_av

      if(kpr.lt.0) return

      nrr=iplas

 2020 format(5e16.9)
 2022 format(2i5)

      return
      end subroutine wrb

!----------------------------------------------------------------
      subroutine tab_efit(tokf, psax, eqdfn, rax, zax, b0, r0)
         
      use ppf_modul
      use bnd_modul
      use sp_parameters, only: pi, amu0, twopi
      use iopath, only: path
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

      write(fname,'(a,a40)') TRIM(path),eqdfn
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

      implicit none

      return
      end subroutine wr_spik
