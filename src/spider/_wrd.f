      SUBROUTINE wrd_tim

      use sp_parameters, only: nstep_p
      use iopath, only: path
      use comtim
      use keys, only: kstep, kpr

      implicit none

      integer :: i, kwr
      character(len=80) :: fname

      if (kpr.ge.0) then
         write(fname,'(a,a)') TRIM(path), '/eq_tim.wr'
         open(1,file=fname)
            kwr=kstep+1
            write(1,*)  kwr
            write(1,*)  (time_t(i),i=1,kwr)
            write(1,*)  (torcur(i),i=1,kwr)
            write(1,*)  (rm_t(i),i=1,kwr)
            write(1,*)  (zm_t(i),i=1,kwr)
            write(1,*)  (rxp_t(i),i=1,kwr)
            write(1,*)  (zxp_t(i),i=1,kwr)
            write(1,*)  (betpol_t(i),i=1,kwr)
            write(1,*)  (bettor_t(i),i=1,kwr)
            write(1,*)  (psim_t(i),i=1,kwr)
            write(1,*)  (psib_t(i),i=1,kwr)      
         close(1)
      endif

      return
      end subroutine wrd_tim

!----------------------------------------------------------------
      subroutine f_wrd

      use sp_parameters, only: nrp
      use iopath, only: path
      use keys, only: kpr, kstep
      use jb, only: Bj_av
      use tim, only: ctim
      use compol_add, only: g, psii, psie, aex, rx0, zx0
      use compol, only: nr, nt, nr1, nt1, nr2, nt2, iplas, r, z,
     &    cur, psi, q, f, rm, zm, psia, dpdpsi, dfdpsi

      implicit none

      integer :: i, j
      real*8 :: ddps
      real*8, dimension(nrp) :: psf, sqtor, b2_av, sigma, cbut_b
      character(len=80) :: fname

      common/compsf/ psf, sqtor
      common /com_b2/ B2_av
      common/com_but/ sigma,cbut_b

      if (kpr.ge.0) then
         write(fname,'(a,a)') TRIM(path), '/outp.wr'
         open(1,file=fname)
            write(1,*) nr,nt,nr1,nt1,nr2,nt2,iplas
            write(1,*) ((r(i,j),i=1,nr),j=1,nt)
            write(1,*) ((z(i,j),i=1,nr),j=1,nt)
            write(1,*) ((g(i,j),i=1,nr),j=1,nt)
            write(1,*) ((cur(i,j),i=1,nr),j=1,nt)
            write(1,*) ((psi(i,j),i=1,nr),j=1,nt)
            write(1,*) ((psii(i,j),i=1,nr),j=1,nt)
            write(1,*) ((psie(i,j),i=1,nr),j=1,nt)
            write(1,*) ((aex(i,j),i=1,nr),j=1,nt)
            write(1,*)  (q(i),i=1,iplas)
            write(1,*)  (f(i),i=1,iplas)
            write(1,*)  rm,zm,rx0,zx0,ctim,kstep
        close(1)

        write(fname,'(a,a)') TRIM(path), '/ddp.wr'
        open(1,file=fname)
           write(1,*) iplas
           write(1,*) (q(i),i=1,iplas)
           write(1,*) (f(i),i=1,iplas)
           write(1,*) (dfdpsi(i),i=1,iplas)
           write(1,*) (psia(i),i=1,iplas)
           write(1,*) (psf(i),i=1,iplas)
           write(1,*) (dpdpsi(i),i=1,iplas)
           write(1,*) (BJ_av(i),i=1,iplas)
           write(1,*) (b2_av(i),i=1,iplas)
           write(1,*) (sigma(i),i=1,iplas)
           write(1,*) (cbut_b(i),i=1,iplas)
           write(1,*)  ctim,kstep
        close(1)
        write(fname,'(a,a)') TRIM(path), '/dps.wr'
        open(1,file=fname)
           do i=1,iplas
              ddps=psi(i,2)-psf(i)
              write(1,*) ddps,i
           enddo
        close(1)
      endif

      return
      end subroutine f_wrd

!----------------------------------------------------------------
      subroutine rdrec

      use sp_parameters, only: nip, njp
      use iopath, only: path
      use comrec

      implicit none

      integer :: i, j
      real*8 :: qcen, b0ax, r0ax
      character(len=80) :: fname

      write(fname,'(a,a)') TRIM(path), '/rect.wr'
      open(1,file=fname)
         read(1,*) ni,nj,ni1,nj1,ni2,nj2,imax,jmax
         read(1,*) (x(i),i=1,ni)
         read(1,*) (y(j),j=1,nj)
         read(1,*) ((u(i,j),i=1,ni),j=1,nj)
         read(1,*) ((ue(i,j),i=1,ni),j=1,nj)
         read(1,*) ((un(i,j),i=1,ni),j=1,nj)
         read(1,*) ((ipr(i,j),i=1,ni),j=1,nj)
         read(1,*) xm,ym,um,xx0,yx0,ux0,up,qcen,b0ax,r0ax
         read(1,*) xx1,yx1,xx2,yx2
         read(1,*) xmax,ymax,xmin,ymin
      close(1)

      return
      end subroutine rdrec

!----------------------------------------------------------------
      subroutine wr_step(numwr, time, istep)

      use iopath, only: path
      use tim, only: dtim, ctim
      use compol, only: nr, nt, iplas, q, ro, r, z, rm, zm,
     &   psi, psia, psin, psim, psip, psi_eav, psibon0,
     &   teta, f, dfdpsi, dpdpsi, tok

      implicit none

      integer, intent(in) :: numwr, istep
      real*8, intent(in) :: time
      integer :: i, j
      character(len=40) :: str, dummy
      character(len=80) :: fname

      write(fname,'(a,a)') TRIM(path), '/nmwr.wr'
      open(1,file=fname)
         write(1,*) numwr
      close(1)

      if(numwr.lt.10) then
         write(str,'(a,a,i1,a)') TRIM(path), '/step',numwr,'.wr'
      elseif(numwr.lt.100) then
         write(str,'(a,a,i2,a)') TRIM(path), '/step',numwr,'.wr'
      elseif(numwr.lt.1000) then
         write(str,'(a,a,i3,a)') TRIM(path), '/step',numwr,'.wr'
      else
         write(str,'(a,a,i4,a)') TRIM(path), '/step',numwr,'.wr'
      endif

      open(1,file=str,form='formatted')
         write(1,*) nr,nt,iplas,istep,dtim,ctim
         write(1,*) ((r(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((z(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((ro(i,j),i=1,iplas),j=1,nt)
         write(1,*) (teta(j),j=1,nt)
         write(1,*) ((psi(i,j),i=1,iplas),j=1,nt)
         write(1,*) ((psin(i,j),i=1,iplas),j=1,nt)
         write(1,*) (psia(i),i=1,iplas)
         write(1,*)  (q(i),i=1,iplas)
         write(1,*)  (f(i),i=1,iplas)
         write(1,*) (dfdpsi(i),i=1,iplas)
         write(1,*) (dpdpsi(i),i=1,iplas)
         write(1,*) psi_eav,rm,zm,psim-psip,psibon0,tok
      close(1)

      write(fname,'(a,a)') TRIM(path), '/wlist.wr'
      open(1,file=fname)
         if(numwr.eq.1) then
            write(1,*) str
         else
            do i=1,numwr-1
               read(1,*) dummy
            enddo
            write(1,*) str
         endif
      close(1)

      return
      end subroutine wr_step
