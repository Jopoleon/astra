      real*8 function xzer(x1, x2, v1, v2)

      implicit none

      real*8, intent(in) :: x1, x2, v1, v2
      xzer = (x1*v2 - x2*v1)/(v2 - v1 + 1.d-8)

      return
      end function xzer

!----------------------------------------------------------------
      subroutine loop

      use sp_parameters, only: nip, njp, nbndp2
      use comblc, only: ni, nj, ni1, imax, jmax, alpnew, un, r, z

      implicit none

      integer :: i, j, ig, ic, jc, ic1, jc1, lin, nxb
      real*8 :: delunb
      real*8, dimension(nbndp2) :: rxb, zxb
      real*8, dimension(nip, njp) :: us
      real*8, external :: xzer

      common/comlop/ rxb, zxb, nxb

      delunb=0.00001d0

      do i=1,ni
         do j=1,nj
            us(i,j) = alpnew * (un(i,j) - 1.d0) + 1.d0
            un(i,j) = un(i,j) - delunb
         enddo
      enddo

      ig=1
      jc=jmax
      j=jmax

      do i=imax,ni1
         if(un(i,j)*un(i+1,j).le.0.) then
            ic=i
            rxb(ig)=xzer(r(i+1),r(i),un(i+1,j),un(i,j))
            zxb(ig)=z(j)
            goto 886
         endif
      enddo
      stop

 886  continue

      ic1=ic
      jc1=jc
      lin=1

      do
         ig=ig+1
         i=ic
         j=jc

         if(un(i+1,j)*un(i,j).le.0. .AND. lin.ne.1) then
            ic  = i
            jc  = j-1
            lin = 3
            rxb(ig) = xzer(r(i+1),r(i),un(i+1,j),un(i,j))
            zxb(ig) = z(j)
         elseif(un(i+1,j+1)*un(i+1,j).le.0..AND. lin.ne.2) then
            ic  = i+1
            jc  = j
            lin = 4
            rxb(ig) = r(i+1)
            zxb(ig) = xzer(z(j+1),z(j),un(i+1,j+1),un(i+1,j))
         elseif(un(i+1,j+1)*un(i,j+1).le.0. .AND. lin.ne.3) then
            ic = i
            jc = j+1
            lin = 1
            rxb(ig) = xzer(r(i+1),r(i),un(i+1,j+1),un(i,j+1))
            zxb(ig) = z(j+1)
         elseif(un(i,j+1)*un(i,j).le.0..AND. lin.ne.4) then
            ic  = i-1
            jc  = j
            lin = 2
            rxb(ig) = r(i)
            zxb(ig) = xzer(z(j+1),z(j),un(i,j+1),un(i,j))
         endif

! control of finish

         if( (jc.eq.jc1) .and. (ic.eq.ic1) ) then
            EXIT
         endif

      enddo

      nxb=ig

!---------------------------------
! End of plasma boundary treatment
!---------------------------------

      do i=1,ni
         do j=1,nj
            un(i,j) = un(i,j) + delunb
         enddo
      enddo

      return
      end subroutine loop

!----------------------------------------------------------------
      subroutine btpol(betpol)

      use sp_parameters, only: pi, amu0
      use comblc, only: ni1, nj1, ipr, um, un, up, dri, dzj, cnor, tok

      implicit none

      real*8, intent(out) :: betpol

      integer :: i, j, iprr
      real*8 :: psi, zpres, pintg
      real*8, external :: funppp

      pintg=0.d0
      do i=2,ni1
         do j=2,nj1
            iprr=ipr(i, j)
            if(iprr.ne.1) CYCLE
            psi=un(i, j)
            zpres=funppp(psi)
            pintg=pintg+zpres*dri(i)*dzj(j)
         enddo
      enddo
      BETpol=8.d0*pi*pintg/(tok*tok)
      BETpol=betpol*(um-up)*cnor/amu0**2

      return
      end subroutine btpol

!----------------------------------------------------------------
      subroutine wrd

      use sp_parameters, only: nbndp2
      use iopath, only: kname, path
      use comblc, only: ni, ni1, ni2, nj, nj1, nj2, ipr, 
     &   r, z, rm, zm, rx0, zx0, r, z, curf,
     &   u, ue, ui, um, un, up, ux0

      implicit none

      integer :: i, j, l, ig, nxb
      real*8, dimension(nbndp2) :: rxb, zxb
      character(len=80) :: fname

      common/comlop/ rxb, zxb, nxb

      write(fname,'(a,a)') path(1:kname),'out.wr'
      open(1,file=fname,form='formatted')
         write(1,*) ni,nj,ni1,nj1,ni2,nj2,nxb
         write(1,*) (r(i),i=1,ni)
         write(1,*) (z(j),j=1,nj)
         write(1,*) ((u(i,j),i=1,ni),j=1,nj)
         write(1,*) ((curf(i,j),i=1,ni),j=1,nj)
         write(1,*) ((ipr(i,j),i=1,ni),j=1,nj)
         write(1,*) rm,zm,um,rx0,zx0,ux0,up
         write(1,*) (rxb(ig),ig=1,nxb)
         write(1,*) (zxb(ig),ig=1,nxb)
         write(1,*) ((ui(i,j),i=1,ni),j=1,nj)
         write(1,*) ((ue(i,j),i=1,ni),j=1,nj)
         write(1,*) ((un(i,j),i=1,ni),j=1,nj)
      close(1)

      write(fname,'(a,a)') path(1:kname),'recbon.wr'
      open(1,file=fname,form='formatted')  
         write(1,*) nxb
         do ig=1,nxb
            write(1,*) rxb(ig),zxb(ig)
         enddo     
      close(1)

      return
      end subroutine wrd

!----------------------------------------------------------------
      subroutine wrrec

      use iopath, only: kname, path
      use comblc, only: ni, ni1, ni2, nj, nj1, nj2, imax, jmax, ipr,
     &   r, z, rx0, zx0, rm, zm, r0ax, rx1, rx2, zx1, zx2, 
     &   rmin, rmax, zmin, zmax, qcen, b0ax, u, ue, um, un, up, ux0

      implicit none

      integer :: i, j, l
      character(len=80) :: fname

      write(fname,'(a,a)') path(1:kname),'rect.wr'
      open(1,file=fname,form='formatted')
         write(1,*) ni,nj,ni1,nj1,ni2,nj2,imax,jmax
         write(1,*) (r(i),i=1,ni)
         write(1,*) (z(j),j=1,nj)
         write(1,*) ((u(i,j),i=1,ni),j=1,nj)
         write(1,*) ((ue(i,j),i=1,ni),j=1,nj)
         write(1,*) ((un(i,j),i=1,ni),j=1,nj)
         write(1,*) ((ipr(i,j),i=1,ni),j=1,nj)
         write(1,*) rm,zm,um,rx0,zx0,ux0,up,qcen,b0ax,r0ax
         write(1,*) rx1,zx1,rx2,zx2
         write(1,*) rmax,zmax,rmin,zmin
      close(1)

      return
      end subroutine wrrec

!----------------------------------------------------------------
      subroutine wrdbnd
 
      use iopath, only: kname, path
      use comblc, only: nbnd, binadg

      implicit none

      integer :: i, j
      character(len=80) :: fname

      write(fname,'(a,a)') path(1:kname),'bnd.wr'
      open(1,file=fname,form='formatted')
         write(1,*) ((binadg(i,j),i=1,nbnd),j=1,nbnd)
      close(1)

      return
      end subroutine wrdbnd

!----------------------------------------------------------------
      subroutine rddbnd

      use iopath, only: kname, path
      use comblc, only: nbnd, binadg

      implicit none

      integer :: i, j
      character(len=80) :: fname

      write(fname,'(a,a)') path(1:kname),'bnd.wr'
      open(1,file=fname,form='formatted')
         read(1,*) ((binadg(i,j),i=1,nbnd),j=1,nbnd)
      close(1)

      return
      end subroutine rddbnd
