      subroutine L_MATR(NC, NCPFC, NTYPE, RC, ZC, VC, HC, NECON, WECON)
! DEFINITION MUTUALS INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  PPIND

      use sp_parameters, only: njlim, nplim, nclim, nilim, npfc0, 
     &   pi, amu0
      use iopath, only: path
      use comevl, only: ppind

      implicit none

      integer, intent(in) :: NC, NCPFC
      integer, intent(in), dimension(NILIM) :: NECON
      integer, intent(in), dimension(NCLIM) :: NTYPE
      real*8, intent(in), dimension(NCLIM) :: RC, ZC, VC, HC
      real*8, intent(in), dimension(NILIM) :: WECON

      integer :: i, j, j1, j2, ncequi, nequi
      real*8, external :: selind, betind
      character(len=80) :: fname

      do J1=1,NC
         do J2=1,NC
            if( J2.LT.J1 ) then
              PPIND(J1,J2) = PPIND(J2,J1)
            endif
            if( J2.EQ.J1 ) then
              PPIND(J1,J2) = SELIND( NTYPE(J1),RC(J1),VC(J1),HC(J1) )
            endif
            if( J2.GT.J1 ) then
              PPIND(J1,J2) = BETIND( NTYPE(J1),
     &                               RC(J1), ZC(J1), VC(J1), HC(J1),
     &                               NTYPE(J2),
     &                               RC(J2), ZC(J2), VC(J2), HC(J2) )
            endif
         enddo
      enddo

      do j1=1,NC
         do j2=1,NC
            ppind(j1,j2) = ppind(j1,j2)*amu0
         enddo
      enddo

      CALL TRAMAT( PPIND, NJLIM, NC, NCPFC, NEQUI,
     &             NECON, WECON )

      NCEQUI = NC - NCPFC + NEQUI

      write(fname,'(a,a)') TRIM(path), '/ppind_mat.wr'
      open(1,file=fname,form='formatted')
         write(1,*) ncequi
         write(1,*) ((ppind(i,j),i=1,ncequi), j=1,ncequi)
      close(1)

      return
      end subroutine L_MATR

!----------------------------------------------------------------
      subroutine rd_ppind

      use sp_parameters, only: njlim, nplim, npfc0
      use iopath, only: path
      use comevl, only: ppind

      implicit none

      integer :: i, j, ncequi
      character(len=80) :: fname

      common/comeqg/ ncequi

      write(fname,'(a,a)') TRIM(path), '/ppind_mat.wr'
        open(1,file=fname,form='formatted')
           read(1,*) ncequi
           read(1,*) ((ppind(i,j),i=1,ncequi), j=1,ncequi)
        close(1)

      return
      end subroutine rd_ppind

!----------------------------------------------------------------
      subroutine PROPNT(NPRO, RPRO, ZPRO, FIPRO)
!--- INPUT OF POSITIONS OF "PF_PROBE" POINTS:

      use iopath, only: path

      implicit none

      integer, intent(out) :: NPRO
      real*8, intent(out), dimension(*) :: RPRO, ZPRO, FIPRO

      integer :: i, l
      character(len=80) :: fname

      write(fname,'(a,a)') TRIM(path), '/pf_probe.dat'
      open(1,file=fname,form='formatted')
         read(1,*)  NPRO
         if( NPRO.NE.0 ) then
            do L=1,NPRO
               read(1,*) RPRO(L), ZPRO(L), FIPRO(L)
            enddo
         endif
      close(1)

      write(fname,'(a,a)') TRIM(path), '/propoi.wr'
      open(1,file=fname,form='formatted')
         write(1,*) npro
         write(1,*) ( rpro(i), i=1,npro)
         write(1,*) ( zpro(i), i=1,npro)
         write(1,*) (fipro(i), i=1,npro)
      close(1)

      return
      end subroutine PROPNT

!----------------------------------------------------------------
      subroutine rd_prob( NPRO, RPRO, ZPRO,  FIPRO )
! READING POSITIONS OF "PF_PROBE" POINTS:

      use iopath, only: path

      implicit none

      integer, intent(out) :: NPRO
      real*8, intent(out), dimension(*) :: RPRO, ZPRO, FIPRO

      integer :: i
      character(len=80) :: fname

      write(fname,'(a,a)') TRIM(path), '/propoi.wr'
      open(1,file=fname,form='formatted')
         read(1,*) npro
         read(1,*) ( rpro(i), i=1,npro)
         read(1,*) ( zpro(i), i=1,npro)
         read(1,*) (fipro(i), i=1,npro)
      close(1)

      return
      end subroutine rd_prob

!-------------------------------------------------------------------
      subroutine LOOPNT(NINFW, NGRA1, NLOO, RLOO, ZLOO)
!--- INPUT OF POSITIONS OF "FL_LOOP" POINTS:

      use iopath, only: path

      implicit none

      integer, intent(out) :: NLOO, NINFW, NGRA1
      real*8, intent(out), dimension(*) :: RLOO, ZLOO

      integer :: i, l
      character(len=80) :: fname
 
      NINFW=1
      ngra1=1
 
      write(fname,'(a,a)') TRIM(path), '/fl_loop.dat'
      open(NINFW,file=fname,form='formatted')
         read(NINFW,*) NLOO
         if( NLOO.NE.0 ) then
            do L=1, NLOO
               read(NINFW,*) RLOO(L), ZLOO(L)
            enddo
         endif
      close(NINFW)
      write(fname,'(a,a)') TRIM(path), '/loopoi.wr'
      open(ngra1,file=fname,form='formatted')
         write(ngra1,*) nloo
         write(ngra1,*) (rloo(i), i=1,nloo)
         write(ngra1,*) (zloo(i), i=1,nloo)
      close(ngra1)

      return
      end subroutine LOOPNT

!----------------------------------------------------------------
      subroutine rd_loop( NLOO, RLOO, ZLOO )
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:

      use iopath, only: path

      implicit none

      integer, intent(out) :: NLOO
      real*8, intent(out), dimension(*) :: RLOO, ZLOO

      integer :: i, ngra1
      character(len=80) :: fname

      ngra1 = 1
      write(fname,'(a,a)') TRIM(path), '/loopoi.wr'
      open(1,file=fname,form='formatted')
         read(ngra1,*) nloo
         read(ngra1,*) (rloo(i), i=1,nloo)
         read(ngra1,*) (zloo(i), i=1,nloo)
      close(ngra1)

      return
      end subroutine rd_loop

!----------------------------------------------------------------
      subroutine PASCUR(NEQUI, NFW, NBP, NVV, PJK, 
     &                  fwcurr, bpcurr, vvcurr )
! TOROIDAL CURRENTS OF PASSIVE CONDUCTOR STRUCTURES

      implicit none

      integer, intent(in) :: NEQUI, NFW, NBP, NVV
      real*8, intent(in), dimension(*) :: PJK
      real*8, intent(out) :: fwcurr, bpcurr, vvcurr

      integer :: i

      fwcurr = 0.d0
      if( NFW .NE. 0 ) then
         do I=1,NFW
            fwcurr = fwcurr + PJK(NEQUI+I)
         enddo
      endif
      bpcurr = 0.d0
      if( NBP .NE. 0 ) then
         do I=1,NBP
            bpcurr = bpcurr + PJK(NEQUI+NFW+I)
         enddo
      endif
      vvcurr = 0.d0
      if( NVV .NE. 0 ) then
         do I=1,NVV
            vvcurr = vvcurr + PJK(NEQUI+NFW+NBP+I)
         enddo
      endif

      return
      end subroutine PASCUR
