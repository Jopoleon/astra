      subroutine DIFFER( AOLD, ANEW, N, ERRLOC, ERRVEC,
     &                   ABSMAX, ABSMIN, NOUT, NTER )

      implicit none

      integer, intent(in) :: N, NOUT, NTER
      real*8, intent(in), dimension(*) :: AOLD, ANEW
      real*8, intent(out) :: ERRLOC, ERRVEC, ABSMAX, ABSMIN

      integer :: I, L
      real*8 :: ABSA, ABSD, ABSR, S1, S2

      ERRLOC = 0.D0
      ABSMAX = dabs( ANEW(1) )
      ABSMIN = dabs( ANEW(1) )

      do I=1,N
         ABSA = dabs( ANEW(I) )
         if( ABSA.LT.ABSMIN ) ABSMIN = ABSA
         if( ABSA.GT.ABSMAX ) ABSMAX = ABSA
         ABSD = dabs( ANEW(I) - AOLD(I) )
         if( ABSA.GT.1.0D-12 ) then
            ABSR = ABSD / ABSA
            if( ABSR.GT.ERRLOC )  ERRLOC = ABSR
	 endif
      enddo

      S1 = 0.D0
      S2 = 0.D0
      do L=1,N
         S1 = S1 + ( ANEW(L) - AOLD(L) )**2
         S2 = S2 +   ANEW(L)**2
      enddo
      S1 = dsqrt(S1)
      S2 = dsqrt(S2)

      if( S2 .GT. 1.0D-12 ) then
         ERRVEC = S1 / S2
      else
         ERRVEC = 0.D0
      endif

      return
      end subroutine DIFFER

!------------------------------------------------------------------
!  NTYPE=1 - CONDUCTOR WITH "LINEAR"      CROSS-SECTION
!  NTYPE=2 - CONDUCTOR WITH "RECTANGULAR" CROSS-SECTION
!  NTYPE=3 - CONDUCTOR WITH "ROUND"       CROSS-SECTION
!  RC      - CILINDR. "R" COORDINATE  OF  CROSS-SECTION CENTRE
!  VC      - VERTICAL (OR RADIUS for NTYPE=3) SIZE OF CROSS-SECTION
!  HC      - HORIZONTAL                       SIZE OF CROSS-SECTION

      real*8 function SELIND( NTYPE, RC, VC, HC )

      use sp_parameters, only: pi, twopi

      implicit none

      integer, intent(in) :: NTYPE
      real*8, intent(in) :: RC, VC, HC

      integer :: i, j, NOUT, NTER
      real*8 :: VCL
      real*8, external :: greeni

      NOUT = 17
      NTER = 6

      if( NTYPE .EQ. 1 ) then

      VCL=dsqrt(HC**2+VC**2)

!--- FROM A.KAVIN
!
!+++      SELIND = RC * (dlog( 8D0*RC/VCL ) - 0.5D0) / twopi/1.0d0

      SELIND = 0.d0
      do i=-1,1
         do j=-1,i
            if(i .eq. j) then
               SELIND = SELIND + (RC+HC*i/3.d0) *
     &          (dlog( 8.D0*(RC+HC*i/3.d0)/(VCL/3.d0) ) - 0.5D0) / twopi
            else
               SELIND = SELIND +
     &            GREENI( RC+HC*i/3.d0, VC*i/3.d0,
     &                    RC+HC*j/3.d0, VC*j/3.d0)*2d0 /pi
            endif
         enddo
      enddo

      SELIND = SELIND/9.d0

      endif

      if( NTYPE .EQ. 2 ) then
         SELIND = RC*(dlog(8.D0*RC/(0.2236d0*(VC+HC)))-2.D0) / twopi
      endif
      if( NTYPE .EQ. 3 ) then
         SELIND = RC*(dlog(8.D0*RC/VC)-1.75D0) / twopi
      endif

      if((NTYPE.NE.1).AND.(NTYPE.NE.2).AND.(NTYPE.NE.3)) then
         STOP
      endif

      return
      end function SELIND

!***************************************************************
!  NTYPE=1 - CONDUCTOR WITH "LINEAR"      CROSS-SECTION
!  NTYPE=2 - CONDUCTOR WITH "RECTANGULAR"  CROSS-SECTION
!  NTYPE=3 - CONDUCTOR WITH "ROUND"        CROSS-SECTION
!  RC,ZC   - CILINDR. "R,Z" COORDINATES OF CROSS-SECTION CENTRE
!  VC      - VERTICAL (OR RADIUS for NTYPE=3) SIZE OF CROSS-SECTION
!  HC      - HORIZONTAL                       SIZE OF CROSS-SECTION
!
!  INDEX "1" MEANS THE FIRST  CONDUCTOR
!  INDEX "2" MEANS THE SECOND CONDUCTOR

      real*8 function BETIND( NTYPE1, RC1, ZC1, VC1, HC1,
     &                        NTYPE2, RC2, ZC2, VC2, HC2 )

      use sp_parameters, only: pi

      implicit none

      integer, intent(in) :: NTYPE1, NTYPE2
      real*8, intent(in) :: RC1, RC2, ZC1, ZC2, VC1, VC2, HC1, HC2

      integer :: i, j
      real*8, external :: greeni

      if( NTYPE1 .EQ. 1 .and. NTYPE2 .EQ. 1) then
         BETIND=0.d0
         do i=-1,1
            do j=-1,1
               BETIND = BETIND +
     &            GREENI( RC1+HC1*i/3.d0, ZC1+VC1*i/3.d0, 
     &                    RC2+HC2*j/3.d0, ZC2+VC2*j/3.d0 )
            enddo
         enddo
         BETIND = BETIND/PI/9.d0
      endif

      if( NTYPE1 .EQ. 1 .and. NTYPE2 .EQ. 2) then
         BETIND=0.d0
         do i=-1,1
            BETIND = BETIND +
     &         GREENI( RC1+HC1*i/3.d0, ZC1+VC1*i/3.d0, RC2, ZC2)
         enddo
         BETIND = BETIND/PI/3.d0
      endif

      if( NTYPE1 .EQ. 2 .and. NTYPE2 .EQ. 1) then
         BETIND=0.d0
         do j=-1,1
            BETIND = BETIND +
     &         GREENI( RC1, ZC1, RC2+HC2*j/3.d0, ZC2+VC2*j/3.d0)
         enddo
         BETIND = BETIND/PI/3.d0
      endif

      if( NTYPE1 .EQ. 2 .and. NTYPE2 .EQ. 2) then
         BETIND = GREENI( RC1, ZC1, RC2, ZC2 ) / PI
      endif

      if( NTYPE1 .EQ. 3 .OR.  NTYPE2 .EQ. 3) then
         BETIND = GREENI( RC1, ZC1, RC2, ZC2 ) / PI
      endif

      return
      end function BETIND

!----------------------------------------------------------
!  TRANSFORM OF "FULL" MATRIX "P" ( "NC*NC" SIZE ) TO
!  "EQUIVALENT" MATRIX "P" ( "NCEQUI*NCEQUI" SIZE ).
!  HERE   NCEQUI = NC - NCPFC + NEQUI.
!----------------------------------------------------------
      subroutine TRAMAT( P, NJLIM, NC, NCPFC, NEQUI, NECON, WECON )

      implicit none

      integer, intent(in) :: NJLIM, NC, NCPFC, NEQUI
      integer, intent(in), dimension(*) :: NECON
      real*8, intent(in), dimension(*) :: WECON
      real*8, intent(inout), dimension(njlim, njlim) :: P

      integer :: I, J, L, ncequi
      real*8, dimension(njlim) :: T

      NCEQUI = NC - NCPFC + NEQUI

!  SUMMING OF COLUMNS !!!
!  ----------------------
      do L=1,NEQUI
         do I=1,NC
            T(I) = 0.D0
         enddo
         do J=1,NCPFC
            if( NECON(J).EQ.L ) then
               do I=1,NC
                  T(I) = T(I) + P(I,J)*WECON(J)
               enddo
            endif
         enddo
         do I=1,NC
            P(I,L) = T(I)
         enddo
      enddo

      do J=NCPFC+1,NC
         do I=1,NC
            P(I,NEQUI + J - NCPFC) = P(I,J)
         enddo
      enddo

!  SUMMING OF ROWS !!!
!  -------------------
      do L=1,NEQUI
         do J=1,NCEQUI
            T(J) = 0.D0
         enddo
         do I=1,NCPFC
            if( NECON(I).EQ.L ) then
               do J=1,NCEQUI
                  T(J) = T(J) + P(I,J)*WECON(I)
               enddo
            endif
         enddo
         do J=1,NCEQUI
            P(L,J) = T(J)
         enddo
      enddo

      do I=NCPFC+1,NC
         do J=1,NCEQUI
            P(NEQUI + I - NCPFC,J) = P(I,J)
         enddo
      enddo

      return
      end subroutine TRAMAT

!------------------------------------------------------------------
      subroutine eq_par( z0cen_e, alp_e, alpnew_e, qcen_e,
     &                   nctrl_e, numlim_e, up_e,
     &                   rm_e, zm_e, rx0_e, zx0_e )

      use comblc, only: nctrl, numlim, alp, alpnew, qcen, up, 
     &   rx0, zx0, rm, zm

      implicit none

      integer, intent(out) :: nctrl_e, numlim_e
      real*8, intent(out) :: z0cen_e, alp_e, alpnew_e, qcen_e, 
     &   up_e, rm_e, zm_e, rx0_e, zx0_e

      real*8 :: bettot, rmx, zzrmx, rmn, zzrmn, zmx, rrzmx, zmn, rrzmn,
     &          r0cen, z0cen, radm, aspect,
     &          Eupper, Elower, DELup, DELlw, Bfvakc

      common /equili/ bettot, rmx, zzrmx, rmn, zzrmn,
     &                zmx, rrzmx, zmn, rrzmn,
     &                r0cen, z0cen, radm, aspect,
     &                Eupper, Elower, DELup, DELlw, Bfvakc

      z0cen_e  = z0cen
      alp_e    = alp
      alpnew_e = alpnew
      qcen_e   = qcen

      nctrl_e  = nctrl
      numlim_e = numlim
      up_e     = up
        
      rm_e  = rm
      zm_e  = zm
      rx0_e = rx0
      zx0_e = zx0

      return
      end subroutine eq_par

!***********************************************************************
!--- INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS

       subroutine CONDUC( NC, NCEQUI, NCPFC, NFW, NBP, NVV,
     &                    RC, ZC, PC, VC, HC, NTYPE,
     &                    RC1, ZC1, RC2, ZC2, RC3, ZC3, RC4, ZC4,
     &                    RES, VOLK, VOLKP1,
     &                    NECON, WECON,
     &                    NOUT, NTER, NINFW, ngra1 )

      use sp_parameters, only: njlim, nplim, npfc0, twopi, nclim
      use iopath, only: path
      use comevl, only: nequi, pfres, pfvol1

      implicit none

      integer, intent(in) :: NCPFC, NFW, NBP, NVV, ngra1
      integer, intent(in), dimension(*) :: NECON
      real*8, intent(in), dimension(*) :: WECON

      integer, intent(out) :: NC, NCEQUI, NOUT, NTER, NINFW
      integer, intent(out), dimension(*) :: NTYPE
      real*8, dimension(nclim) :: VC, HC
      real*8, intent(out), dimension(*) :: PC, RC, RC1, RC2, RC3, RC4,
     &   ZC, ZC1, ZC2, ZC3, ZC4, VOLK, VOLKP1
      real*8, intent(out), dimension(njlim, njlim) :: RES

      integer :: i, j, l, nsegbp, nsegvv, nsegfw
      real*8 :: BBB
      integer, dimension(nplim) :: KDFW, KDBP, KDVV, 
     &   NTYBP, NTYFV, NTYVV, NTYFW
      real*8, dimension(nplim) :: RP1FW, ZP1FW, RP2FW, ZP2FW, RFWSEG, 
     &   CFWSEG, RP1BP, ZP1BP, RP2BP, ZP2BP, RBPSEG, CBPSEG, RP1VV, 
     &   ZP1VV, RP2VV, ZP2VV, RVVSEG, CVVSEG, RFW, ZFW, DFW, HFW, 
     &   RESFW, CURFW, RFW1, ZFW1, RFW2, ZFW2, RBP, ZBP, DBP, HBP, 
     &   RESBP, CURBP, RBP1, ZBP1, RBP2, ZBP2, RVV, ZVV, DVV, HVV, 
     &   RESVV, CURVV, RVV1, ZVV1, RVV2, ZVV2
      character(len=80) :: fname

      BBB = 1.D0 / twopi

      do j=1,njlim
         do i=1,njlim
            res(i,j)=0.d0
         enddo
      enddo

!--- INPUT OF PFC SYSTEM PARAMETERS

      CALL TRECUR( NCPFC, RC, ZC, PC, NTYPE, NECON, WECON,
     &             HC, VC, NOUT, NTER )

      do l=1,NEQUI
         do j=1,NEQUI
            RES(L,j) = PFRES(L, j) * BBB
         enddo
      enddo

      do L=1,NEQUI
         VOLK(L) = PFVOL1(L) * BBB
      enddo

      NC     = NCPFC
      NCEQUI = NEQUI

!--- INPUT OF PASSIVE CONDUCTORS PARAMETERS

      call TREFW( NOUT, NTER, NINFW, NSEGFW,
     &            KDFW, RP1FW, ZP1FW, RP2FW, ZP2FW, RFWSEG, CFWSEG,
     &            NFW,  NTYFW,
     &            RFW,  ZFW,  DFW,   HFW,  RESFW, CURFW,
     &            RFW1, ZFW1, RFW2,  ZFW2 )

      if( NFW.NE.0 )  then
         do I=1,NFW
            NTYPE(NC+I) = NTYFW(I)
            RC(NC+I)  = RFW(I)
            ZC(NC+I)  = ZFW(I)
            PC(NC+I)  = CURFW(I)
            VC(NC+I)  = DFW(I)
            HC(NC+I)  = HFW(I)
            RC1(NC+I) = RFW1(I)
            ZC1(NC+I) = ZFW1(I)
            RC2(NC+I) = RFW2(I)
            ZC2(NC+I) = ZFW2(I)
            j=NCEQUI+I
            RES(j,j) = RESFW(I) * BBB 
            VOLK(NCEQUI+I) = 0.00D0   * BBB
         enddo
         NC     = NC     + NFW
         NCEQUI = NCEQUI + NFW
      endif

      call TREBP( NOUT, NTER, NINFW, NSEGBP,
     &            KDBP, RP1BP, ZP1BP, RP2BP, ZP2BP, RBPSEG, CBPSEG,
     &            NBP,  NTYBP,
     &            RBP,  ZBP,  DBP,   HBP,  RESBP,  CURBP,
     &            RBP1, ZBP1, RBP2,  ZBP2 )

      if( NBP.NE.0 )  then
         do I=1,NBP
            NTYPE(NC+I) = NTYBP(I)
            RC(NC+I)  = RBP(I)
            ZC(NC+I)  = ZBP(I)
            PC(NC+I)  = CURBP(I)
            VC(NC+I)  = DBP(I)
            HC(NC+I)  = HBP(I)
            RC1(NC+I) = RBP1(I)
            ZC1(NC+I) = ZBP1(I)
            RC2(NC+I) = RBP2(I)
            ZC2(NC+I) = ZBP2(I)
            j=NCEQUI+I
            RES(j,j) = RESbp(I) * BBB
            VOLK(NCEQUI+I) = 0.00D0   * BBB
         enddo
         NC     = NC     + NBP
         NCEQUI = NCEQUI + NBP
      endif

      call TREVV( NOUT, NTER, NINFW, NSEGVV,
     &            KDVV, RP1VV, ZP1VV, RP2VV, ZP2VV, RVVSEG, CVVSEG,
     &            NVV,  NTYVV,
     &            RVV,  ZVV,  DVV,   HVV,  RESVV, CURVV,
     &            RVV1, ZVV1, RVV2,  ZVV2 )

      if( NVV.NE.0 )  then
         do I=1,NVV
            NTYPE(NC+I) = NTYVV(I)
            RC(NC+I)  = RVV(I)
            ZC(NC+I)  = ZVV(I)
            PC(NC+I)  = CURVV(I)
            VC(NC+I)  = DVV(I)
            HC(NC+I)  = HVV(I)
            RC1(NC+I) = RVV1(I)
            ZC1(NC+I) = ZVV1(I)
            RC2(NC+I) = RVV2(I)
            ZC2(NC+I) = ZVV2(I)
            j=NCEQUI+I
            RES(j,j) = RESvv(I) * BBB 
            VOLK(NCEQUI+I) = 0.00D0 * BBB
         enddo
         NC     = NC     + NVV
         NCEQUI = NCEQUI + NVV
      endif

      write(fname,'(a,a)') TRIM(path), '/pascon.wr'
      open(1,file=fname,form='formatted')
         write(1,*) nc, ncpfc, nfw, nbp, nvv
         write(1,*) (rc(i), i=1,nc), (zc(i), i=1,nc)
      close(1)
      do i=1,NCEQUI
         VOLKP1(i) = VOLK(i)
      enddo

      return
      end subroutine CONDUC
