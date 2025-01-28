!----------------------------------------------------------------
!  SOLVES THE MATRIX CIRCUIT EQUATION
!----------------------------------------------------------------
!  INPUT DATE:
!  -----------
!  NLES, NREG, TSTEP0, TSTEP, SIGM, NJ,
!  VOLK(NJ), VOLKP1(NJ), RES(NJ),
!  PSPK(NJ), PSPKP1(NJ), CRPK(NJ), CRPKP1(NJ), CRPKP(NJ),
!  NOUT, NTER, KEYPRI,
!
!  PPIND(NJ,NJ)  - from 'comevl.inc'
!  DPSIDJ(NJ,NJ) - from 'comevl.inc' only for NREG=1
!
!  OUTPUT DATE:
!  -----------
!  CRPKP1(NJ) - solution array - currents
!  EREVE      - "nevyazka" of circuit equation solution
!               ( without DPSIDJ*(CRPKP1 - CRPKP) term )
!----------------------------------------------------------------

      subroutine EVSLV(NLES, NREG, TSTEP, SIGM, NJ, VOLK, VOLKP1, RES,
     &                 PSPK, PSPKP1, CRPK, CRPKP1, CRPKP,
     &                 EREVE, tstep2)

      use sp_parameters, only: njlim, nplim, npfc0, nnlim, nsp
      use comevl, only: ppind, dpsidj

      implicit none

      integer, intent(in) :: NLES, NREG, NJ
      real*8, intent(in) :: tstep2, tstep, sigm
      real*8, intent(in), dimension(njlim) :: VOLK, VOLKP1, PSPK,
     &    PSPKP1, CRPK, CRPKP1, CRPKP
      real*8, intent(in), dimension(njlim, njlim) :: RES

      real*8, intent(out) :: EREVE

      integer :: cmnd_dioh2s, cmnd_dioh2u, i, j, nesp, nflag, NN,
     &    npath, ia1, ia2, jaj
      integer :: IA(NJLIM+1), JA(NNLIM), ISP(NSP)
      integer, dimension(njlim) :: PP, P, IP
      real*8 :: BZZ, RS1, RS2
      real*8, dimension(njlim) :: BZ, A, RSP, B

      COMMON /SPA_MAT/ A, IA, JA
      COMMON /SPArsp/ RSP, PP, P, IP
      common /commanddioh/ cmnd_dioh2s,cmnd_dioh2u
	
      EQUIVALENCE (RSP(1), ISP(1))

      nflag = 0
      NN = NJ**2

      if (NLES == 0) then
          NPATH = 1
          do I=1, NJ
              PP(I) = I
              P(I)  = I
              IP(I) = I
          enddo

          IA(1) = 1
          do I=2, NJ
              IA(I) = IA(I-1) + NJ
          enddo
          IA(NJ+1) = NN + 1

          do I=1, NJ
              IA1 = IA(I)
              IA2 = IA(I+1) - 1
              do J=IA1, IA2
                  JAJ = J - IA1 + 1
                  JA(J) = JAJ
                  A(J) = PPIND(I, JAJ) + TSTEP*SIGM*RES(I, JAJ)
                  IF (NREG ==  1) A(J) = A(J) + DPSIDJ(I, JAJ)
              enddo
          enddo
      endif

      if (NLES == 1) then
          NPATH = 1
          do I=1, NJ
              IA1 = IA(I)
              IA2 = IA(I+1) - 1
              do J=IA1, IA2
                  JAJ   = JA(J)
                  A(J)  = PPIND(I, JAJ) + TSTEP*SIGM*RES(I, JAJ)
                  if (NREG ==  1) A(J) = A(J) + DPSIDJ(I, JAJ)
              enddo
          enddo
      endif

      if (NLES == 2) then
          NPATH = 3
      endif

      if ((NLES /= 0).AND.(NLES /= 1).AND.(NLES /= 2)) then
          WRITE(*, *) 'PARAMETER NLES = ',  NLES
          WRITE(*, *) 'IT IS WRONG. PROGRAM INTERRUPT'
          STOP
      endif

      do I=1, NJ
          RS1 = 0.d0
          RS2 = 0.d0
          do J=1, NJ
              RS1 = RS1 + CRPK(J)*(PPIND(I, J) -
     &            TSTEP*RES(I, j)*(1.D0-SIGM))
              IF (NREG == 1) RS2 = RS2 + CRPKP(J) * DPSIDJ(I, J)
          enddo

          B(I)  = TSTEP*(SIGM*VOLKP1(I)+(1.D0-SIGM)*VOLK(I)) +
     &            RS1 - (PSPKP1(I) - PSPK(I)) * TSTEP/TSTEP2

          if (NREG == 1) B(I) = B(I) + RS2
          BZ(I) = B(I)
      enddo

      call sdrvd(nj, p, ip, ia, ja, a, b, CRPKP1, nsp,
     *           isp, rsp, nesp, npath, nflag)

      EREVE = 0.0D0

      do I=1, NJ
          RS1 = 0.d0
          do J=1, NJ
              RS1 = RS1 + (CRPKP1(J) - CRPK(J)) * PPIND(I, J) +
     &        (SIGM*CRPKP1(J) + (1.D0 - SIGM)*CRPK(J))*TSTEP*res(I, J)
          enddo
          BZZ = TSTEP*(1.D0 - SIGM)* VOLK(I) + TSTEP*SIGM*VOLKP1(I) -
     &          RS1 - (PSPKP1(I) - PSPK(I)) * TSTEP / TSTEP2
          EREVE = max(EREVE,  ABS(BZZ))
      enddo

      if (NESP < 0) then
          WRITE(*, *) 'ATTENTION ESP(SDRVD) = ',  NESP
          WRITE(*, *) 'IT IS WRONG. PROGRAM INTERRUPT'
          STOP
      endif
      if (NFLAG /= 0) then
          WRITE(*, *) 'ATTENTION FLAG(SDRVD) = ',  NFLAG
          WRITE(*, *) 'IT IS WRONG. PROGRAM INTERRUPT'
          STOP
      endif

      return
      end subroutine EVSLV
