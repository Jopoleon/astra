      subroutine SSTEPON(KLUCH, k_auto, nstep, dt, time, voltpf,
     &    key_dmf)

      use durs_d_modul, only: n_tht, n_psi, igdf, nurs, i_eqdsk, i_betp, 
     &    keyctr, epsro, betplx, tokf, 
     &    b0, r0, rax, zax, psax, rxpnt, zxpnt,
     &    alf0, alf1, alf2, bet0, bet1, bet2
      use sp_parameters, only: nilim, njlim, nclim, nloopp, nprobp,
     &    twopi
      use iopath, only: path
      use comevl, only: nequi, pfceqw, rm_ef, zm_ef, rxpnt_ef, zxpnt_ef,
     &    tstep_ef, dt_ef, time_ef
      use keys, only: kstep, kpr, kastr, kastr2
      use e_nels, only: enels
      use ndmf, only: n_dmf

      implicit none

      integer, intent(in) :: kluch, k_auto, nstep, key_dmf
      real*8, intent(in) :: dt, time
      real*8, intent(in), dimension(*) :: voltpf

      integer :: i, j, l, ipsmk, eq_cmd, i_bsh, ngrid, nprob, nvv,
     &   kstop, k_step, kstepr, keypri, knel, knels, 
     &   n_ctrl, nctrl, numlim, it_dmf, nreg, nles,
     &   npro, nloop, nc, ncequi, ncpfc, nfw, nbp, ngav, ngav1, nursb,
     &   nfrpr1, nfrwr1, nout, nter, ninfw, ninev, ngra1, ngra2, numwr
      integer, dimension(nilim) :: necon
      integer, dimension(nclim) :: ntype

      real*8 :: timev, tstep, sigm, tstepr, tstart, tstop, time_fin,
     &   time_sta, bbb, qcen, betap0, z0cen, zcen0, zlold, 
     &   alw0, alw1, alw2, helinp, helout, diftok, platok, ftok, 
     &   psi_bnd, psi0_bnd, psicen, psiax, psiout, psi0_ax, pspl_av,
     &   psex_av, psidel, psibou, psax_sta, e_psi, erps, betpol, zli3,
     &   ereve, ztok_n, zpsim_n, beold, tokout,
     &   alp, alp_b, alpnew, up, f_wes, zpsim, ztok, flx_fi, difpsi,
     &   curmin, curmax, sgmcur, errcu1, errcu2,
     &   rx0, zx0, rmax0, zmax0, rm, zm, rxpold, zxpold, rxppr, zxppr,
     &   raxpr, zaxpr, rmaold, zmaold
      real*8, dimension(nilim) :: wecon
      real*8, dimension(nclim) :: PC, PSIP, VC, HC, RC, RC1, RC2,
     &   RC3, RC4, ZC, ZC1, ZC2, ZC3, ZC4
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, 
     &   pjkd, psk, pskp1, pskp, pskm1
      real*8, dimension(njlim, njlim) :: res
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      double precision :: dteqz
      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/ ncequi
      common /comhel/ helinp, helout
      common /comsta/ platok, eqdfn, i_bsh 
      common /com234/ betpol, tokout, psiout
      common /comst0/ RES, VOLK, VOLKP1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob
      common /comst2/ PSK, PSKP1, PSKP, PSKM1
      common /timingcmdipsmk/ ipsmk, eq_cmd, dteqz

      kstep = nstep
      timev = time
      tstep = dt

      if (ipsmk.lt.1) dteqz=dt !Efable

      BBB  = 1.d0 / twopi
      SIGM   = 1.d0  !0.50d0

      TSTEPR = TSTEP
      TSTART = 0.00d0
      TSTOP  = 20.0d0
      KSTOP  = 1000
      KNELS  = 500

      NFRPR1 = 5
      NFRWR1 = 0
      NOUT  = 17
      NTER  = 6
      NINFW = 7
      NINEV = 10
      ngra1 = 14
      ngra2 = 15

      if ( KLUCH .NE. 0 ) goto 1111

      numwr  = 0

!----------
! ITER, SOF

      KSTEPR = KSTEP
      KEYPRI = 1

!-----------------------------------------------------------------
! INPUT OF "BASIC" EQULIBRIUM PARAMETERS FROM FILE "durs.dat"

      if (k_auto.eq.0) goto 2005
      if (kastr.ne.1) then
         write(fname, '(a, a)') TRIM(path), '/durs.dat'
         open(1, file=fname)
            read(1, *) n_tht
            read(1, *) n_psi
            read(1, *) igdf
            read(1, *) epsro
            read(1, *) nurs
            read(1, *) keyctr
            read(1, *)    i_eqdsk
            read(1, '(A40)') eqdfn 
            read(1, *) i_betp
            read(1, *) betplx
            read(1, *) tokf
            read(1, *) psax
            read(1, *) b0, r0
            read(1, *) alf0
            read(1, *) alf1 ! P'=alf0*(1-(1-psi)**alf1)**alf2
            read(1, *) alf2
            read(1, *) bet0
            read(1, *) bet1 ! FF'=bet0*(1-(1-psi)**bet1)**bet2
            read(1, *) bet2
            read(1, *) alw0
            read(1, *) alw1 ! (ro*w**2)'=alw0*(1-(1-psi)**alw1)**alw2
            read(1, *) alw2
            read(1, *) rax
            read(1, *) zax
         close(1)
         if (i_eqdsk.eq.1) then
            call tab_efit(tokf, psax, eqdfn, rax, zax, b0, r0)
            nurs   = -3999
            i_betp = 0
         endif
         write(fname, '(a, a)') TRIM(path), '/durs_d.dat'
         open(1, file=fname, form='formatted')
            write(1, *) n_tht, n_psi, igdf, nurs, keyctr, i_eqdsk, 
     &                  i_betp
            write(1, *) epsro, betplx, tokf, psax, b0, r0, rax, zax
            write(1, *) alf0, alf1, alf2, bet0, bet1, bet2
         close(1)
      endif

!-----------------------------------------
! INPUT OF POSITIONS OF "PF_PROBE" POINTS:

      call PROPNT(NOUT, NTER, NINFW, NGRA1, NPRO, RPROb, ZPROb, FIPROb)

!----------------------------------------
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:
      call LOOPNT(NOUT, NTER, NINFW, NGRA1, NLOOp, RLOOp, ZLOOp)

!-----------------------------------------------------
! INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS
      call CONDUC(NC, NCEQUI, NCPFC, NFW, NBP, NVV, 
     &            RC, ZC, PC, VC, HC, NTYPE, 
     &            RC1, ZC1, RC2, ZC2, RC3, ZC3, RC4, ZC4, 
     &            RES, VOLK, VOLKP1, 
     &            NECON, WECON, 
     &            NOUT, NTER, NINFW, ngra1)

!----------------------------------------------------
! DEFINITION INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  "PPIND" from COMMON /PPIDPS/

      call L_MATR(NOUT, NTER, NC, NCPFC, NTYPE, RC,  ZC, VC, HC, 
     &            NECON, WECON)

!---------------------------------------------------
! Initial condition (currents) for circuit equations

      do L=1, NCEQUI
	 if ( L.LE.NEQUI ) then
            PJK(L)  = PFCEQW(L)
            PJKP(L) = PJK(L)
         else
            PJK(L)  = PC(NCPFC+L-NEQUI)
            PJKP(L) = PJK(L)
         endif
      enddo

      if (kastr2.eq.1) then
         open(32, file='currents.wr')
            read(32, *) nequi, ncequi
            read(32, *) (pjk(i), i=1, ncequi)
         close(32)
         pjk(nequi+1:ncequi)=0.d0
         do L=1, NCEQUI
            if ( L.LE.NEQUI ) then
               PFCEQW(L)  = PJK(L)
               PJKP(L) = PJK(L)
            else
               PC(NCPFC+L-NEQUI)=PJK(L)
               PJKP(L) = PJK(L)
            endif
         enddo
      endif

      write(fname, '(a, a)') TRIM(path), '/currents.wr'
      open(1, file=fname, form='formatted')
         write(1, *) NEQUI, NCEQUI
         write(1, *)(pjk(j), j=1, NCEQUI)
      close(1)

      write(fname, '(a, a)') TRIM(path), '/res_mat.wr'
      open(1, file=fname, form='formatted')
         write(1, *) ((res(i, j), j=1, ncequi), i=1, ncequi)
      close(1)

      HELOUT = HELINP
      call wrcoil(nc, ncpfc, rc, zc, pc, necon, wecon)
      ngrid=1

      call auto(rc, zc, pc, nc, nstep, ngrid, 
     &          rc1, zc1, rc2, zc2, rc3, zc3, rc4, zc4, 
     &          ntype, necon, wecon )

 2005 continue

      if (kastr.eq.0 .AnD. i_eqdsk.eq.0) then 
         write(fname, '(a, a)') TRIM(path), '/inpol.dat'
         open(1, file=fname, form='formatted')
            read(1, *) i_bsh
         close(1)    
      else
         i_bsh=1
      endif

      if (kastr.eq.0 ) then 
         write(fname, '(a, a)') TRIM(path), '/durs_d.dat'
         open(1, file=fname, form='formatted')
            read(1, *) n_tht, n_psi, igdf, nurs, keyctr, i_eqdsk, i_betp
            read(1, *) epsro, betplx, tokf, psax, b0, r0, rax, zax
            read(1, *) alf0, alf1, alf2, bet0, bet1, bet2
         close(1)
      endif

      call noauto
      call rd_ppind
      call rd_prob( NPROb, RPROb, ZPROb, FIPROb )
      call rd_loop( NLOOp, RLOOp, ZLOOp )

      write(fname, '(a, a)') TRIM(path), '/currents.wr'
      open(1, file=fname, form='formatted')
         read(1, *) nequi, ncequi
         read(1, *)(pjk(j), j=1, ncequi)
      close(1)

      write(fname, '(a, a)') TRIM(path), '/res_mat.wr'
      open(1, file=fname, form='formatted')
         read(1, *) ((res(i, j), j=1, ncequi), i=1, ncequi)
      close(1)

      do L=1, NCEQUI
         PJKP(L) = PJK(L)
         PJKP1(L) = PJK(L)
      enddo

! voltage for zero time-level

      do i=1, nequi
         VOLK(i) = voltpf(i)*BBB
         VOLKp1(i) = VOLK(i)
      enddo
      do i=nequi+1, NCEQUI
         VOLK(i) = 0.0d0
         VOLKp1(i) = VOLK(i)
      enddo

      platok = tokf
      kstep  = 0
      k_step=0
      ngav=keyctr

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         ngav, k_step, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      ngrid=1
      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf

      call eq_0(pjk, psip, ncequi, kstep, ngrid, 
     &          alf0, alf1, alf2, bet0, bet1, bet2, 
     &          betplx, ngav1, 
     &          ftok, tokout, psiax, nursb, 
     &          psi_bnd, alp_b, rax, zax, n_ctrl, b0, r0 )

      call rdexf(ncequi)

      call eq(pjk, psk, ncequi, kstep, ngrid, 
     &        alf0, alf1, alf2, bet0, bet1, bet2, 
     &        betpol, betplx, ngav1, 
     &        tokout, psiout, 
     &        nursb, psi_bnd, alp_b, rax, zax )

      do L=1, NCEQUI
         PSKP1(L) = PSK(L)
      enddo

      call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &            rm, zm, rx0, zx0)

! --- FOR diagnostic

      RMAX0  = RM
      ZMAX0  = ZM
      ZCEN0  = Z0CEN
      BEOLD  = BETPOL
      ZLOLD  = ZLI3
      RAXPR  = RM
      ZAXPR  = ZM
      RMAOLD = RM
      ZMAOLD = ZM
      RXPPR  = RX0
      ZXPPR  = ZX0
      RXPOLD = RX0
      ZXPOLD = ZX0

!--------------------------------------------------------
! Definition of input parameters from "basic" equilibrium
      BETPLX = BETPOL
      BETAP0 = BETPOL
      FTOK   = TOKOUT
      PSIAX  = PSIOUT
      HELINP = HELOUT

!----------------------------------------------------
! For PSI_axis value time scenario (for NGAV1=2 or 3)
      time_sta = TSTART
      time_fin = TSTOP
      psax_sta = PSIOUT
      PSIBOU = UP
      PSIDEL = PSIOUT - PSIBOU

! If we need to compute the "basic" equilibrium only -> KSTOP=0

      if (KSTOP.EQ.0) then
         return
      endif

      i_bsh  = -1
      i_eqdsk  = 0
      kstep  =  0
      keyctr =  0
      e_psi  =  0.d0

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      call get_par(psi_bnd)
      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      ztok_n  = platok
      zpsim_n = psax
      rxpnt = rx0
      zxpnt = zx0

      return

!-------------------

 1111 continue

      keyctr = key_dmf
      if (keyctr .eq. 0) platok = tokf
      call get_par(psi_bnd)
      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      ztok_n  = platok
      zpsim_n = psax
      NGAV1 = 0
      KEYPRI = 0 ! Parameter for printing

! Definition of input equil. paramters: BETPLX, PSIAX, FTOK, HELINP

      if (NGAV1.EQ.0 .OR. NGAV1.EQ.2 .OR. NGAV1.EQ.4) then
         BETPLX = BETPOL
      endif
      if ( NGAV1.EQ.0 .OR. NGAV1.EQ.1 ) then
         PSIAX  = PSIOUT
         HELINP = HELOUT
      endif
      if ( NGAV1.EQ.2 .OR. NGAV1.EQ.3 ) then
         FTOK   = TOKOUT
         HELINP = HELOUT
      endif
      if ( NGAV1.EQ.4 .OR. NGAV1.EQ.5 ) then
         FTOK   = TOKOUT
         PSIAX  = PSIOUT
      endif

!-----------------------------------------------------
! Computing of "VOLKP1"  voltage  for CIRCUIT EQUATION 
! for new time level KSTEP

      do i=1, nequi
         VOLKP1(i) = voltpf(i)*BBB
      enddo
      do i=nequi+1, NCEQUI
         VOLKP1(i) = 0.0d0
      enddo

      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)
      call get_flfi(flx_fi)

      psi0_bnd=pspl_av+psex_av
      psi0_ax=pspl_av+psex_av+psax
      ztok=platok
      zpsim=psax
      i_bsh=-1

      if (keyctr.eq.0) n_dmf=1

      do it_dmf=1, n_dmf ! it_dmf <<< mag.field diffusion iter. loop
         if (kpr.eq.1) print *, ' it_dfm==', it_dmf
         if (keyctr.ne.0) then
            call eqb(alf0, alf1, alf2, bet0, bet1, bet2, 
     &               alw0, alw1, alw2, betplx, i_betp, keyctr,
     &               kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &               n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &               psi_bnd, psi0_bnd )
         endif
         f_wes=1.0d0
         if (it_dmf.eq.1)then
            ftok=platok
         else
            ftok=f_wes*platok+(1.d0-f_wes)*ftok
         endif
         if ( KSTEP.EQ.1 ) then
            NGRID  = 1
            do L=1, NCEQUI
               PSKP1(L) = PSK(L)
            enddo
            NREG = 0
            NLES = 0
            call EVSLV(NLES,  NREG, TSTEP, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSK,   PSKP1, PJK,   PJKP1, PJKP, 
     &                 NOUT,  NTER, KEYPRI, EREVE, tstep )
         else
!----------------------------------------------------------------
!  The initial approximation for the case of closed-loop evolution

            do L=1, NCEQUI
               PJKP(L) = PJK(L)
            enddo
            NREG = 0
            NLES = 1
            call EVSLV(NLES,  NREG, TSTEPR, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSKM1, PSK, PJK,   PJKP1, PJKP, 
     &                 NOUT,  NTER, KEYPRI, EREVE , tstep)

         endif

	 if (eq_cmd.eq.1.and.ipsmk.ge.1) then
            PSKM1 = PSK
            PSK   = PSKP1
            call EQ_AX(pjkp1, PSkp1, NCequi,  KSTEP, NGRID, 
     &                 ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                 BETPOL,  BETPLX,  ZLI3, 
     &                 NGAV1, 
     &                 FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                 ENELS, ERPS, 
     &                 psi_bnd, alp_b, rax, zax, 0)

            call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &                  rm, zm, rx0, zx0)
	 endif

!---------------------
! Start tieration loop

         do KNEL=1, KNELS

            do L=1, NCEQUI
               PSKP(L) = PSKP1(L)
               PJKP(L) = PJKP1(L)
            enddo
            NREG = 0
            NLES = 2
            call EVSLV(NLES,  NREG, TSTEP, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSK,   PSKP1, PJK,   PJKP1, PJKP, 
     &                 NOUT,  NTER, KEYPRI, EREVE , tstep)
            call DIFFER(PJKP, PJKP1, NCEQUI, ERRCU1, ERRCU2, 
     &                  CURMAX, CURMIN, NOUT, NTER)
            SGMCUR = 0.75D0
            do L=1, NCEQUI
               PJKP1(L) = SGMCUR*PJKP1(L) + (1.0D0 - SGMCUR)*PJKP(L)
            enddo

            if (kpr.ge.-1) then !efable for skip fsim
            call EQ_AX(pjkp1, PSkp1, NCequi, KSTEP, NGRID, 
     &                 ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                 BETPOL, BETPLX, ZLI3, 
     &                 NGAV1,
     &                 FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                 ENELS, ERPS, 
     &                 psi_bnd, alp_b, rax, zax, 0)
            if (kpr.gt.0) write(*, *) 'stepon:EQ_AX done, erru=', ERPS
               call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up,
     &                     rm, zm, rx0, zx0)
            endif     
            if (NGAV1.EQ.2 .OR. NGAV1.EQ.3)  FTOK = TOKOUT
            if (NGAV1.EQ.4 .OR. NGAV1.EQ.5)  FTOK = TOKOUT
            if (kpr.gt.0) write(*, *) 'stepon:DIFTIM, erru=', ERPS
            if (kpr.lt.-1) erps=ereve !for skip EFable fsim
            if (ERPS .LT. ENELS) EXIT
            if (KNEL .EQ. KNELS) then
               KEYPRI = 1
               write(*, *) 'spider, sstepon: no covergence ' 
               write(*, *) 'KNEL .GE. KNELS  ', KNEL, KNELS 
            endif
         enddo

         if (kpr.ge.0) then
            write(fname, '(a, a)') TRIM(path), '/knel_iters.wr'
            open(1, file=fname)
               write(1, *) KNEL
            close(1)
         endif

         diftok=dabs(platok-ztok)/(dabs(platok-ztok_n)+1.d-8)
         difpsi=dabs(psax-zpsim)/(dabs(psax-zpsim_n)+1.d-8)

         if ((diftok.lt.5.0d-3 .AnD. it_dmf.gt.0) .OR. keyctr.eq.0) then
            call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1,
     &               alw2, betplx, i_betp, 
     &               keyctr, kstep, platok, rax, zax, b0, r0, psax,
     &               igdf, n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &               psi_bnd, psi0_bnd)
            EXIT
         endif
         ztok  = platok
         zpsim = psax
      enddo

      ztok_n=platok
      zpsim_n=psax
      PSIBOU = UP
      PSIDEL = PSIOUT - PSIBOU

      RM_EF = RM
      ZM_EF = ZM
      TSTEP_EF = TSTEP
      DT_EF = DT
      TIME_EF = TIME
      RXPNT_EF = 0.
      ZXPNT_EF = 0.
      rxpnt = rx0
      zxpnt = zx0

      return
      end subroutine SSTEPON

!---------------------------------------------------------------------
      subroutine SSTEPON_BKDW(KLUCH, k_auto, nstep, dt, time, voltpf)

      use durs_d_modul, only: n_tht, n_psi, igdf, nurs, i_eqdsk, i_betp,
     &    keyctr, betplx, tokf, b0, r0, epsro,
     &    rax, zax, psax, rxpnt, zxpnt,
     &    alf0, alf1, alf2, bet0, bet1, bet2
      use sp_parameters, only: nilim, njlim, nclim, nloopp, nprobp,
     &    twopi
      use iopath, only: path
      use comevl, only: nequi, rm_ef, zm_ef, tstep_ef, dt_ef, time_ef,
     &    rxpnt_ef, zxpnt_ef, pfceqw
      use keys, only: kstep, kpr, kastr
      use e_nels, only: enels
      use ndmf, only: n_dmf

      implicit none

      integer, intent(in) :: kluch, k_auto, nstep
      real*8, intent(in) :: dt, time
      real*8, intent(in), dimension(*) :: voltpf

      integer :: i, j, l, ipsmk, eq_cmd, i_bsh, ngrid, nprob, nvv,
     &   kstop, k_step, kstepr,  keypri, knel, knels, 
     &   n_ctrl, nctrl, numlim, it_dmf, nreg, nles,
     &   npro, nloop, nc, ncequi, ncpfc, nfw, nbp, ngav, ngav1, nursb,
     &   nfrpr1, nfrwr1, nout, nter, ninfw, ninev, ngra1, ngra2, numwr,
     &   j_calll, ispid_contour, res_trigts06
      integer, dimension(nilim) :: necon
      integer, dimension(nclim) :: ntype

      real*8 :: timev, tstep, sigm, tstepr, tstart, tstop, time_fin,
     &   time_sta, bbb, qcen, betap0, z0cen, zcen0, zlold, 
     &   alw0, alw1, alw2, helinp, helout, diftok, platok, ftok, 
     &   psi_bnd, psi0_bnd, psicen, psiax, psiout, psi0_ax, pspl_av,
     &   psex_av, psidel, psibou, psax_sta, e_psi, erps, betpol, zli3,
     &   ereve, ztok_n, zpsim_n, beold, tokout,
     &   alp, alp_b, alpnew, up, f_wes, zpsim, ztok, flx_fi, difpsi,
     &   curmin, curmax, sgmcur, errcu1, errcu2,
     &   rx0, zx0, rmax0, zmax0, rm, zm, rxpold, zxpold, rxppr, zxppr,
     &   raxpr, zaxpr, rmaold, zmaold,
     &   cmnd_dioh2s, cmnd_dioh2u
      real*8, dimension(nilim) :: wecon
      real*8, dimension(nclim) :: PC, PSIP, VC, HC, RC, RC1, RC2,
     &   RC3, RC4, ZC, ZC1, ZC2, ZC3, ZC4
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, 
     &   pjkd, psk, pskp1, pskp, pskm1
      real*8, dimension(njlim, njlim) :: res
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      double precision :: dteqz, dt_nowww, resres_oh6
      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/  ncequi
      common /comhel/ helinp, helout
      common /comsta/ platok, eqdfn, i_bsh
      common /com234/ betpol, tokout, psiout
      common /comst0/ RES, VOLK, VOLKP1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob
      common /comst2/ PSK, PSKP1, PSKP, PSKM1
      common /timingcmdipsmk/ ipsmk, eq_cmd, dteqz
      common /commanddioh/ cmnd_dioh2s, cmnd_dioh2u
      common /spidcontouring/ ispid_contour

      data j_calll/0/
      save dt_nowww, j_calll, resres_oh6, res_trigts06

      kstep = nstep
      timev = time
      tstep = dt

      if (j_calll.eq.0) then
         dt_nowww=tstep		
         j_calll=1
      endif

      if (ipsmk.lt.1) dteqz=dt !Efable

      BBB  = 1.d0 / twopi
      SIGM   = 1.d0  !0.50d0
      TSTEPR = TSTEP
      TSTART = 0.00d0
      TSTOP  = 20.0d0
      KSTOP  = 1000
      KNELS  = 500
      NFRPR1 = 5
      NFRWR1 = 0
      NOUT  = 17
      NTER  = 6
      NINFW = 7
      NINEV = 10
      ngra1 = 14
      ngra2 = 15

      if ( KLUCH .NE. 0 ) goto 1111

      numwr  = 0

!----------
! ITER, SOF

      dt_nowww=TSTEP

      KSTEPR = KSTEP
      KEYPRI = 1

!------------------------------------------------------------
! INPUT OF "BASIC" EQULIBRIUM PARAMETERS FROM FILE "durs.dat"

      if (k_auto.eq.0) goto 2005

      if (kastr.ne.1) then
         write(fname, '(a, a)') TRIM(path), '/durs.dat'
         open(1, file=fname)
            read(1, *) n_tht
            read(1, *) n_psi
            read(1, *) igdf
            read(1, *) epsro
            read(1, *) nurs
            read(1, *) keyctr
            read(1, *)    i_eqdsk
            read(1, '(A40)') eqdfn
            read(1, *) i_betp
            read(1, *) betplx
            read(1, *) tokf
            read(1, *) psax
            read(1, *) b0, r0
            read(1, *) alf0
            read(1, *) alf1 ! P'=alf0*(1-(1-psi)**alf1)**alf2
            read(1, *) alf2
            read(1, *) bet0
            read(1, *) bet1 ! FF'=bet0*(1-(1-psi)**bet1)**bet2
            read(1, *) bet2
            read(1, *) alw0
            read(1, *) alw1 ! (ro*w**2)'=alw0*(1-(1-psi)**alw1)**alw2
            read(1, *) alw2
            read(1, *) rax
            read(1, *) zax
         close(1)

         if (i_eqdsk.eq.1) then
            call tab_efit(tokf, psax, eqdfn, rax, zax, b0, r0)
            nurs   = -3999
            i_betp = 0
         endif

         write(fname, '(a, a)') TRIM(path), '/durs_d.dat'
         open(1, file=fname, form='formatted')
            write(1, *) n_tht, n_psi, igdf, nurs, keyctr, i_eqdsk, 
     &         i_betp
            write(1, *) epsro, betplx, tokf, psax, b0, r0, rax, zax
            write(1, *) alf0, alf1, alf2, bet0, bet1, bet2
         close(1)
      endif

!-----------------------------------------
! INPUT OF POSITIONS OF "PF_PROBE" POINTS:
      call PROPNT(NOUT, NTER, NINFW, NGRA1,
     &            NPRO, RPROb, ZPROb, FIPROb)

!----------------------------------------
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:
      call LOOPNT(NOUT, NTER, NINFW, NGRA1, 
     &            NLOOp, RLOOp, ZLOOp)

!-----------------------------------------------------
! INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS

      call CONDUC(NC, NCEQUI, NCPFC, NFW, NBP, NVV, 
     &            RC, ZC, PC, VC, HC, NTYPE, 
     &            RC1, ZC1, RC2, ZC2, RC3, ZC3, RC4, ZC4, 
     &            RES, VOLK, VOLKP1, 
     &            NECON, WECON, 
     &            NOUT, NTER, NINFW, ngra1)

!-----------------------------------------------------
! DEFINITION INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  "PPIND" from COMMON /PPIDPS/

      call L_MATR(NOUT, NTER, NC, NCPFC, NTYPE, RC,  ZC, VC, HC,
     &            NECON, WECON )

!---------------------------------------------------
! Initial condition (currents) for circuit equations

      do L=1, NCEQUI
         if ( L.LE.NEQUI ) then
            PFCEQW(L)  = PJK(L)
         else
            PC(NCPFC+L-NEQUI)=PJK(L)
         endif
      enddo
      do L=1, NCEQUI
         if ( L.LE.NEQUI ) then
            PJKP(L) = PJK(L)
         else
            PJKP(L) = PJK(L)
         endif
      enddo

      write(fname, '(a, a)') TRIM(path), '/currents.wr'
      open(1, file=fname, form='formatted')
         write(1, *) NEQUI, NCEQUI
         write(1, *)(pjk(j), j=1, NCEQUI)
      close(1)

      write(fname, '(a, a)') TRIM(path), '/res_mat.wr'
      open(1, file=fname, form='formatted')
         write(1, *) ((res(i, j), j=1, ncequi), i=1, ncequi)
      close(1)

      HELOUT = HELINP

      call wrcoil(nc, ncpfc, rc, zc, pc, necon, wecon)

      ngrid=1

      call auto(rc, zc, pc, nc, nstep, ngrid, 
     &          rc1, zc1, rc2, zc2, rc3, zc3, rc4, zc4, 
     &          ntype, necon, wecon )

 2005 continue

      if (kastr.eq.0 .AnD. i_eqdsk.eq.0) then 
         write(fname, '(a, a)') TRIM(path), '/inpol.dat'
         open(1, file=fname, form='formatted')
            read(1, *) i_bsh
         close(1)    
         else
            i_bsh=1
         endif

         if (kastr.eq.0 ) then 
         write(fname, '(a, a)') TRIM(path), '/durs_d.dat'
         open(1, file=fname, form='formatted')
            read(1, *) n_tht, n_psi, igdf, nurs, keyctr, i_eqdsk, i_betp
            read(1, *) epsro, betplx, tokf, psax, b0, r0, rax, zax
            read(1, *) alf0, alf1, alf2, bet0, bet1, bet2
         close(1)
      endif

      call noauto
      call rd_ppind
      call rd_prob( NPROb, RPROb, ZPROb, FIPROb )
      call rd_loop( NLOOp, RLOOp, ZLOOp )

      write(fname, '(a, a)') TRIM(path), '/currents.wr'
      open(1, file=fname, form='formatted')
         read(1, *) nequi, ncequi
         read(1, *)(pjk(j), j=1, ncequi)
      close(1)

      write(fname, '(a, a)') TRIM(path), '/res_mat.wr'
      open(1, file=fname, form='formatted')
         read(1, *) ((res(i, j), j=1, ncequi), i=1, ncequi)
      close(1)

      do L=1, NCEQUI
         PJKP(L) = PJK(L)
         PJKP1(L) = PJK(L)
      enddo

! Voltage for zero time-level

      do i=1, nequi
         VOLK(i) = voltpf(i)*BBB
         VOLKp1(i) = VOLK(i)
      enddo
      do i=nequi+1, NCEQUI
         VOLK(i) = 0.0d0
         VOLKp1(i) = VOLK(i)
      enddo

      platok = tokf
      kstep  = 0
      k_step=0
      ngav=keyctr

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         ngav, k_step, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      ngrid=1
      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf

      call eq_0(pjk, psip, ncequi, kstep, ngrid, 
     &          alf0, alf1, alf2, bet0, bet1, bet2, 
     &          betplx, ngav1, 
     &          ftok, tokout, psiax, nursb, 
     &          psi_bnd, alp_b, rax, zax, n_ctrl, b0, r0)
      call rdexf(ncequi)
      call eq(pjk, psk, ncequi, kstep, ngrid, 
     &        alf0, alf1, alf2, bet0, bet1, bet2, 
     &        betpol, betplx, ngav1, 
     &        tokout, psiout, 
     &        nursb, psi_bnd, alp_b, rax, zax)

      do L=1, NCEQUI
         PSKP1(L) = PSK(L)
      enddo

      call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &            rm, zm, rx0, zx0)

! FOR diagnostic

      RMAX0  = RM
      ZMAX0  = ZM
      ZCEN0  = Z0CEN
      BEOLD  = BETPOL
      ZLOLD  = ZLI3
      RAXPR  = RM
      ZAXPR  = ZM
      RMAOLD = RM
      ZMAOLD = ZM
      RXPPR  = RX0
      ZXPPR  = ZX0
      RXPOLD = RX0
      ZXPOLD = ZX0

!--------------------------------------------------------
! Definition of input parameters from "basic" equilibrium
      BETPLX = BETPOL
      BETAP0 = BETPOL
      FTOK   = TOKOUT
      PSIAX  = PSIOUT
      HELINP = HELOUT

!-----------------------------------------------------
!  For PSI_axis value time scenario (for NGAV1=2 or 3)

      time_sta = TSTART
      time_fin = TSTOP
      psax_sta = PSIOUT
      PSIBOU = UP
      PSIDEL = PSIOUT - PSIBOU

! If we need to compute the "basic" equilibrium only -> KSTOP=0

      if ( KSTOP.EQ.0 ) then
         return
      endif

      i_bsh  = -1
      i_eqdsk  = 0
      kstep  =  0
      keyctr =  0
      e_psi  =  0.d0

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      call get_par(psi_bnd)
      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      ztok_n  = platok
      zpsim_n = psax
      rxpnt = rx0
      zxpnt = zx0

      return

 1111 continue

      keyctr = 0
      if (keyctr .eq. 0) platok = tokf
      call get_par(psi_bnd)
      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      ztok_n  = platok
      zpsim_n = psax

      NGAV1 = 0
      KEYPRI = 0 ! Parameter for printing

!-------------------------------------------------------------------
! DEFINITION OF INPUT EQUIL. PARAMETERS: BETPLX, PSIAX, FTOK, HELINP

      if (NGAV1.EQ.0 .OR. NGAV1.EQ.2 .OR. NGAV1.EQ.4) BETPLX = BETPOL
      if ( NGAV1.EQ.0 .OR. NGAV1.EQ.1 ) then
         PSIAX  = PSIOUT
         HELINP = HELOUT
      endif
      if ( NGAV1.EQ.2 .OR. NGAV1.EQ.3 ) then
         FTOK   = TOKOUT
         HELINP = HELOUT
      endif
      if ( NGAV1.EQ.4 .OR. NGAV1.EQ.5 ) then
         FTOK   = TOKOUT
         PSIAX  = PSIOUT
      endif

!-----------------------------------------------------
! Computing of "VOLKP1"  voltage  for CIRCUIT EQUATION 
! for new time level KSTEP

      do i=1, nequi
         VOLKP1(i) = voltpf(i)*BBB
      enddo
      do i=nequi+1, NCEQUI
         VOLKP1(i) = 0.0d0
      enddo

      if (res_trigts06.eq.1) then
         volkp1(1)=volkp1(1)-resres_oh6*pjkp(1)
      	 volkp1(2)=volkp1(2)-resres_oh6*pjkp(1)
      	 volkp1(3)=volkp1(3)-resres_oh6*pjkp(1)
      endif

      call get_flfi(flx_fi)

      psi0_ax=pspl_av+psex_av+psax
      ztok=platok
      zpsim=psax
      i_bsh=-1

      if (keyctr.eq.0) n_dmf=1

      do it_dmf=1, n_dmf ! it_dmf <<< mag.field diffusion iter. loop

         if (kpr.eq.1)print *, ' it_dfm==', it_dmf

         f_wes=1.0d0
         if (it_dmf.eq.1)then
            ftok=platok
         else
            ftok=f_wes*platok+(1.d0-f_wes)*ftok
         endif

         if ( KSTEP.EQ.1 ) then
            NGRID  = 1
            do L=1, NCEQUI
               PSKP1(L) = PSK(L)
            enddo

!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
      	          pjkp1(2)=pjkp1(1)	
      	          pjkp(2)=pjkp(1)	
      	          pjk(2)=pjk(1)
      	          pjkp1(3)=pjkp1(1)	
      	          pjkp(3)=pjkp(1)	
      	          pjk(3)=pjk(1)	
               endif
               if (cmnd_dioh2u.eq.0) then
      	          pjkp1(3)=pjkp1(2)	
      	          pjkp(3)=pjkp(2)	
      	          pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end

            if (ipsmk.ge.1) then
      	       PSKM1= PSK    
      	       PSK   = PSKP1
            endif
            do L=1, NCEQUI
               PSKP1(L) = PSK(L)
            enddo
            NREG = 0
            NLES = 0
            call EVSLV(NLES,  NREG, TSTEP, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSK,   PSKP1, PJK,   PJKP1, PJKP, 
     &                 NOUT,  NTER, KEYPRI, EREVE, dteqz)

!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
      	          pjkp1(2)=pjkp1(1)	
      	          pjkp(2)=pjkp(1)	
      	          pjk(2)=pjk(1)	

      	          pjkp1(3)=pjkp1(1)	
      	          pjkp(3)=pjkp(1)	
      	          pjk(3)=pjk(1)	
               endif

               if (cmnd_dioh2u.eq.0) then
      	          pjkp1(3)=pjkp1(2)	
      	          pjkp(3)=pjkp(2)	
      	          pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end
         endif    !  if ( KSTEP.EQ.1 )

!-------------------------------------------------
         if (eq_cmd.eq.1.and.ipsmk.ge.1.and.kpr.eq.-2) then
            PSKM1 = PSK
            PSK   = PSKP1
            call EQ_AX0(pjk, PSkp1, NCequi, KSTEP, NGRID, 
     &                  ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                  BETPLX, 
     &                  NGAV1, 
     &                  FTOK, PSIAX, 
     &                  ENELS,
     &                  psi_bnd)
      	    if (ispid_contour.eq.0) then
               call eqb(alf0, alf1, alf2, 
     &             bet0, bet1, bet2, alw0, alw1, alw2, 
     &             betplx, i_betp, 
     &             keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &             n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &             psi_bnd, psi0_bnd)
            else
               call eqb_contour(alf0, alf1, alf2, 
     &             bet0, bet1, bet2, alw0, alw1, alw2, 
     &             betplx, i_betp, 
     &             keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &             n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &             psi_bnd, psi0_bnd)
            endif
            call EQ_AX2( pjk, PSkp1, NCequi,  KSTEP, NGRID, 
     &         ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &         BETPOL, BETPLX,
     &         NGAV1, 
     &         FTOK, TOKOUT, PSIAX, PSIOUT, 
     &         ENELS, ERPS, 
     &         psi_bnd, alp_b, rax, zax, 0)

            call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &                  rm, zm, rx0, zx0)
         endif
         if (eq_cmd.eq.1.and.ipsmk.ge.1.and.kpr.ge.-1) then
            PSKM1 = PSK
            PSK   = PSKP1
            call EQ_AX(pjk, PSkp1, NCequi,  KSTEP, NGRID, 
     &                 ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                 BETPOL,  BETPLX,  ZLI3, 
     &                 NGAV1, 
     &                 FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                 ENELS, ERPS, 
     &                 psi_bnd, alp_b, rax, zax, 0)

            call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &                  rm, zm, rx0, zx0)
         endif
         if (ipsmk.lt.1) then
            call EQ_AX(pjk, PSkp1, NCequi,  KSTEP, NGRID, 
     &                 ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                 BETPOL,  BETPLX,  ZLI3, 
     &                 NGAV1, 
     &                 FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                 ENELS, ERPS, 
     &                 psi_bnd, alp_b, rax, zax, 0)

            call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &                  rm, zm, rx0, zx0)
         endif

!----------------------------
! Start of the iteration loop 

         if (KSTEP.NE.1) then

!----------------------------------------------------------------
! Initial approximation for the case of closed-loop evolution
            do L=1, NCEQUI
               PJKP(L) = PJK(L)
            enddo

            if (dt_nowww.eq.TSTEP) then
               NREG = 0
               NLES = 1
            else
               NREG = 0
               NLES = 0
            endif
           
!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
      	          pjkp1(2)=pjkp1(1)
      	          pjkp(2)=pjkp(1)
      	          pjk(2)=pjk(1)
                  pjkp1(3)=pjkp1(1)
     	          pjkp(3)=pjkp(1)
                  pjk(3)=pjk(1)
               endif
               if (cmnd_dioh2u.eq.0) then
      	          pjkp1(3)=pjkp1(2)	
      	          pjkp(3)=pjkp(2)	
      	          pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end

            call EVSLV(NLES, NREG, TSTEPR, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSK, PSKP1, PJK, PJKP1, PJKP, 
     &                 NOUT, NTER, KEYPRI, EREVE, dteqz)

!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
      	          pjkp1(2)=pjkp1(1)
      	          pjkp(2)=pjkp(1)
      	          pjk(2)=pjk(1)
      	          pjkp1(3)=pjkp1(1)	
      	          pjkp(3)=pjkp(1)	
      	          pjk(3)=pjk(1)	
               endif
               if (cmnd_dioh2u.eq.0) then
                  pjkp1(3)=pjkp1(2)	
      	          pjkp(3)=pjkp(2)	
                  pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end

         endif  ! if ( KSTEP.NE.1 )

         do KNEL=1, KNELS

            do L=1, NCEQUI
               PSKP(L) = PSKP1(L)
               PJKP(L) = PJKP1(L)
            enddo

            if (dt_nowww.eq.TSTEP) then
               NREG = 0
               NLES = 2
            else
               NREG = 0
               NLES = 0
               dt_nowww=TSTEP
            endif
!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
      	          pjkp1(2)=pjkp1(1)	
      	          pjkp(2)=pjkp(1)	
      	          pjk(2)=pjk(1)	
      	          pjkp1(3)=pjkp1(1)	
      	          pjkp(3)=pjkp(1)	
      	          pjk(3)=pjk(1)	
               endif

               if (cmnd_dioh2u.eq.0) then
                  pjkp1(3)=pjkp1(2)	
                  pjkp(3)=pjkp(2)	
      	          pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end

            call EVSLV(NLES, NREG, TSTEP, TSTEP, SIGM, 
     &                 NCEQUI, VOLK, VOLKP1, RES, 
     &                 PSK, PSKP1, PJK, PJKP1, PJKP, 
     &                 NOUT, NTER, KEYPRI, EREVE, dteqz)

            call DIFFER(PJKP, PJKP1, NCEQUI, ERRCU1, ERRCU2, 
     &                  CURMAX, CURMIN, NOUT, NTER)

            SGMCUR = 0.25D0

            do L=1, NCEQUI
               PJKP1(L) = SGMCUR*PJKP1(L) + (1.0D0 - SGMCUR)*PJKP(L)
            enddo

!Efable this block below is to keep dioh2u and dioh2s to zero for some time
            if (ipsmk.ge.1) then
               if (cmnd_dioh2s.eq.0) then
                  pjkp1(2)=pjkp1(1)	
      	          pjkp(2)=pjkp(1)	
      	          pjk(2)=pjk(1)	
       	          pjkp1(3)=pjkp1(1)	
      	          pjkp(3)=pjkp(1)	
      	          pjk(3)=pjk(1)	
               endif
               if (cmnd_dioh2u.eq.0) then
      	          pjkp1(3)=pjkp1(2)	
      	          pjkp(3)=pjkp(2)	
      	          pjk(3)=pjk(2)	
               endif
            endif
!Efable dioh block end

            if (kpr.ge.-1) then !efable for skip fsim
               if (eq_cmd.eq.1.and.ipsmk.ge.1) then
      
                  call EQ_AX(pjkp1, PSkp1, NCequi,  KSTEP, NGRID, 
     &                       ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                       BETPOL, BETPLX, ZLI3, 
     &                       NGAV1, 
     &                       FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                       ENELS, ERPS, 
     &                       psi_bnd, alp_b, rax, zax, 0)
                  if (kpr.gt.0) write(*, *) 'stepon:EQ_AX done, erru=', 
     &               ERPS
                  call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim,
     &                        up, rm, zm, rx0, zx0)
               endif
            endif     
            if (kpr.ge.-1) then !efable for skip fsim
      	       if (ipsmk.lt.1) then
                  call EQ_AX(pjkp1, PSkp1, NCequi, KSTEP, NGRID, 
     &                       ALF0, ALF1, ALF2, BET0, BET1, BET2, 
     &                       BETPOL, BETPLX, ZLI3, 
     &                       NGAV1, 
     &                       FTOK, TOKOUT, PSIAX, PSIOUT, 
     &                       ENELS, ERPS,
     &                       psi_bnd, alp_b, rax, zax, 0)

                  if (kpr.gt.0) write(*, *) 'stepon:EQ_AX done, erru=', 
     &               ERPS
                  call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, 
     &                        up, rm, zm, rx0, zx0)
               endif
            endif     

            if (NGAV1.EQ.2 .OR. NGAV1.EQ.3) FTOK = TOKOUT
            if (NGAV1.EQ.4 .OR. NGAV1.EQ.5) FTOK = TOKOUT
            if (kpr.gt.0) write(*, *) 'stepon:DIFTIM, erru=', ERPS
            if (kpr.lt.-1) erps=ereve !for skip EFable fsim
            if (ERPS .LT. ENELS) then
               EXIT
            endif
            if (KNEL .EQ. KNELS) then
               KEYPRI = 1
               write(*, *) 'spider, sstepon: no covergence ' 
               write(*, *) 'KNEL .GE. KNELS  ', KNEL, KNELS
            endif
         enddo

         if (kpr.ge.0) then
            write(fname, '(a, a)') TRIM(path), '/knel_iters.wr'
            open(1, file=fname)
               write(1, *) KNEL
            close(1)
         endif

         diftok=dabs(platok-ztok)/(dabs(platok-ztok_n)+1.d-8)
         difpsi=dabs(psax-zpsim)/(dabs(psax-zpsim_n)+1.d-8)

         if ((diftok.lt.5.0d-3 .AnD. it_dmf.gt.0) .OR. keyctr.eq.0) then
            if (ipsmk.lt.1.or.kpr.ge.-1) then
               call eqb( alf0, alf1, alf2, bet0, bet1, bet2, 
     &             alw0, alw1, alw2, betplx, i_betp, 
     &             keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &             n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &             psi_bnd, psi0_bnd )
            endif
            EXIT
         endif

         ztok  = platok
         zpsim = psax

      enddo ! End of huge iteration loop

      ztok_n=platok
      zpsim_n=psax
      PSIBOU = UP
      PSIDEL = PSIOUT - PSIBOU

!EFable
      RM_EF = RM
      ZM_EF = ZM
      TSTEP_EF = TSTEP
      DT_EF = DT
      TIME_EF = TIME
      RXPNT_EF = 0.
      ZXPNT_EF = 0.
      rxpnt = rx0
      zxpnt = zx0

      return
      end subroutine SSTEPON_BKDW

!---------------------------------------------------------------------
      subroutine SSTEPON_BKDW000(dt, time)

      use durs_d_modul, only: n_tht, n_psi, keyctr, i_betp, igdf, nurs,
     &    tokf, i_eqdsk, epsro, betplx,
     &    rax, zax, b0, r0, psax, rxpnt, zxpnt,
     &    alf0, alf1, alf2, bet0, bet1, bet2 
      use sp_parameters, only: nilim, njlim, nclim, twopi
      use iopath, only: path
      use keys, only: kstep

      implicit none

      real*8, intent(in) :: dt, time

      integer :: i, j, ngra1, ngra2, numwr, k_step, kstop, knels, 
     &   ngav, ngav1, nursb, ipsmk, nfrpr1, nfrwr1, nout, nter, 
     &   ninfw, ninev, nequi, ncequi, ngrid, n_ctrl, nctrl, i_bsh,
     &   numlim, eq_cmd
      real*8 :: alw0, alw1, alw2, alp, alp_b, alpnew, helinp, helout, 
     &   platok, ftok, tokout, qcen,
     &   psiax, psicen, psi_bnd, psi0_bnd, psiout, psax_sta, psibou,
     &   psidel, e_psi, pspl_av, psex_av, zpsim_n, ztok_n,
     &   timev, tstep, tstepr, tstart, tstop, time_sta, time_fin,
     &   bbb, sigm, betap0, betpol, zli3, up, beold, zlold,
     &   rm, zm, rx0, zx0, rmax0, zmax0, zcen0, z0cen, raxpr, zaxpr,
     &   rmaold, zmaold, rxppr, zxppr, rxpold, zxpold
      real*8, dimension(nclim) :: psip
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, 
     &   pjkd
      real*8, dimension(njlim, njlim) :: res
      double precision :: dteqz 
      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/ ncequi
      common /comhel/ helinp, helout
      common /comsta/ platok, eqdfn, i_bsh
      common /com234/ betpol, tokout, psiout
      common /comst0/ res, volk, volkp1
      common /comst1/ pjk, pjkp1, pjkp, pjkd
      common /timingcmdipsmk/ ipsmk, eq_cmd, dteqz

      kstep = 0
      timev = time
      tstep = dt
      write(*, *) TRIM(path)

      if (ipsmk.lt.1) dteqz=dt !Efable

      BBB  = 1.d0 / twopi
      SIGM   = 1.d0
      TSTEPR = TSTEP
      TSTART = 0.00d0
      TSTOP  = 20.0d0
      KSTOP  = 1000
      KNELS  = 500
      NFRPR1 = 5
      NFRWR1 = 0
      NOUT  = 17
      NTER  = 6
      NINFW = 7
      NINEV = 10
      ngra1 = 14
      ngra2 = 15

      numwr  = 0
      HELOUT = HELINP

      platok = tokf
      kstep  = 0
      k_step = 0
      ngav=keyctr

      write(fname, '(a, a)') TRIM(path), '/currents.wr'
      open(1, file=fname, form='formatted')
         write(1, *) NEQUI, NCEQUI
         write(1, *)(pjk(j), j=1, NCEQUI)
      close(1)

      write(fname, '(a, a)') TRIM(path), '/res_mat.wr'
      open(1, file=fname, form='formatted')
         write(1, *) ((res(i, j), j=1, ncequi), i=1, ncequi)
      close(1)

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         ngav, k_step, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      ngrid=1
      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf

      call eq_0(pjk, psip, ncequi, kstep, ngrid, 
     &          alf0, alf1, alf2, bet0, bet1, bet2, 
     &          betplx, ngav1, 
     &          ftok, tokout, psiax, nursb, 
     &          psi_bnd, alp_b, rax, zax, n_ctrl, b0, r0)
      call rdexf(ncequi)

      call eq(pjk, psip, ncequi, kstep, ngrid, 
     &        alf0, alf1, alf2, bet0, bet1, bet2, 
     &        betpol, betplx, ngav1, 
     &        tokout, psiout, 
     &        nursb, psi_bnd, alp_b, rax, zax)

      call eq_par(z0cen, alp, alpnew, qcen, nctrl, numlim, up, 
     &            rm, zm, rx0, zx0)

      RMAX0  = RM
      ZMAX0  = ZM
      ZCEN0  = Z0CEN
      BEOLD  = BETPOL
      ZLOLD  = ZLI3
      RAXPR  = RM
      ZAXPR  = ZM
      RMAOLD = RM
      ZMAOLD = ZM
      RXPPR  = RX0
      ZXPPR  = ZX0
      RXPOLD = RX0
      ZXPOLD = ZX0

!--------------------------------------------------------
! Definition of input parameters from "basic" equilibrium
      BETPLX = BETPOL
      BETAP0 = BETPOL
      FTOK   = TOKOUT
      PSIAX  = PSIOUT
      HELINP = HELOUT

!----------------------------------------------------
! For PSI_axis value time scenario (for NGAV1=2 or 3)

      time_sta = TSTART
      time_fin = TSTOP
      psax_sta = PSIOUT

!***********************************************************************
      PSIBOU = UP
      PSIDEL = PSIOUT - PSIBOU

! --- If we need to compute the "basic" equilibrium only -> KSTOP=0

      i_bsh  = -1
      i_eqdsk  = 0
      kstep  =  0
      keyctr =  0
      e_psi  =  0.d0

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2, 
     &         betplx, i_betp, 
     &         keyctr, kstep, platok, rax, zax, b0, r0, psax, igdf, 
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, 
     &         psi_bnd, psi0_bnd)

      call get_par(psi_bnd)
      call bongri
      call psib_pla(pspl_av)
      call psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av

      ztok_n  = platok
      zpsim_n = psax
      rxpnt = rx0
      zxpnt = zx0

      return
      end subroutine SSTEPON_BKDW000
