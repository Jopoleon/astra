      subroutine B_STEPON( KLUCH, k_auto, nstep, dt, key_dmf, dpsdt)

      use durs_d_modul, only: n_tht, n_psi, igdf, nurs, keyctr,
     &    i_eqdsk, i_betp, epsro, betplx, tokf, rax, zax, psax, b0, r0,
     &    alf0, alf1, alf2, bet0, bet1, bet2
      use sp_parameters, only: nstep_p
      use iopath, only: path
      use keys, only: kstep,  kastr
      use curpl, only: cur_pl

      implicit none

      integer, intent(in) :: KLUCH, k_auto, nstep, key_dmf
      real*8, intent(in) :: dt, dpsdt

      integer :: i_bsh
      real*8 :: alw0, alw1, alw2, psi_bnd, psi0_bnd, plat_ok,
     &   pspl_av, psex_bnd, flfi_m
      real*8, dimension(300) :: psi_ext_bon
      character(len=80) :: fname

      save plat_ok,psex_bnd,eqdfn,i_bsh

      common/psi_test/ psi_ext_bon

      character(len=40) :: eqdfn

      kstep = nstep

      if (KLUCH .NE. 0) goto 1111

      if (k_auto.eq.0 .or. kastr .eq. 1) goto 2005

      write(fname,'(a,a)') TRIM(path), '/durs.dat'
      open(1,file=fname)
         read(1,*) n_tht
         read(1,*) n_psi
         read(1,*) igdf
         read(1,*) epsro
         read(1,*) nurs
         read(1,*) keyctr
         read(1,*) i_eqdsk
         read(1,'(A40)') eqdfn 
         read(1,*) i_betp
         read(1,*) betplx
         read(1,*) tokf
         read(1,*) psax
         read(1,*) b0,r0
         read(1,*) alf0
         read(1,*) alf1 !! P'=alf0*(1-(1-psi)**alf1)**alf2
         read(1,*) alf2
         read(1,*) bet0
         read(1,*) bet1 !! FF'=bet0*(1-(1-psi)**bet1)**bet2
         read(1,*) bet2
         read(1,*) alw0
         read(1,*) alw1 !! (ro*w**2)'=alw0*(1-(1-psi)**alw1)**alw2
         read(1,*) alw2
         read(1,*) rax
         read(1,*) zax
      close(1)

      if (i_eqdsk.eq.1) then
         call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
         nurs   = -3999
         i_betp = 0
      endif

      cur_pl = tokf

 2005 continue

      if (kastr.eq.0) then 
         write(fname,'(a,a)') TRIM(path), '/inpol.dat'
         open(1,file=fname,form='formatted')
            read(1,*) i_bsh
         close(1)    
      else
         i_bsh=1
      endif

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, tokf, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,psi0_bnd)

      plat_ok = tokf
      call bongri
      call psib_pla(pspl_av)
      psex_bnd=-pspl_av

      return

1111  continue

      keyctr=key_dmf
      plat_ok = cur_pl
      call bongri
      call psib_pla(pspl_av)
      call get_flfi(flfi_m)
 
      psi0_bnd=psex_bnd+pspl_av
      psex_bnd=psex_bnd+dpsdt*dt
      psi_bnd=psex_bnd

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &          betplx, i_betp,
     &          keyctr, nstep, plat_ok, rax,zax, b0,r0, psax, igdf,
     &          n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &          psi_bnd,psi0_bnd)

      return
      end subroutine B_STEPON

!---------------------------------------------------------------------
      subroutine F_STEPON(KLUCH, k_auto, nstep, dt, time, voltpf, 
     &                    key_dmf)

      use durs_d_modul, only: n_tht, n_psi, igdf, nurs, keyctr,
     &    i_eqdsk, i_betp, epsro, betplx, tokf, rax, zax, psax, b0, r0,
     &    alf0, alf1, alf2, bet0, bet1, bet2,
     &    rxpnt, zxpnt, psbo
      use sp_parameters, only: nloopp, nprobp, njlim, nplim, npfc0, 
     &    nilim, nclim, pi, twopi, nstep_p
      use iopath, only: path
      use comevl, only: pfceqw, nequi
      use keys, only: kstep, kpr
      use e_nels, only: enels

      implicit none

      integer, intent(in) :: kluch, k_auto, nstep, key_dmf
      real*8, intent(in) :: dt, time
      real*8, intent(in), dimension(*) :: voltpf

      integer :: i, j, k, l, nout, nter, ninfw, ninev, nursb, istep,
     &   ngra1, ngra2, i_bsh, nflag, nbp, nvv, nreg, nles,
     &   kstop, knels, kstepr, keypri, knel,
     &   ipsmk, eq_cmd, nfrpr1, nfrwr1,
     &   numwr, cmnd_dioh2s, cmnd_dioh2u, nprob, nloop, nc, ncequi,
     &   ncpfc, nfw
      integer :: NECON(NILIM), NTYPE(NCLIM)

      real*8 :: timev, sigm, tstep, bbb, tstepr, tstart, tstop, 
     &   alw0, alw1, alw2, psi_bnd, platok, psi0_bnd, psdel, sgmcur,
     &   erro, betpol, pspl_av, psex_av, dpsdt,
     &   ereve, psi_eav, psi_eav_n
      real*8, dimension(nclim) :: RC, ZC, PC, PSIP, VC, HC,
     &   RC1, RC2, RC3, RC4, ZC1, ZC2, ZC3, ZC4
      real*8, dimension(NILIM) :: wecon
      real*8, dimension(njlim) :: VOLK, VOLKP1, PJK, PJKP1, PJKP, PJKD,
     &   PSK, PSKP1, PSKP, PSKM1
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      real*8, dimension(njlim, njlim) :: RES
      real*8 :: errarr(10)

      double precision :: dteqz

      character(len=80) :: fname
      character(len=40) :: eqdfn

      common /commanddioh/ cmnd_dioh2s, cmnd_dioh2u
      common /comeqg/ ncequi
      common /comsta/ platok, eqdfn, i_bsh
      common /comst0/ RES, VOLK, VOLKP1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comst2/ PSK, PSKP1, PSKP, PSKM1
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob
      common /timingcmdipsmk/ ipsmk, eq_cmd, dteqz

      save numwr
      save psi_eav_n

      kstep = nstep
      timev = time
      tstep = dt
      if (ipsmk.lt.1) dteqz=dt !Efable

      BBB  = 1.d0 / twopi

      NOUT  = 17
      NTER  = 6
      NINFW = 7
      NINEV = 10
      ngra1 = 14
      ngra2 = 15
      SIGM   = 1.0d0  !0.d0
      TSTEPR = TSTEP
      TSTART = 0.00d0
      TSTOP  = 20.0d0
      KSTOP  = 1000
      KNELS  = 500
      NFRPR1 = 5
      NFRWR1 = 0
      KSTEPR = KSTEP
      KEYPRI = 1

      if (KLUCH .NE. 0) goto 1111

      numwr  = 0

!-----------------------------------------------------------------
! --- INPUT OF "BASIC" EQULIBRIUM PARAMETERS FROM FILE "durs.dat"

      if (k_auto.ne.0) then
         write(fname,'(a,a)') TRIM(path), '/durs.dat'
         open(1,file=fname,form='formatted')
            read(1,*) n_tht
            read(1,*) n_psi
            read(1,*) igdf
            read(1,*) epsro
            read(1,*) nurs
            read(1,*) keyctr
            read(1,*)    i_eqdsk
            read(1, '(A40)') eqdfn
            read(1,*) i_betp
            read(1,*) betplx
            read(1,*) tokf
            read(1,*) psax
            read(1,*) b0,r0
            read(1,*) alf0
            read(1,*) alf1 ! P'=alf0*(1-(1-psi)**alf1)**alf2
            read(1,*) alf2
            read(1,*) bet0
            read(1,*) bet1 ! FF'=bet0*(1-(1-psi)**bet1)**bet2
            read(1,*) bet2
            read(1,*) alw0
            read(1,*) alw1 ! (ro*w**2)'=alw0*(1-(1-psi)**alw1)**alw2
            read(1,*) alw2
            read(1,*) rax
            read(1,*) zax
         close(1)

!-----------------------------------------
! INPUT OF POSITIONS OF "PF_PROBE" POINTS:
         call PROPNT(NOUT, NTER, NINFW, NGRA1,
     &               NPROb, RPROb, ZPROb, FIPROb)

!----------------------------------------
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:

         call LOOPNT(NOUT, NTER, NINFW, NGRA1,
     &               NLOOp, RLOOp, ZLOOp)

!-----------------------------------------------------
! INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS

         call CONDUC(NC, NCEQUI, NCPFC, NFW, NBP, NVV,
     &               RC, ZC, PC, VC, HC, NTYPE,
     &               RC1, ZC1, RC2, ZC2, RC3, ZC3, RC4, ZC4,
     &               RES, VOLK, VOLKP1,
     &               NECON, WECON,
     &               NOUT, NTER, NINFW, ngra1)

!-----------------------------------------------------
! DEFINITION INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  "PPIND" from COMMON /PPIDPS/

         call L_MATR(NOUT, NTER, NC, NCPFC,
     &               NTYPE, RC, ZC, VC, HC,
     &               NECON,WECON)

!---------------------------------------------------
! Initial condition (currents) for circuit equations

         do L=1,NCEQUI
	    if (L.LE.NEQUI) then
               PJK(L)  = PFCEQW(L)
               PJKP(L) = PJK(L)
            else
               PJK(L)  = PC(NCPFC+L-NEQUI)
               PJKP(L) = PJK(L)
            endif
         enddo

!Efable this block is to keep dioh2u and dioh2s to zero for some time
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

        write(fname,'(a,a)') TRIM(path), '/currents.wr'
        open(1,file=fname,form='formatted')
           write(1,*) NEQUI,NCEQUI
           write(1,*)(pjk(j),j=1,NCEQUI)
        close(1)

        write(fname,'(a,a)') TRIM(path), '/res_mat.wr'
        open(1,file=fname,form='formatted')
           write(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
        close(1)

        call wrcoil(nc, ncpfc, rc, zc, pc, necon, wecon)

! TOROIDAL CURRENTS OF PASSIVE CONDUCTOR STRUCTURES 
      endif

      call rd_ppind
      call rd_prob(NPROb, RPROb, ZPROb,  FIPROb)
      call rd_loop(NLOOp, RLOOp, ZLOOp)

      write(fname,'(a,a)') TRIM(path), '/currents.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi,ncequi
         read(1,*)(pjk(j),j=1,ncequi)
      close(1)

      write(fname,'(a,a)') TRIM(path), '/res_mat.wr'
      open(1,file=fname,form='formatted')
         read(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

      call rdcoil(nc, ncpfc, rc, zc, pc, necon, wecon)

      do L=1,NCEQUI
         PJKP(L) = PJK(L)
         PJKP1(L) = PJK(L)
      enddo

! voltage for zero time-level

      do i=1,NEQUI
         VOLK(i) = voltpf(i)*BBB
         VOLKp1(i) = VOLK(i)
      enddo
      do i=NEQUI+1,NCEQUI
         VOLK(i) = 0.0d0
         VOLKp1(i) = VOLK(i)
      enddo

      kstep = 0
      nursb=nurs
      platok=tokf
      i_eqdsk=0
      i_bsh=-1

      call eqb(alf0, alf1,alf2, bet0, bet1, bet2, alw0, alw1, alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax, zax, b0, r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh,
     &         psi_bnd, psi0_bnd)

      call eqa_in(alf0, alf1, alf2, bet0, bet1, bet2, platok,  
     &            pjk, ncequi, b0, r0)

      if (kpr.eq.1) then
         write(*,*) 'START  OF BASIC FREE BOUNDARY EQUILIBRIUM'
      endif

      call eqa(keyctr, igdf, platok, psax, i_betp, betplx, 
     &         rax, zax, rxpnt, zxpnt, psbo, psdel,
     &         ncequi, psip, betpol, nflag, errarr)
      call f_wrd
      call renet
      call f_bndmat
      call f_wrd

      call eqa(keyctr, igdf, platok, psax, i_betp, betplx, 
     &         rax, zax, rxpnt, zxpnt, psbo, psdel,
     &         ncequi, psip, betpol, nflag, errarr)

      do k=1,ncequi
         psk(k)=psip(k)
         pskp1(k)=psk(k)
      enddo

      istep=nstep
      call bongri
      call f_psib_pla(pspl_av)
      call f_psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      call put_psib0(psi0_bnd)
      call wr_step(numwr, time, istep)

      return

1111  continue

      keyctr=key_dmf

      if (keyctr .eq. 0) platok = tokf
      call bongri
      call f_psib_pla(pspl_av)
      call f_psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      call put_psib0(psi0_bnd)
      psi_eav=psex_av
      dpsdt=(psi_eav-psi_eav_n)/dt
      psi_eav_n=psi_eav

!-----------------------------------------------------
! Computing of "VOLKP1"  voltage  for CIRCUIT EQUATION 
! for new time level KSTEP

      do i=1,NEQUI
         VOLKP1(i) = voltpf(i)*BBB
      enddo
      do i=NEQUI+1,NCEQUI
         VOLKP1(i) = 0.0d0
      enddo

!--------------------------------------------------------
! preparations for currrents - equilibrium iteration loop

      if (KSTEP.EQ.1) then
         do L=1,NCEQUI
            PSKP1(L) = PSK(L)
         enddo
         NREG = 0
         NLES = 0

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
         call EVSLV(NLES, NREG, TSTEP, TSTEP, SIGM,
     &              NCEQUI, VOLK, VOLKP1, RES,
     &              PSK, PSKP1, PJK, PJKP1, PJKP,
     &              NOUT, NTER, KEYPRI, EREVE, dteqz)

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

      else !if (KSTEP.NE.1)

!----------------------------------------------------------------
!  The initial approximation for the case of closed-loop evolution

         do L=1,NCEQUI
            PJKP(L) = PJK(L)
         enddo
         NREG = 0
         NLES = 1

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
         call EVSLV(NLES, NREG, TSTEPR, TSTEP, SIGM,
     &              NCEQUI, VOLK, VOLKP1, RES,
     &              PSKM1, PSK, PJK, PJKP1, PJKP,
     &              NOUT, NTER, KEYPRI, EREVE, dteqz)

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
      endif   ! if KSTEP.NE.1 

	 if (eq_cmd.eq.1.and.ipsmk.ge.1) then
         call eqa_ax(dt, time,
     &             keyctr, igdf, nstep, platok, psax, i_betp, betplx, 
     &             rax, zax, rxpnt, zxpnt, psbo, psdel,
     &             pjkp1, ncequi, psip, betpol, nflag, errarr)

         PSKM1 = PSK
         PSK   = PSKP1
         do k=1,ncequi
            pskp1(k)=psip(k)
         enddo
      endif

! currrents - equilibrium iteration loop begining
 
      do KNEL=1, KNELS

         do L=1,NCEQUI
            PSKP(L) = PSKP1(L)
            PJKP(L) = PJKP1(L)
         enddo
         NREG = 0
         NLES = 2

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

         call EVSLV(NLES, NREG, TSTEP, TSTEP, SIGM,
     &                 NCEQUI, VOLK,  VOLKP1, RES,
     &                 PSK, PSKP1, PJK, PJKP1, PJKP,
     &                 NOUT, NTER, KEYPRI, EREVE, dteqz)

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

         SGMCUR = 0.5D0
         do L=1,NCEQUI
            PJKP1(L) = SGMCUR*PJKP1(L) + (1.0D0 - SGMCUR)*PJKP(L)
         enddo

!added kpr = -2 this for skip EFable fsim
         if (kpr.ge.-1) then
            call eqa_ax(dt, time,
     &               keyctr, igdf, nstep, platok, psax, i_betp, betplx, 
     &               rax, zax, rxpnt, zxpnt, psbo, psdel,
     &               pjkp1, ncequi, psip, betpol, nflag, errarr)
            do k=1,ncequi
               pskp1(k)=psip(k)
            enddo
         endif

         erro=errarr(1)
	    if (kpr.lt.-1) erro=ereve !for skip EFable fsim
         if (kpr.eq.1) then
            write(*,*) ' '
            write(*,*) ' '
            write(*,*) 'stepon:erro',erro
            write(*,*) ' '
            write(*,*) ' '
         endif

         if (erro .LT. ENELS) then
            if (kpr.ge.0) then
               write(fname,'(a,a)') TRIM(path), '/knel_iters.wr'
               open(1,file=fname)
                  write(1,*) KNEL
               close(1)
            endif
            EXIT
         endif

         if (KNEL .EQ. KNELS) then
            if (kpr.eq.1) then
               write(*,*) ' '
               write(*,*) 'stepon:limit of iterations is exceded'
               write(*,*) 'knel=',knel
               write(*,*) ' '
            endif
         endif

      enddo ! currrents - equilibrium iteration loop finish

      call f_wrd

      return
      end  subroutine F_STEPON

!----------------------------------------------------------------
      subroutine F_STEPON_BKDW(KLUCH, nstep, dt, time, voltpf, key_dmf)

      use durs_d_modul, only: n_tht, n_psi, igdf, nurs, keyctr,
     &    i_eqdsk, i_betp, epsro, betplx, tokf, rax, zax, psax, b0, r0,
     &    alf0, alf1, alf2, bet0, bet1, bet2,
     &    rxpnt, zxpnt, psbo
      use sp_parameters, only: nloopp, nprobp, njlim, nplim, npfc0, 
     &    nilim, nclim, pi, twopi, nstep_p
      use iopath, only: path
      use comevl, only: nequi
      use keys, only: kstep, kpr
      use e_nels, only: enels

      implicit none

      integer, intent(in) :: kluch, nstep, key_dmf
      real*8, intent(in) :: dt, time
      real*8, intent(in), dimension(*) :: voltpf

      integer :: i, j, k, l, nout, nter, ninfw, ninev, nursb, istep,
     &   ngra1, ngra2, i_bsh, nflag, nreg, nles,
     &   kstop, knels, kstepr, keypri, knel,
     &   ipsmk, eq_cmd, nfrpr1, nfrwr1,
     &   numwr, cmnd_dioh2s, cmnd_dioh2u, nprob, nloop, nc, ncequi,
     &   ncpfc, max_iteri, max_max_iteri
      integer :: NECON(NILIM), NTYPE(NCLIM)

      real*8 :: timev, sigm, tstep, bbb, tstepr, tstart, tstop, 
     &   alw0, alw1, alw2, psi_bnd, platok, psi0_bnd, psdel, sgmcur,
     &   erro, betpol, betful, pspl_av, psex_av, dpsdt,
     &   ereve, psi_eav, psi_eav_n
      real*8, dimension(nclim) :: RC, ZC, PC, PSIP
      real*8, dimension(NILIM) :: wecon
      real*8, dimension(njlim) :: VOLK, VOLKP1, PJK, PJKP1, PJKP, PJKD,
     &   PSK, PSKP1, PSKP, PSKM1
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      real*8, dimension(njlim, njlim) :: RES
      real*8 :: errarr(10)

      double precision :: dteqz, dt_nowww

      character(len=80) :: fname
      character(len=40) :: eqdfn

      common /commanddioh/ cmnd_dioh2s, cmnd_dioh2u
      common /comeqg/ ncequi
      common /comsta/ platok, eqdfn, i_bsh
      common /comst0/ RES, VOLK, VOLKP1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comst2/ PSK, PSKP1, PSKP, PSKM1
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob
      common /timingcmdipsmk/ ipsmk, eq_cmd, dteqz
      common /maxiterii/ max_max_iteri

      save numwr
      save psi_eav_n, dt_nowww

      kstep = nstep
      timev = time
      tstep = dt
      if (ipsmk.lt.1) dteqz=dt

      BBB  = 1.d0 / twopi

      NOUT  = 17
      NTER  = 6
      NINFW = 7
      NINEV = 10
      ngra1 = 14
      ngra2 = 15
      SIGM   = 1.0d0  !0.d0
      TSTEPR = TSTEP
      TSTART = 0.00d0
      TSTOP  = 20.0d0
      KSTOP  = 1000
      KNELS  = 500
      NFRPR1 = 5
      NFRWR1 = 0
      KSTEPR = KSTEP
      KEYPRI = 1

      if (KLUCH .NE. 0) goto 1111

      call rd_ppind

      write(*,*) TRIM(path)
      call rdcoil(nc,ncpfc,rc,zc,pc,necon,wecon)


      write(fname,'(a,a)') TRIM(path), '/currents.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi,ncequi
         read(1,*)(pjk(j),j=1,ncequi)
      close(1)
      write(*,*) ncequi
      write(fname,'(a,a)') TRIM(path), '/res_mat.wr'
      open(1,file=fname,form='formatted')
         read(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

      dt_nowww=TSTEP

      do L=1,NCEQUI
         PJKP(L) = PJK(L)
         PJKP1(L) = PJK(L)
      enddo

      do i=1,NEQUI
         VOLK(i) = voltpf(i)*BBB
         VOLKp1(i) = VOLK(i)
      enddo
      do i=NEQUI+1,NCEQUI
         VOLK(i) = 0.0d0
         VOLKp1(i) = VOLK(i)
      enddo

      kstep = 0
      nursb=nurs
      platok=tokf
      i_eqdsk=0
      i_bsh=-1
      write(*,*) 'eqb'
      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,psi0_bnd)

      write(*,*) 'eqain'

      call eqa_in(alf0,alf1,alf2,bet0,bet1,bet2,platok,  
     &            pjk,ncequi, b0,r0)

      write(*,*) 'eqaout'

      if (kpr.eq.1) then
         write(*,*) 'START  OF BASIC FREE BOUNDARY EQUILIBRIUM'
      endif

      write(*,*) 'eqa'
      call eqa(keyctr, igdf, platok, psax, i_betp, betplx, 
     &         rax, zax, rxpnt,zxpnt, psbo, psdel,
     &         ncequi, psip, betpol, nflag, errarr)
      write(*,*) 'eqa out'

      call f_wrd
      call renet
      call f_bndmat
      call f_wrd

      call eqa(keyctr, igdf, platok, psax, i_betp, betplx, 
     &         rax, zax, rxpnt,zxpnt, psbo, psdel,
     &         ncequi, psip, betpol, nflag, errarr)

      istep=nstep
      call bongri
      call f_psib_pla(pspl_av)
      call f_psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      call put_psib0(psi0_bnd)
      call wr_step(numwr,time,istep)

      return

1111  continue

      keyctr=key_dmf
      if (keyctr .eq. 0) platok = tokf
      call bongri
      call f_psib_pla(pspl_av)
      call f_psib_ext(psex_av)

      psi0_bnd=pspl_av+psex_av
      call put_psib0(psi0_bnd)
      psi_eav=psex_av
      dpsdt=(psi_eav-psi_eav_n)/dt
      psi_eav_n=psi_eav

!  Computing "VOLKP1" voltage for CIRCUIT EQUATION 
!  for new time level KSTEP

      do i=1,NEQUI
         VOLKP1(i) = voltpf(i)*BBB
      enddo
      do i=NEQUI+1,NCEQUI
         VOLKP1(i) = 0.0d0
      enddo

! preparations for currrents - equilibrium iteration loop

      if (KSTEP.EQ.1) then
         do L=1,NCEQUI
            PSKP1(L) = PSK(L)
         enddo
         NREG = 0
         NLES = 0

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

         call EVSLV(NLES, NREG, TSTEP, TSTEP, SIGM,
     &              NCEQUI, VOLK, VOLKP1, RES,
     &              PSK, PSKP1, PJK, PJKP1, PJKP,
     &              NOUT, NTER, KEYPRI, EREVE, dteqz)

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

      endif    ! if (KSTEP.EQ.1)

      if (eq_cmd.eq.1.and.ipsmk.ge.1) then
         PSKM1 = PSK
         PSK   = PSKP1
         call eqa_ax(dt, time,
     &               keyctr, igdf, nstep,platok, psax, i_betp, betplx, 
     &               rax, zax, rxpnt,zxpnt, psbo, psdel,
     &               pjk, ncequi, psip, betpol, nflag, errarr)

         do k=1,ncequi
            pskp1(k)=psip(k)
         enddo
      endif

      if (KSTEP.NE.1) then

!----------------------------------------------------------------
!  The initial approximation for the case of closed-loop evolution

         do L=1,NCEQUI
            PJKP(L) = PJK(L)
         enddo

         if (dt_nowww.eq.TSTEP) then
            NREG = 0
            NLES = 1
         else
            NREG = 0
            NLES = 0
         endif

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

         call EVSLV(NLES, NREG, TSTEPR, TSTEP, SIGM,
     &              NCEQUI, VOLK, VOLKP1, RES,
     &              PSK, PSKP1, PJK, PJKP1, PJKP,
     &              NOUT, NTER, KEYPRI, EREVE, dteqz)

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

      endif   !  if (KSTEP.NE.1)

      max_iteri = 0 !Efable

      do KNEL=1, KNELS

         do L=1,NCEQUI
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

         call EVSLV(NLES, NREG, TSTEP, TSTEP, SIGM,
     &              NCEQUI, VOLK, VOLKP1, RES,
     &              PSK, PSKP1, PJK, PJKP1, PJKP,
     &              NOUT, NTER, KEYPRI, EREVE ,tstep)

         SGMCUR = 0.25D0

         do L=1,NCEQUI
            PJKP1(L) = SGMCUR*PJKP1(L) + (1.0D0 - SGMCUR)*PJKP(L)
         enddo

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

	 if (kpr.ge.-1) then
            if (eq_cmd.eq.1.and.ipsmk.ge.1) then
               call eqa_ax(dt, time,
     &               keyctr, igdf, nstep, platok, psax, i_betp, betplx, 
     &               rax, zax, rxpnt, zxpnt, psbo, psdel,
     &               pjkp1, ncequi, psip, betpol, nflag, errarr)
               do k=1,ncequi
                 pskp1(k)=psip(k)
               enddo
            endif
	 endif
	 if (kpr.eq.-2.and.max_iteri.le.max_max_iteri) then
            if (eq_cmd.eq.1.and.ipsmk.ge.1) then
               call eqa_ax(dt, time,
     &                 keyctr, igdf, nstep, platok, psax,i_betp, betplx, 
     &                 rax, zax, rxpnt,zxpnt, psbo, psdel,
     &                 pjkp1,ncequi,betpol,betful,     
     &                 necon,wecon,ntype , nflag, errarr)
               do k=1,ncequi
                  pskp1(k)=psip(k)
               enddo
               max_iteri=max_iteri+1
            endif
         endif

         erro=errarr(1)
	 if (kpr.lt.-1.and.max_iteri.gt.max_max_iteri) erro=ereve !for skip EFable fsim
         if (kpr.eq.1) then
            write(*,*) ' '
            write(*,*) ' '
            write(*,*) 'stepon:erro',erro
            write(*,*) ' '
            write(*,*) ' '
         endif

         if (erro .LT. ENELS) then
            if (kpr.ge.0) then
               write(fname,'(a,a)') TRIM(path), '/knel_iters.wr'
               open(1,file=fname)
                  write(1,*) KNEL
               close(1)
            endif
            EXIT
         endif

         if (KNEL .EQ. KNELS) then
            if (kpr.eq.1) then
               write(*,*) ' '
               write(*,*) 'stepon:limit of iterations is exceded'
               write(*,*) 'knel=',knel
               write(*,*) ' '
            endif
         endif

      enddo ! currrents - equilibrium iteration loop finish

      call f_wrd

      return
      end subroutine F_STEPON_BKDW
