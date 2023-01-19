      subroutine cf_init(k_auto, nstep, dt, time,
     &                   voltpf, d_pf_mat, d_tcam_mat)

      use durs_d_modul
      use sp_parameters
      use iopath, only: kname, path
      use parcur, only: curref
      use comevl
      use keys

      implicit none

      integer, intent(in) :: k_auto
      integer, intent(out) :: nstep
      real*8, intent(in), dimension(*) :: voltpf, d_pf_mat, d_tcam_mat

      integer :: i, j, l, ncequi, i_bsh, ik, iq, ncpfc, ngrid, nvv, 
     &   nbp, nfw, nc, nloop, n_ctrl, ngav1, ngra1, nprob, nursb,
     &   nout, nter, isymm, ninfw, ninf
      integer, dimension(nclim) :: ntype
      integer, dimension(nilim) :: necon
      real*8 :: alp_b, psi_bnd, psi0_bnd, platok, dt, time, ereve0,
     &   erps, e_psi, pscout, tokout, zli3, betpol, ftok, psicen
      real*8, dimension(nclim) :: pc, psip, vc, hc, ccurx, ccury,
     &   rc, rc1, rc2, rc3, rc4, zc, zc1, zc2, zc3, zc4, 
     &   alw0, alw1, alw2 
      real*8, dimension(nilim) :: wecon
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, pjkd,
     &   psk, pskp1, pskp, pskm1
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      real*8, dimension(njlim, njlim) :: res

      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/ ncequi
      common /comsta/ platok, eqdfn, i_bsh
      common /comst0/ res, volk, volkp1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comst2/ PSK,  PSKP1, PSKP, PSKM1
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob

      kstep = nstep
      write(*,*) kstep, k_auto, kastr

      ninf=1       
      kxwx=1  !#prescrib. x-points;
      ksnf=0  ! if ksnf=1,then kxwx must be =1

      if(k_auto.eq.0) goto 2005

      if(kastr.ne.1) then
         write(fname,'(a,a)') path(1:kname),'durs.dat'
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
         if(i_eqdsk.eq.1) then
            call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
            nurs   = -3999
            i_betp = 0
         endif
      endif

!-----------------------------------------
! INPUT OF POSITIONS OF "PF_PROBE" POINTS:
      call PROPNT(NOUT, NTER, NINFW, NGRA1,
     &            NPROb, RPROb, ZPROb, FIPROb)

!----------------------------------------
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:
      call LOOPNT(NOUT, NTER, NINFW, NGRA1,
     &            NLOOp, RLOOp, ZLOOp)

!-----------------------------------------------------
! INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS

      call CONDUC(NC, NCEQUI, NCPFC, NFW, NBP, NVV,
     &            RC,ZC, PC, VC,HC, NTYPE,
     &            RC1,ZC1,  RC2,ZC2,  RC3,ZC3,  RC4,ZC4,
     &            RES, VOLK, VOLKP1,
     &            NECON, WECON,
     &            NOUT, NTER, NINFW, ngra1 )

!-----------------------------------------------------
! DEFINITION INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  "PPIND" from COMMON /PPIDPS/

      call L_MATR(NOUT,  NTER, NC, NCPFC,
     &            NTYPE, RC,   ZC, VC, HC,
     &            NECON,WECON )

      do L=1,NCEQUI
         if ( L.LE.NEQUI ) then
            PJK(L)  = PFCEQW(L)
         else
            PJK(L)  = PC(NCPFC+L-NEQUI)
         endif
      enddo

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'pfcurr.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI
         write(1,*)(PFCEQW(j),j=1,NEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'res_mat.wr'
         open(1,file=fname,form='formatted')
         write(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

!  TOROIDAL CURRENTS OF PASSIVE CONDUCTOR STRUCTURES 

      call wrcoil(nc, ncpfc, rc, zc, pc, necon, wecon)
      ngrid=1

      call auto(rc, zc, pc, nc, nstep, ngrid, 
     &          rc1, zc1,  rc2, zc2, 
     &          rc3, zc3,  rc4, zc4, 
     &          ntype, necon, wecon)

 2005 continue

      if(kastr.eq.0) then 
         write(fname,'(a,a)') path(1:kname),'inpol.dat'
         open(1,file=fname,form='formatted')
            read(1,*) i_bsh
         close(1)    
      else
         i_bsh=1
      endif

      call noauto
      call rd_ppind
      call rd_prob( NPROb, RPROb, ZPROb,  FIPROb )
      call rd_loop( NLOOp, RLOOp, ZLOOp )

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi,ncequi
         read(1,*)(pjk(j),j=1,ncequi)
      close(1)

      write(fname,'(a,a)') path(1:kname),'pfcurr.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi
         read(1,*)(PFCEQW(j),j=1,nequi)
      close(1)

      write(fname,'(a,a)') path(1:kname),'res_mat.wr'
      open(1,file=fname,form='formatted')
         read(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

      call rdcoil(nc,ncpfc,rc,zc,pc,necon,wecon)

      platok = tokf
      nstep  = 0

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax, zax, b0, r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh,
     &         psi_bnd, psi0_bnd)

      call prefit(rc, zc, ncpfc, NECON, WECON, rax, zax, alp_b, psi_bnd)

      do iq=1,NEQUI
         curref(iq)=PFCEQW(iq)
      enddo

      call curfit_(rc,zc,ncpfc,NECON,WECON,psi_bnd)

      do ik=1,nequi
         ccurx(ik)=PFCEQW(ik)
         pjk(ik)=PFCEQW(ik)
      enddo

      isymm=0 ! symmetric case

      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf
      nstep=0
      ngrid=1

      call eq_0(pjk, psk, ncequi, nstep, ngrid,
     &          alf0, alf1, alf2, bet0, bet1, bet2,
     &          betpol, betplx, zli3,
     &          ngav1,
     &          ftok, tokout, psicen, pscout, 
     &          nursb, psi_bnd, alp_b, rax, zax, n_ctrl, b0, r0 )
      call rdexf(ncequi)
      call eq_ax(pjk, psk, ncequi, nstep, ngrid,
     &           alf0, alf1, alf2, bet0, bet1, bet2,
     &           betpol,  betplx, zli3,
     &           ngav1,
     &           ftok, tokout, psicen, pscout,
     &           EREVE0, ERPS,
     &           psi_bnd, alp_b, rax, zax, isymm)
      call wrd
      call prefit(rc,zc,ncpfc,NECON,WECON,rax,zax,alp_b,psi_bnd)

      ngrid=1
      EREVE0=1.d-7


      do

         nstep=nstep+1
         if(kxwx.eq.1)then
            call precal(rc,zc,ncpfc,NECON, WECON )
         elseif(kxwx.eq.0)then
            call precal_wx(rc,zc,ncpfc,NECON, WECON )
         elseif(kxwx.eq.2)then
            call precal(rc,zc,ncpfc,NECON, WECON )
         endif

         if(n_ctrl.eq.1 .OR. kastr.eq.1) then  !limiter point
            call curfit_L(rc,zc,ncpfc,NECON,WECON,psi_bnd)
         else
            call curfit(rc,zc,ncpfc,NECON,WECON,psi_bnd)
         endif

         do ik=1,nequi
            ccurx(ik)=PFCEQW(ik)
            pjk(ik)=PFCEQW(ik)
         enddo

         call eq_ax(pjk, psk, ncequi, nstep, ngrid,
     &              alf0, alf1, alf2, bet0, bet1, bet2,
     &              betpol, betplx, zli3,
     &              ngav1,
     &              ftok, tokout, psicen, pscout,
     &              EREVE0, ERPS,
     &              psi_bnd, alp_b, rax, zax, isymm)

         do  L=1,NCEQUI
            PJKp1(L) = PJK(L)
            PSKP1(L) = PSK(L)
            PSKm1(L) = PSK(L)
         enddo

         call wrd

         if(ERPS.le.EREVE0) EXIT
      enddo

      call wrd

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'psiplcoils.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
         write(1,*) psi_bnd,tokout,ftok
      close(1)
       
      write(fname,'(a,a)') path(1:kname),'pfc_curr.wr'
      open(1,file=fname,form='formatted')
         write(1,*) 'coil currents [mA]'
         write(1,'(a6,e13.5)') 'PF1',pjk(1) !*1.d3
         write(1,'(a6,e13.5)') 'PF2',pjk(2) !*1.d3
         write(1,'(a6,e13.5)') 'F3,4',pjk(3) !*1.d3
         write(1,'(a6,e13.5)') 'F5,6',pjk(4) !*1.d3
         write(1,'(a6,e13.5)') 'CS',pjk(5) !*1.d3
      close(1)

      write(fname,'(a,a)') path(1:kname),'tcurrs.wr'
      open(1,file=fname,form='formatted')
         do ik=1,nequi
	    do j=1,NPFC
	       if(ik .eq. NEPFC(j)) then
                  write(1,'(1E25.11)') PFCEQW(ik)
                  EXIT
	       endif
            enddo
         enddo
      close(1)

      i_bsh=-1
      nstep  = 0
      keyctr =0
      i_betp =0
      i_eqdsk =0

      call eqb(alf0, alf1, alf2, bet0, bet1, bet2, alw0, alw1, alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax, zax, b0, r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh,
     &         psi_bnd, e_psi)

      return 
      end subroutine cf_init

!-----------------------------------------------------------
      subroutine cf_init_full(k_auto, nstep, dt, time,
     &                        voltpf, d_pf_mat, d_tcam_mat)

      use durs_d_modul
      use sp_parameters
      use iopath, only: kname, path
      use parcur, only: curref
      use comevl
      use keys

      implicit none

      integer, intent(in) :: k_auto
      integer, intent(out) :: nstep
      real*8, intent(in), dimension(*) :: voltpf, d_pf_mat, d_tcam_mat

      integer :: i, j, l, ncequi, i_bsh, ik, iq, ncpfc, ngrid, nvv, 
     &   nbp, nfw, nc, nloop, n_ctrl, ngav1, ngra1, nprob, nursb,
     &   nout, nter, isymm, ninfw, ninf
      integer, dimension(nclim) :: ntype
      integer, dimension(nilim) :: necon
      real*8 :: alp_b, psi_bnd, psi0_bnd, platok, dt, time, ereve0,
     &   erps, e_psi, pscout, tokout, zli3, betpol, ftok, psicen
      real*8, dimension(nclim) :: pc, psip, vc, hc, ccurx, ccury,
     &   rc, rc1, rc2, rc3, rc4, zc, zc1, zc2, zc3, zc4, 
     &   alw0, alw1, alw2 
      real*8, dimension(nilim) :: wecon
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, pjkd,
     &   psk, pskp1, pskp, pskm1
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      real*8, dimension(njlim, njlim) :: res

      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/ ncequi
      common /comsta/ platok, eqdfn, i_bsh
      common /comst0/ res, volk, volkp1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comst2/ PSK,  PSKP1, PSKP, PSKM1
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob

      kstep = nstep
      ninf=1     
      kxwx=1  ! #prescribed x-points;
      ksnf=0  ! if ksnf=1,then kxwx must be =1

      if(k_auto.eq.0) goto 2005 

      if(kastr.ne.1) then
         write(fname,'(a,a)') path(1:kname),'durs.dat'
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
         if(i_eqdsk.eq.1) then
            call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
            nurs   = -3999
            i_betp = 0
         endif
      endif

!-----------------------------------------
! INPUT OF POSITIONS OF "PF_PROBE" POINTS:
      call PROPNT(NOUT, NTER, NINFW, NGRA1,
     &            NPROb, RPROb, ZPROb, FIPROb)

!----------------------------------------
! INPUT OF POSITIONS OF "FL_LOOP" POINTS:
      call LOOPNT(NOUT, NTER, NINFW, NGRA1,
     &            NLOOp, RLOOp, ZLOOp )

!-----------------------------------------------------
! INPUT PARAMETERS OF PFC SYSTEM AND PASSIV CONDUCTORS
!
      call CONDUC(NC, NCEQUI, NCPFC, NFW, NBP, NVV,
     &            RC, ZC, PC, VC,HC, NTYPE,
     &            RC1, ZC1, RC2, ZC2, RC3, ZC3, RC4, ZC4,
     &            RES, VOLK, VOLKP1,
     &            NECON, WECON,
     &            NOUT, NTER, NINFW, ngra1 )

!-----------------------------------------------------
! DEFINITION INDUCT. AND  SELFINDUCT. MATRIX
! FOR "EDDY" CONDUCTORS:  "PPIND" from COMMON /PPIDPS/
      call L_MATR(NOUT, NTER, NC, NCPFC,
     &            NTYPE, RC, ZC, VC, HC,
     &            NECON, WECON )

      write(*,*) '**'
      do L=1,NCEQUI
         if ( L.LE.NEQUI ) then
            PJK(L)  = PFCEQW(L)
         else
            PJK(L)  = PC(NCPFC+L-NEQUI)
         endif
      enddo

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'pfcurr.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI
         write(1,*)(PFCEQW(j),j=1,NEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'res_mat.wr'
      open(1,file=fname,form='formatted')
         write(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

!--------------------------------------------------
! TOROIDAL CURRENTS OF PASSIVE CONDUCTOR STRUCTURES 

      call wrcoil(nc,ncpfc,rc,zc,pc,necon,wecon)

      ngrid=1

      call auto(rc,zc,pc,nc,nstep,ngrid,
     &          rc1,zc1, rc2,zc2,
     &          rc3,zc3, rc4,zc4,
     &          ntype, necon, wecon )

 2005   continue

      if(kastr.eq.0) then 
         write(fname,'(a,a)') path(1:kname),'inpol.dat'
         open(1,file=fname,form='formatted')
            read(1,*) i_bsh
         close(1)
      else
         i_bsh=1
      endif

      call noauto
      call rd_ppind
      call rd_prob(NPROb, RPROb, ZPROb, FIPROb)
      call rd_loop(NLOOp, RLOOp, ZLOOp)

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi,ncequi
         read(1,*)(pjk(j),j=1,ncequi)
      close(1)

      write(fname,'(a,a)') path(1:kname),'pfcurr.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nequi
         read(1,*)(PFCEQW(j),j=1,nequi)
      close(1)

      write(fname,'(a,a)') path(1:kname),'res_mat.wr'
      open(1,file=fname,form='formatted')
         read(1,*) ((res(i,j),j=1,ncequi),i=1,ncequi)
      close(1)

      call rdcoil(nc,ncpfc,rc,zc,pc,necon,wecon)

      platok = tokf
      nstep  = 0

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,psi0_bnd)

      call prefit(rc,zc,ncpfc,NECON,WECON,rax,zax,alp_b,psi_bnd)

      do iq=1,NEQUI
         curref(iq)=PFCEQW(iq)
      enddo

      call curfit_(rc,zc,ncpfc,NECON,WECON,psi_bnd)

      do ik=1,nequi
         ccurx(ik)=PFCEQW(ik)
         pjk(ik)=PFCEQW(ik)
      enddo

      isymm=0 ! symmetric case
      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf
      nstep=0
      ngrid=1

      call eq_0(pjk,psk,ncequi,nstep,ngrid,
     &          alf0, alf1, alf2, bet0, bet1, bet2,
     &          betpol, betplx, zli3,
     &          ngav1,
     &          ftok, tokout, psicen, pscout, 
     &          nursb,psi_bnd,alp_b,rax,zax,n_ctrl,b0,r0)
      call rdexf(ncequi)
      call eq_ax(pjk, psk, ncequi, nstep,ngrid,
     &           alf0, alf1, alf2, bet0, bet1, bet2,
     &           betpol,  betplx, zli3,
     &           ngav1,
     &           ftok,tokout,psicen,pscout,
     &           EREVE0, ERPS,
     &           psi_bnd,alp_b,rax,zax,isymm)
      call wrd
      call prefit(rc,zc,ncpfc,NECON,WECON,rax,zax,alp_b,psi_bnd)

      ngrid=1
      EREVE0=1.d-7

      do

         nstep=nstep+1

         if(kxwx.eq.1)then
            call precal(rc,zc,ncpfc,NECON, WECON )
         elseif(kxwx.eq.0)then
            call precal_wx(rc,zc,ncpfc,NECON, WECON )
         elseif(kxwx.eq.2)then
            call precal(rc,zc,ncpfc,NECON, WECON )
         endif

         if(n_ctrl.eq.1 .OR. kastr.eq.1) then  !limiter point
            call curfit_L(rc,zc,ncpfc,NECON,WECON,psi_bnd)
         else
            call curfit(rc,zc,ncpfc,NECON,WECON,psi_bnd)
         endif

         do ik=1,nequi
            ccurx(ik)=PFCEQW(ik)
            pjk(ik)=PFCEQW(ik)
         enddo

         call eq_ax(pjk, psk, ncequi, nstep,ngrid,
     &              alf0, alf1, alf2, bet0, bet1, bet2,
     &              betpol,  betplx, zli3,
     &              ngav1,
     &              ftok,tokout,psicen,pscout,
     &              EREVE0, ERPS,
     &              psi_bnd,alp_b,rax,zax,isymm )

         do L=1,NCEQUI
            PJKp1(L) = PJK(L)
            PSKP1(L) = PSK(L)
            PSKm1(L) = PSK(L)
         enddo

         call wrd

        if(ERPS.le.EREVE0) EXIT

      enddo

      call wrd

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'psiplcoils.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
         write(1,*) psi_bnd,tokout,ftok
      close(1)

      write(fname,'(a,a)') path(1:kname),'pfc_curr.wr'
      open(1,file=fname,form='formatted')
         write(1,*) 'coil currents [mA]'
         write(1,'(a6,e13.5)') 'PF1',pjk(1) !*1.d3
         write(1,'(a6,e13.5)') 'PF2',pjk(2) !*1.d3
         write(1,'(a6,e13.5)') 'F3,4',pjk(3) !*1.d3
         write(1,'(a6,e13.5)') 'F5,6',pjk(4) !*1.d3
         write(1,'(a6,e13.5)') 'CS',pjk(5) !*1.d3
      close(1)

      write(fname,'(a,a)') path(1:kname),'tcurrs.wr'
      open(1,file=fname,form='formatted')
         do ik=1,nequi
            do j=1,NPFC
               if(ik .eq. NEPFC(j)) then
                  write(1,'(1E25.11)') PFCEQW(ik)
                  EXIT
               endif
            enddo
         enddo
      close(1)

      i_bsh=-1
      nstep  = 0
      keyctr =0
      i_betp =0
      i_eqdsk =0

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,e_psi)

      return
      end subroutine cf_init_full

!-------------------------------------------------------------------
      subroutine cf_init_bkwd(k_auto, nstep, dt, time,
     &                        voltpf, d_pf_mat, d_tcam_mat)

      use durs_d_modul
      use sp_parameters
      use iopath, only: kname, path
      use parcur, only: curref
      use comevl
      use keys

      implicit none

      integer, intent(in) :: k_auto
      integer, intent(out) :: nstep
      real*8, intent(in), dimension(*) :: voltpf, d_pf_mat, d_tcam_mat

      integer :: i, j, k, l, ncequi, i_bsh, ik, iq, ncpfc, ngrid, nvv, 
     &   nbp, nfw, nc, nloop, n_ctrl, ngav1, ngra1, nprob, nursb,
     &   nout, nter, isymm, ninfw, ninf, nequiold, nbpold
      integer, dimension(nclim) :: ntype
      integer, dimension(nilim) :: necon
      real*8 :: alp_b, psi_bnd, psi0_bnd, platok, dt, time, ereve0,
     &   erps, e_psi, pscout, tokout, zli3, betpol, ftok, psicen
      real*8, dimension(nclim) :: pc, psip, vc, hc, ccurx, ccury,
     &   rc, rc1, rc2, rc3, rc4, zc, zc1, zc2, zc3, zc4, 
     &   alw0, alw1, alw2 
      real*8, dimension(nilim) :: wecon
      real*8, dimension(njlim) :: volk, volkp1, pjk, pjkp1, pjkp, pjkd,
     &   psk, pskp1, pskp, pskm1
      real*8, dimension(nloopp) :: rloop, zloop
      real*8, dimension(nprobp) :: rprob, zprob, fiprob
      real*8, dimension(njlim, njlim) :: res

      character(len=40) :: eqdfn
      character(len=80) :: fname

      common /comeqg/ ncequi
      common /comsta/ platok, eqdfn, i_bsh
      common /comst0/ res, volk, volkp1
      common /comst1/ PJK, PJKP1, PJKP, PJKD
      common /comst2/ PSK,  PSKP1, PSKP, PSKM1
      common /comloo/ rloop, zloop, rprob, zprob, fiprob, nloop, nprob

      kstep = nstep
      ninf=1    
      kxwx=1  ! #prescribed x-points;
      ksnf=0  ! if ksnf=1,then kxwx must be =1
      i_bsh=1

      call noauto
      call rd_ppind
      call rd_prob(NPROb, RPROb, ZPROb,  FIPROb)
      call rd_loop(NLOOp, RLOOp, ZLOOp)
      call rdcoil(nc,ncpfc,rc,zc,pc,necon,wecon)

      write(*,*) 'inizia fit', nequi, ncequi, ncpfc, necon, wecon, nc

      do L=1,NCEQUI
         if ( L.LE.NCEQUI ) then
            PFCEQW(L)=PJK(L)  
         else
            PC(NCPFC+L-NEQUI)= PJK(L) 
         endif
      enddo
      nequiold=nequi
      nbpold=ncequi-nequiold
      nequi=ncequi	
      do k=ncpfc+1,nc
         necon(k)=k-ncpfc+nequiold
         wecon(k)=1.
      enddo
      ncpfc=nc

      write(*,*) 'cf_ini_bkdw', nequi, pjk(1:ncequi)
      ngrid=1

      write(*,*) 'inizia fit 2', nequi, ncequi, ncpfc, necon, wecon, nc

      platok = tokf
      nstep  = 0

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,psi0_bnd)

      write(*,*) nequi,pjk(1:ncequi)
      call prefit(rc,zc,nc,NECON,WECON,rax,zax,alp_b,psi_bnd)
      write(*,*) nequi,pjk(1:ncequi)

      do iq=1,NEQUI
         curref(iq)=PFCEQW(iq)
      enddo

      write(*,*) nequi,curref(1:nequi)
      call curfit_(rc,zc,nc,NECON,WECON,psi_bnd)
      write(*,*) nequi,PFCEQW(1:nequi)

      do ik=1,nequi
         ccurx(ik)=PFCEQW(ik)
         pjk(ik)=PFCEQW(ik)
      enddo

      write(*,*) nequi,pjk(1:ncequi)

      isymm=0 ! symmetric case

      nursb=nurs
      psicen=psax
      ngav1=0
      ftok=tokf
      nstep=0
      ngrid=1

      call eq_0(pjk,psk,ncequi,nstep,ngrid,
     &          alf0, alf1, alf2, bet0, bet1, bet2,
     &          betpol, betplx, zli3,
     &          ngav1,
     &          ftok, tokout, psicen, pscout, 
     &          nursb,psi_bnd,alp_b,rax,zax,n_ctrl,b0,r0 )

      write(*,*) nequi, ncequi
      call rdexf(ncequi)
      nequi=ncequi
      write(*,*) nequi, ncequi
      call eq_ax(pjk, psk, ncequi, nstep,ngrid,
     &           alf0, alf1, alf2, bet0, bet1, bet2,
     &           betpol,  betplx, zli3,
     &           ngav1,
     &           ftok,tokout,psicen,pscout,
     &           EREVE0, ERPS,
     &           psi_bnd,alp_b,rax,zax,isymm )

      call wrd
      call prefit(rc,zc,nc,NECON,WECON,rax,zax,alp_b,psi_bnd)

      ngrid=1
      EREVE0=1.d-7

      do

         nstep=nstep+1

         if(kxwx.eq.1)then
            call precal(rc,zc,nc,NECON, WECON )
         elseif(kxwx.eq.0)then
            call precal_wx(rc,zc,nc,NECON, WECON )
         elseif(kxwx.eq.2)then
            call precal(rc,zc,nc,NECON, WECON )
         endif

         if(n_ctrl.eq.1 .OR. kastr.eq.1) then  !limiter point
            call curfit_L(rc,zc,nc,NECON,WECON,psi_bnd)
         else
            call curfit(rc,zc,nc,NECON,WECON,psi_bnd)
         endif

         do ik=1,nequi
            ccurx(ik)=PFCEQW(ik)
            pjk(ik)=PFCEQW(ik)
         enddo

         write(*,*) 'iter', nequi,pjk(1:ncequi)

         call eq_ax(pjk, psk, ncequi, nstep,ngrid,
     &              alf0, alf1, alf2, bet0, bet1, bet2,
     &              betpol,  betplx, zli3,
     &              ngav1,
     &              ftok,tokout,psicen,pscout,
     &              EREVE0, ERPS,
     &              psi_bnd,alp_b,rax,zax,isymm )

         do L=1,NCEQUI
            PJKp1(L) = PJK(L)
            PSKP1(L) = PSK(L)
            PSKm1(L) = PSK(L)
         enddo

         call wrd

         if(ERPS.le.EREVE0) EXIT

      enddo

      call wrd

      write(fname,'(a,a)') path(1:kname),'currents.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
      close(1)

      write(fname,'(a,a)') path(1:kname),'psiplcoils.wr'
      open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
         write(1,*) psi_bnd
      close(1)
       
      write(fname,'(a,a)') path(1:kname),'pfc_curr.wr'
      open(1,file=fname,form='formatted')
         write(1,*) 'coil currents [mA]'
         write(1,'(a6,e13.5)') 'PF1',pjk(1) !*1.d3
         write(1,'(a6,e13.5)') 'PF2',pjk(2) !*1.d3
         write(1,'(a6,e13.5)') 'F3,4',pjk(3) !*1.d3
         write(1,'(a6,e13.5)') 'F5,6',pjk(4) !*1.d3
         write(1,'(a6,e13.5)') 'CS',pjk(5) !*1.d3
      close(1)

      write(fname,'(a,a)') path(1:kname),'tcurrs.wr'
      open(1,file=fname,form='formatted')
         do ik=1,nequi
            do j=1,NPFC
               if(ik .eq. NEPFC(j)) then
                  write(1,'(1E25.11)') PFCEQW(ik)
                  EXIT
               endif
            enddo
      	 enddo
      close(1)

      i_bsh=-1
      nstep  = 0
      keyctr =0
      i_betp =0
      i_eqdsk =0

      call eqb(alf0,alf1,alf2, bet0,bet1,bet2, alw0,alw1,alw2,
     &         betplx, i_betp,
     &         keyctr, nstep, platok, rax,zax, b0,r0, psax, igdf,
     &         n_tht, n_psi, epsro, nurs, i_eqdsk,i_bsh,
     &         psi_bnd,e_psi)

      nequi=nequiold
      ncpfc=nc-nbpold
      do k=ncpfc+1,nc
         necon(k)=0
         wecon(k)=0.
      enddo

      write(*,*) nequi, ncequi, pjk(1:ncequi)

      return
      end subroutine cf_init_bkwd

!--------------------------------------------
      subroutine aspid_flag(k_astr)

      use keys, only: kastr

      implicit none

      integer, intent(in) :: k_astr

      kastr=k_astr

      return   
      end subroutine aspid_flag
