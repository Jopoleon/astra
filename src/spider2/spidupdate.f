!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       SUBROUTINE  SPIDUPDATE(machine,coilzzz,time_nowz,nccc)

         use durs_d_modul       
	use fenix_params
C---------------------------------------------------------------
       INCLUDE 'Descript.inc'
C  ----------
C
       include 'double.inc'
C
       INCLUDE 'prm.inc'
       INCLUDE 'comevl.inc'
       !INCLUDE 'contro.inc'
       !parameter(nstep_p=3000)
C-----------------------------------------------------------------------
       INCLUDE 'comtim.inc'

	 common/comeqg/  ncequi
       common /nostep/ kstep
       common /comhel/ helinp, helout
       
!       common /comsta/ epsro,platok,
!     *                 eqdfn,betplx,tokf,psax,b0,r0,
!     *                 alf0,alf1,alf2,bet0,bet1,bet2,rax,zax,
!     *                 nurs,igdf,n_psi,n_tht,i_bsh,keyctr,i_betp,i_eqdsk
     
       common /comsta/ platok,eqdfn,i_bsh
       
       common /com234/ betpol,tokout,psiout
       common /com_cam/ camtok
	 common /c_kpr/ kpr
       common /key_for_pet/ key_fixbon
       common/com_flag/kastr
C-----------------------------------------------------------------------
C
       DIMENSION  RC(NCLIM),  ZC(NCLIM),  PC(NCLIM), PSIP(NCLIM),
     *            VC(NCLIM),  HC(NCLIM)
C
       DIMENSION  RC1(NCLIM), ZC1(NCLIM), RC2(NCLIM), ZC2(NCLIM),
     *            RC3(NCLIM), ZC3(NCLIM), RC4(NCLIM), ZC4(NCLIM)
C
       INTEGER    NECON(NILIM), NTYPE(NCLIM),nccc
       DIMENSION  WECON(NILIM)
C
       !common/comst0/ RES(NJLIM),  VOLK(NJLIM),  VOLKP1(NJLIM)
       common/comst0/ RES(NJLIM,NJLIM),  VOLK(NJLIM),  VOLKP1(NJLIM)
C
       common/comst1/ PJK(NJLIM),PJKP1(NJLIM),PJKP(NJLIM),PJKD(NJLIM)
       include 'parloo.inc'
       common/comloo/ rloop(nloopp),zloop(nloopp),
     &          rprob(nprobp),zprob(nprobp),fiprob(nprobp),nloop,nprob
       common/comst2/ PSK(NJLIM), PSKP1(NJLIM),PSKP(NJLIM),PSKM1(NJLIM)
       common /com_kout/ key_out
       common /com_enels/ ENELS

       !DIMENSION  cp_com_0(n_cp_m)
C***********************************************************************

       DIMENSION  psiplb(nstep_p),psiexb(nstep_p),psimag(nstep_p),
     *            flu_tor(nstep_p)  
	character*4 machine
	double precision coilzzz(nccc),time_nowz

!	include 'includes_fsim.dat'

	

C***********************************************************************

C       real*8 voltpf(*)
C       real*8 d_pf_mat(*)
C       real*8 d_tcam_mat(*)

C***********************************************************************


	 real*4       a_print(200)
	 character*30 apr
       character*40 eqdfn

C         if(key_out.eq.0) return             

  222        CONTINUE
C
             ERRCU1 = 0.D0
             ERRCU2 = 0.D0
             ERRPS1 = 0.D0
             ERRPS2 = 0.D0
             DELPS1 = 0.D0
             DELPS2 = 0.D0
             EREVE  = 0.D0
             ERPS   = 0.D0
c
             RAXPR  = RM_EF
             ZAXPR  = ZM_EF
            ! RXPPR  = RX0
            ! ZXPPR  = ZX0
            KSTEPR  = KSTEP
            TSTEPR  = TSTEP_EF
C
          DO 293 L=1,NPFC
             PFCUR1(L) = PFCUR2(L)
             PFCW1(L)  = PFCW2(L)
             PFCD1(L)  = PFCD2(L)
  293     CONTINUE

C          do i=1,nequi
C             d_pf_mat(i) = PJKP1(i)*1.0d6  !!!  in [A]
C          enddo
 
             camtok=0.d0
         do i=nequi+1,ncequi
C             d_tcam_mat(i-nequi) = PJKP1(i)*1.0d6  !!!  in [A]
             camtok=camtok+PJKP1(i)
        enddo

      call get_tim(dt,time)
CEFable
	DT_EF = dt
        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'dpsipldt.wr'
            open(1,file=fname,form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)( (pskp1(j)-psk(j))/dt_EF ,j=1,NCEQUI)
            close(1)
        endif


CEfable  mutual inductances plasma to coils
        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'psi_to_coils.wr'
            open(1,file=fname,form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
            close(1)
        endif
CEfable

!	write(6542,'(533E25.11)') time_nowz,
!   & PJK(1:ncequi)


            
          DO 260 L=1,NCEQUI
             PJK(L)   = PJKP1(L)
             VOLK(L)  = VOLKP1(L)
	if (kpr.eq.-2) then
!             PSKM1(L) = PSK(L)
!             PSK(L)   = PSKP1(L)
	endif
	if (kpr.ge.-1) then
             PSKM1(L) = PSK(L)
             PSK(L)   = PSKP1(L)
	endif
  260     CONTINUE


!	write(*,*) 'psicoil flux',PSK(1:ncequi),platok
  
        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'currents.wr'
            open(1,file=fname,form='formatted')
           !open(1,file='currents.wr',form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)(pjk(j),j=1,NCEQUI)
            close(1)
        endif
	if (machine.eq.'aug '.or.machine.eq.'aug_'
     & .or.machine.eq.'aug5') then 
	coilzzz(1)=pjk(1)
	coilzzz(2)=pjk(2)-pjk(1)
	coilzzz(3)=pjk(3)-pjk(2)
	coilzzz(4)=pjk(4)
	coilzzz(5)=pjk(5)
	coilzzz(6)=pjk(6)
	coilzzz(7)=pjk(7)
	coilzzz(8)=pjk(8)
	coilzzz(9)=pjk(9)
	coilzzz(10)=pjk(10)
	coilzzz(11)=pjk(11)
	coilzzz(12)=pjk(12)

	coilzzz(1:12)=coilzzz(1:12)*1.e3
	else
	
	coilzzz(1:nccc)=pjk(1:nccc)
	coilzzz(1:nccc)=coilzzz(1:nccc)*1.e3
	

	endif


C
C***********************************************************************
!        if(kstep.lt.nstep_p .and. kstep.gt.0)then
        if(kstep.lt.nstep_p)then
           i_tim=kstep+1
           
      call get_psib(r_ax,z_ax,psi_b)
      call psib_ext(psi_eav)

!         time_t(i_tim)=time
!         torcur(i_tim)=platok
!         rm_t(i_tim)=r_ax
!         zm_t(i_tim)=z_ax
!         rxp_t(i_tim)=rxpnt_ef
!         zxp_t(i_tim)=zxpnt_ef
!         rxp_t(i_tim)=rxpnt
!         zxp_t(i_tim)=zxpnt
!         betpol_t(i_tim)=betpol
!         bettor_t(i_tim)=psi_eav
!         psim_t(i_tim)=psax+psi_b
!         psib_t(i_tim)=psi_b
        endif
C-----------------------------------------------------------------------
C For writing for "motion picture"
C
C       psiplb(kstep)=pspl_av
C       psiexb(kstep)=psex_av
C       psimag(kstep)=psi0_ax
C       flu_tor(kstep)=flx_fi
C
C       if(kstep.gt.0 .AnD. kstep.le.2000) then
C       kskw=2
C       if(kstep/kskw*kskw .eq. kstep) then
C         numwr=numwr+1
C         call wrdfmv(numwr,timev)
C         call wrdump(numwr,timev,kstep,psiplb,psiexb,psimag,flu_tor)
C       endif
C       endif
C 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
!!!!                   GO TO 100  
C 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100
C-----------------------------------------------------------------------
  103  FORMAT(2X,8E12.5)
 1915  FORMAT(3X,'LEVEL K =',I3,2X,'ITERATION KNEL =',I3,2X,
     *        'TIME(K) =',E12.5)
C********************************


      RETURN
      END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       SUBROUTINE  f_SPIDUPDATE(machine,coilzzz,time_nowz,nccc)

         use durs_d_modul       
	use fenix_params
C---------------------------------------------------------------
       INCLUDE 'Descript.inc'
C  ----------
C
       include 'double.inc'
C
       INCLUDE 'prm.inc'
       INCLUDE 'comevl.inc'
       !INCLUDE 'contro.inc'
       INCLUDE 'comtim.inc'
C-----------------------------------------------------------------------
C
	 common/comeqg/  ncequi
       common /nostep/ kstep
       
!       common /comsta/ epsro,platok,
!     *                 eqdfn,betplx,tokf,psax,b0,r0,
!     *                 alf0,alf1,alf2,bet0,bet1,bet2,rax,zax,
!     *                 nurs,igdf,n_psi,n_tht,i_bsh,keyctr,i_betp,i_eqdsk

       common /comsta/ platok,eqdfn,i_bsh
C-----------------------------------------------------------------------
C
       DIMENSION  RC(NCLIM),  ZC(NCLIM),  PC(NCLIM), PSIP(NCLIM),
     *            VC(NCLIM),  HC(NCLIM)
C
       DIMENSION  RC1(NCLIM), ZC1(NCLIM), RC2(NCLIM), ZC2(NCLIM),
     *            RC3(NCLIM), ZC3(NCLIM), RC4(NCLIM), ZC4(NCLIM)
C
       INTEGER    NECON(NILIM), NTYPE(NCLIM)
       DIMENSION  WECON(NILIM)
C
       !common/comst0/ RES(NJLIM),  VOLK(NJLIM),  VOLKP1(NJLIM)
       common/comst0/ RES(NJLIM,njlim),  VOLK(NJLIM),  VOLKP1(NJLIM)
C
       common/comst1/ PJK(NJLIM),PJKP1(NJLIM),PJKP(NJLIM),PJKD(NJLIM)
       common/comst2/ PSK(NJLIM), PSKP1(NJLIM),PSKP(NJLIM),PSKM1(NJLIM)
       include 'parloo.inc'
       common/comloo/ rloop(nloopp),zloop(nloopp),
     &           rprob(nprobp),zprob(nprobp),fiprob(nprobp),nloop,nprob

        common /com_enels/ ENELS
      !DIMENSION  cp_com_0(n_cp_m)
C***********************************************************************

!       DIMENSION  psiplb(nstep_p),psiexb(nstep_p),
!     *            flu_tor(nstep_p)  
       real*8 errarr(10)

C***********************************************************************

!       real*8 voltpf(*)
!       real*8 d_pf_mat(*)
!       real*8 d_tcam_mat(*)

C***********************************************************************
!	include 'includes_fsim.dat'

       common /com_cam/ camtok
       common /key_for_pet/ key_fixbon

	 common /c_kpr/ kpr
	character*4 machine
	double precision coilzzz(12),time_nowz

	 real*4       a_print(200)
	 character*30 apr
       character*40 eqdfn
       save numwr
       save psi_eav_n

!	write(6542,'(533E25.11)') time_nowz,
!     & PJK(1:52),PJKP1(1:52)-PJK(1:52),
!     & PSK(1:52),PSKP1(1:52)-PSK(1:52)

          do L=1,NCEQUI
             PJK(L)   = PJKP1(L)
             VOLK(L)  = VOLKP1(L)


	if (eq_cmd.eq.1.and.ipsmk2.ge.1) then
!             PSKM1(L) = PSK(L)
!             PSK(L)   = PSKP1(L)
	endif
	if (ipsmk2.lt.1) then
             PSKM1(L) = PSK(L)
             PSK(L)   = PSKP1(L)
	endif

          enddo

          camtok=0.d0
          do i=nequi+1,ncequi
             camtok=camtok+PJKP1(i)
          enddo

      call get_tim(dt,time)
CEFable
	DT_EF = dt
        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'dpsipldt.wr'
            open(1,file=fname,form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)( (pskp1(j)-psk(j))/dt_EF ,j=1,NCEQUI)
            close(1)
        endif


CEfable  mutual inductances plasma to coils
        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'psi_to_coils.wr'
            open(1,file=fname,form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
            close(1)
        endif
CEfable

        if(kpr .ge. 0) then
            write(fname,'(a,a)') path(1:kname),'currents.wr'
            open(1,file=fname,form='formatted')
           !open(1,file='currents.wr',form='formatted')
            write(1,*) NEQUI,NCEQUI
            write(1,*)(pjk(j),j=1,NCEQUI)
            close(1)
        endif

        if(kpr .gt. 0) then
         numwr=numwr+1
         call wr_step(numwr,time,kstep)
        endif

	if (machine.eq.'aug '.or.machine.eq.'aug_'
     & .or.machine.eq.'aug5') then 
	coilzzz(1)=pjk(1)
	coilzzz(2)=pjk(2)-pjk(1)
	coilzzz(3)=pjk(3)-pjk(2)
	coilzzz(4)=pjk(4)
	coilzzz(5)=pjk(5)
	coilzzz(6)=pjk(6)
	coilzzz(7)=pjk(7)
	coilzzz(8)=pjk(8)
	coilzzz(9)=pjk(9)
	coilzzz(10)=pjk(10)
	coilzzz(11)=pjk(11)
	coilzzz(12)=pjk(12)

	coilzzz(1:12)=coilzzz(1:12)*1.e3
	else
	
	coilzzz(1:12)=pjk(1:12)
	coilzzz(1:12)=coilzzz(1:12)*1.e3
	

	endif

!!!!!!!!!!time dependent arrays initialization
 757        continue

        if(kstep.lt.nstep_p)then
           i_tim=kstep+1
!         time_t(i_tim)=time
!         torcur(i_tim)=platok
!         rm_t(i_tim)=rax
!         zm_t(i_tim)=zax
!         rxp_t(i_tim)=rxpnt
!         zxp_t(i_tim)=zxpnt
!         betpol_t(i_tim)=betpol
!         bettor_t(i_tim)=betful
!         psim_t(i_tim)=psax
!         psib_t(i_tim)=psbo
        endif

!!!!!!!!!!time dependent arrays initialization

      RETURN
      END
