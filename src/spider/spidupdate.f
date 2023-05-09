      subroutine SPIDUPDATE(machine, coilzzz, time_nowz, nccc)

      use durs_d_modul
      use sp_parameters
      use iopath, only: path
      use comtim
      use comevl
      use keys, only: kstep, kpr

      implicit none

      integer, intent(in) :: nccc
      character(len=4), intent(in) :: machine
      double precision, intent(in) :: time_nowz
      double precision, intent(out), dimension(nccc) :: coilzzz

      integer :: i, j, l, ncequi, i_bsh, i_tim, ipsmk2, eq_cmd
      real*8 :: camtok, platok, dt, time, r_ax, z_ax, psi_b,
     &   psi_eav, betpol, tokout, psiout
      real*8, dimension(njlim) :: VOLK, VOLKP1, PJK, PJKP1, PJKP, PJKD,
     &   PSK, PSKP1, PSKP, PSKM1
      real*8, dimension(njlim, njlim) :: RES
      character(len=40) :: eqdfn
      character(len=80) :: fname

      common/comeqg/  ncequi
      common /comsta/ platok,eqdfn,i_bsh
      common /com234/ betpol,tokout,psiout
      common /com_cam/ camtok
      common/comst0/ RES, VOLK,  VOLKP1
      common/comst1/ PJK, PJKP1, PJKP, PJKD
      common/comst2/ PSK, PSKP1, PSKP, PSKM1
      common /timingcmdipsmk/ ipsmk2,eq_cmd

      do L=1,NPFC
         PFCUR1(L) = PFCUR2(L)
         PFCW1(L)  = PFCW2(L)
         PFCD1(L)  = PFCD2(L)
      enddo

      camtok=0.d0
      do i=nequi+1,ncequi
         camtok=camtok+PJKP1(i)
      enddo

      call get_tim(dt,time)

      DT_EF = dt
      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/dpsipldt.wr'
         open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)-psk(j))/dt_EF ,j=1,NCEQUI)
         close(1)
      endif

! mutual inductances plasma to coils
      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/psi_to_coils.wr'
         open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
         close(1)
      endif

      do L=1,NCEQUI
         PJK(L)   = PJKP1(L)
         VOLK(L)  = VOLKP1(L)

	 if (ipsmk2.lt.1) then
            PSKM1(L) = PSK(L)
            PSK(L)   = PSKP1(L)
	 endif
      enddo

      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/currents.wr'
         open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)(pjk(j),j=1,NCEQUI)
         close(1)
      endif
      if (machine.eq.'aug '.or.machine.eq.'aug_'
     & .or.machine.eq.'aug5') then 
         coilzzz(1)=pjk(1)
         coilzzz(2)=pjk(7)
         coilzzz(3)=pjk(8)
         coilzzz(4)=pjk(6)
         coilzzz(5)=pjk(5)
         coilzzz(6)=pjk(4)
         coilzzz(7)=pjk(2)-pjk(1)
         coilzzz(8)=pjk(3)-pjk(2)
         coilzzz(9)=pjk(9)
         coilzzz(10)=pjk(10)
         coilzzz(11)=pjk(11)
         coilzzz(12)=pjk(12)

         coilzzz(1:12)=coilzzz(1:12)*1.e3
      else	
         coilzzz(1:nccc)=pjk(1:nccc)
         coilzzz(1:nccc)=coilzzz(1:nccc)*1.e3
      endif

      if(kstep.lt.nstep_p) then
         i_tim=kstep+1

         call get_psib(r_ax,z_ax,psi_b)
         call psib_ext(psi_eav)

         time_t(i_tim)=time
         torcur(i_tim)=platok
         rm_t(i_tim)=r_ax
         zm_t(i_tim)=z_ax
         rxp_t(i_tim)=rxpnt
         zxp_t(i_tim)=zxpnt
         betpol_t(i_tim)=betpol
         bettor_t(i_tim)=psi_eav
         psim_t(i_tim)=psax+psi_b
         psib_t(i_tim)=psi_b
      endif

      return
      end subroutine SPIDUPDATE

!----------------------------------------------------------------
      subroutine f_SPIDUPDATE(machine,coilzzz,time_nowz,nccc)

      use durs_d_modul       
      use sp_parameters
      use iopath, only: path
      use comtim
      use comevl, only: dt_ef, nequi
      use keys, only: kstep, kpr

      implicit none

      integer, intent(in) :: nccc
      character(len=4), intent(in) :: machine
      double precision, intent(in) :: time_nowz
      double precision, intent(out), dimension(nccc) :: coilzzz

      integer :: i, j, l, ncequi, i_bsh, i_tim, ipsmk2, 
     &   eq_cmd, numwr
      real*8 :: camtok, platok, dt, time, r_ax, z_ax, psi_b,
     &   betpol, betful, tokout, psiout
      real*8, dimension(njlim) :: VOLK, VOLKP1, PJK, PJKP1, PJKP, PJKD,
     &   PSK, PSKP1, PSKP, PSKM1
      real*8, dimension(njlim, njlim) :: RES
      character(len=40) :: eqdfn 
      character(len=80) :: fname

      common/comeqg/  ncequi
      common /comsta/ platok,eqdfn,i_bsh
      common/comst0/ RES, VOLK, VOLKP1
      common/comst1/ PJK, PJKP1, PJKP, PJKD
      common/comst2/ PSK, PSKP1, PSKP, PSKM1
      common /timingcmdipsmk/ ipsmk2,eq_cmd
      common /com_cam/ camtok

      save numwr

      do L=1,NCEQUI
         PJK(L)  = PJKP1(L)
         VOLK(L) = VOLKP1(L)
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

      DT_EF = dt
      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/dpsipldt.wr'
         open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)-psk(j))/dt_EF ,j=1,NCEQUI)
         close(1)
      endif

! mutual inductances plasma to coils
      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/psi_to_coils.wr'
         open(1,file=fname,form='formatted')
         write(1,*) NEQUI,NCEQUI
         write(1,*)( (pskp1(j)+psk(j))/2. ,j=1,NCEQUI)
         close(1)
      endif

      if(kpr .ge. 0) then
         write(fname,'(a,a)') TRIM(path), '/currents.wr'
         open(1,file=fname,form='formatted')
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
         coilzzz(2)=pjk(7)
         coilzzz(3)=pjk(8)
         coilzzz(4)=pjk(6)
         coilzzz(5)=pjk(5)
         coilzzz(6)=pjk(4)
         coilzzz(7)=pjk(2)-pjk(1)
         coilzzz(8)=pjk(3)-pjk(2)
         coilzzz(9)=pjk(9)
         coilzzz(10)=pjk(10)
         coilzzz(11)=pjk(11)
         coilzzz(12)=pjk(12)
         coilzzz(1:12)=coilzzz(1:12)*1.e3
      else   
         coilzzz(1:12)=pjk(1:12)
         coilzzz(1:12)=coilzzz(1:12)*1.e3
      endif

      if(kstep.lt.nstep_p)then
         i_tim=kstep+1
         time_t(i_tim)=time
         torcur(i_tim)=platok
         rm_t(i_tim)=rax
         zm_t(i_tim)=zax
         rxp_t(i_tim)=rxpnt
         zxp_t(i_tim)=zxpnt
         betpol_t(i_tim)=betpol
         bettor_t(i_tim)=betful
         psim_t(i_tim)=psax
         psib_t(i_tim)=psbo
      endif

      return
      end subroutine f_SPIDUPDATE
