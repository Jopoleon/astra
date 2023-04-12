      subroutine rrsudg(time_ext,dt_smlk)

	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
      use fenix_params
			
      implicit none


      real*8 time_ext,dt_smlk,gvcoil(15)
      real*8 vcoiltmp(10),dumz81

      save vcoiltmp



      if (MACHINE.eq.'dem_') then
          call shmr( &
     &        CPEL1, CIMP3, CV4,  CBND3, &
     &        CSCL1, ZRD71, ZRD73, CDYM3, &
     &        CDWM5, CDWM6, &
     &        CDJM5, CDJM6, CDJM7,  CDJM8, ZRD70, &
     &        CV13, CSOL1, &
     &        CHE1, CDHJ1, CDHJ2, &
     &        CHE3, CDHJ3, CDHJ4, &
     &        CV3, CDHJ5, &
     &        CDMJ1, CDMJ2, CDMJ3, CDMJ4,gvcoil(1:15), &
     &        time_ext, CV6, CDVM7)

		!	 gvcoil(12:15)=0.
!	CV3 = 30. !test EFable, uncontrolled Pfus
!	CPEL1=CPEL1/10
!	CIMP3=CIMP3/10
	
          CV13 = MAX(1., CV13)  ! finite pump speed to avoid NaN
          dt_smlk = CDVM7  ! simulink tau defined in equ log


      endif  ! (MACHIN.eq.'dem_')





      if (MACHINE.eq.'aug_') then
          call shmr( &
     &        CPEL1,CV13, &
     &        CDMJ1, CDMJ2, CDMJ3, CDMJ4, &
     &        ZRD84, &
     &        CAR32(1:8), CAR32(9:16), CAR32(17:24), CAR32(25:26), &
     &        vcoiltmp(1:10),CAR33(1:24),BTOR, &
     &        time_ext, CV6,CDVM7)


!	write(9871,*) 'rrsug', &
!     &        CPEL1,CV13, &
!     &        CDMJ1, CDMJ2, CDMJ3, CDMJ4, &
!     &        ZRD84, &
!     &        CAR32(1:8), CAR32(9:16), CAR32(17:24), CAR32(25:26), &
!     &        vcoiltmp(1:10),CAR33(1:24),BTOR, &
!     &        time_ext, CV6,CDVM7


		 BTOR=1.9 !TEST EFABLE to be removed
		 write(*,*) 'remove this 2 lines'

	BTOR=abs(BTOR) !Btor defined here absolute value. sign has to be given separatly
!	BTORX      =   BTOR
          dt_smlk = CDVM7 ! simulink tau defined in equ log
          ! geom1d(64)=WTOZR(ROC) !total Wmhd including fast ions  in MJ
          ! geom1d(72)=sum(NE(1:NA1))/NA1 !H-1 1019 m-3
          ! car32(25:32)=0.
          ! car32(32)=min(7.5,1.5*1.3*geom1d(64)*geom1d(72)/1e19/1e6)

          vcoil(1:10)=vcoiltmp(1:10)
          if (btipdirec.eq.-1) then
              ! vcoil=-vcoil
          endif
          vcoil(11:12)=0.

          ZRD93 = time_ext+dt_smlk
      endif  ! (MACHIN.eq.'aug_')

      end subroutine rrsudg
