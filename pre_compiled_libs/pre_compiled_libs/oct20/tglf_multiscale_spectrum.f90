!
!--------------------------------------------------------------
!
      SUBROUTINE get_multiscale_spectrum
!
!***********************************************************************
!  questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!
!  TGLF multiscale saturation rule: sat_rule_in = 1 option
!  April 8, 2016: G.M. Staebler, J. Candy, N. T. Howard, and C. Holland
!                  Physics of Plasmas, 23 (2016) 062518
!  June 22, 2017: Retuned after coding error to Laplacian terms in Ampere 
!                 and Poisson equations was fixed
!  TGLF multiscale saturation rule: sat_rule_in = 2 option
!  Aug. 6, 2020 based on papers submitted to PPCF & NF
!
!***********************************************************************
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_xgrid
      IMPLICIT NONE
      !
      LOGICAL :: USE_MIX=.TRUE.
      LOGICAL :: USE_X3=.FALSE.
      LOGICAL :: USE_X4=.FALSE.
      LOGICAL :: first_pass = .TRUE.
      LOGICAL :: USE_SUB1=.FALSE.
      LOGICAL :: USE_PRESSURE = .TRUE.
      INTEGER :: i,is,k,j,j1,j2,jmax1,jmax2
      
      INTEGER :: expsub=2, exp_ax, jmax_mix
      REAL :: test,testmax1,testmax2
      REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
!      REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
      REAL :: gamma0,gamma1,gammaeff,dgamma
      REAL :: cnorm, phinorm, czf, cz1, cz2, kyetg
      REAL :: kycut
      REAL :: cky,sqcky,delta,ax,ay,kx
      REAL :: mix1,mix2,mixnorm,gamma_ave
      REAL :: vzf,dvzf,vzf1,vzf2
      REAL :: vzf_mix, kymax_mix
      REAL :: etg_streamer
      REAL :: kxzf, sat_geo_factor
      REAL :: b0,b1,b2,b3
      REAL :: d1,d2,Gq,kx_width,dlnpdr,ptot
      REAL :: measure
      REAL,DIMENSION(nkym) :: gamma_net=0.0
      REAL,DIMENSION(nkym) :: gamma=0.0
      REAL,DIMENSION(nkym) :: gamma_kymix=0.0
      REAL,PARAMETER :: small=1.0E-10
      !
       
      if(jmax_out.eq.0)then
         first_pass = .TRUE.
      else
        first_pass = .FALSE.
      endif
!
      if(sat_rule_in .lt. 2)then
        measure = SQRT(taus(1)*mass(2))
      else
        measure = SQRT(taus(2)*mass(2)/zs(2)**2)/grad_r0_out
      endif
!
      if(first_pass)then  ! first pass for spectral shift model or only pass for quench rule
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1)
        CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix)
        vzf_out = vzf_mix
        kymax_out = kymax_mix
        jmax_out = jmax_mix
        gamma_net(:) = eigenvalue_spectrum_out(1,:,1)
!        write(*,*)"vzf_out = ",vzf_out,"  jmax_out = ",jmax_out
      endif
      ! model fit parameters
      ! need to set alpha_zf_in = 1.0
      ! Miller geometry values igeo=1
      if(xnu_model_in.eq.3)USE_X3=.TRUE.
      if(xnu_model_in.eq.4)USE_X4=.TRUE.
      if(USE_PRESSURE)then
         dlnpdr = 0.0
         ptot = 0.0
         do is=1,nstotal_in   ! include all species even non-kinetic ones like fast ions
           ptot = ptot + as(is)*taus(is)
           dlnpdr = dlnpdr + as(is)*taus(is)*(rlns(is)+rlts(is))
         enddo
         dlnpdr = rmaj_input*dlnpdr/MAX(ptot,0.01)
         if(dlnpdr .gt. 20.0)dlnpdr = 20.0
         if(dlnpdr .lt. 4.0)dlnpdr = 4.0
      else
         dlnpdr = 12.0
      endif
!         write(*,*)"dlnpdr = ",dlnpdr
!
      czf = ABS(alpha_zf_in)
!
! coefficients for SAT_RULE = 1
      if(sat_rule_in.eq.1) then
         kyetg=1.28
         cnorm=14.21
        if(USE_SUB1)then
          cnorm=12.12
          expsub=1
        endif
        cz1=0.48*czf
        cz2=1.0*czf
        cnorm = 14.29
        if(USE_X3)then
         cnorm = 12.94  ! note this is normed to GASTD CGYRO units
         cz1 = 0.0
         cz2=1.4*czf
         etg_streamer = 1.0
         kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
         cky=3.0
         sqcky=SQRT(cky)
        endif
        if(USE_MIX)then
!original        kyetg=1.9*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
          cky=3.0
          sqcky=SQRT(cky)
          if(USE_SUB1)cnorm=12.12
          cz1=0.48*czf
!original         cz2=0.92*czf
! retuned June 22,2017
          cz2 = 1.0*czf
          etg_streamer=1.05
          if(alpha_quench_in .ne. 0.0)etg_streamer=2.1
          kyetg=etg_streamer*ABS(zs(2))/SQRT(taus(2)*mass(2))  ! fixed to ion gyroradius
        endif
        if(igeo.eq.0)then ! s-alpha
           cnorm=14.63
           cz1=0.90*czf
           cz2=1.0*czf
        endif
     endif
! coefficents for SAT_RULE = 2
     if(sat_rule_in.eq.2)then
       ! SAT2 fit for CGYRO linear modes
        b0 = 0.72
        b1 = 1.22
!        b2 = 24.52  ! note this is b2**2 in PPCF paper 2020
        b2 = 11.21
!        b3 = 0.88
        b3 = 2.4
!        d1 = 0.5475*((Bt0_out/B_geo0_out)**4)/grad_r0_out**2
        d1 = (Bt0_out/B_geo0_out)**4    ! PPCF paper 2020
        Gq = B_geo0_out/grad_r0_out
        d2 = b3/Gq**2
        cnorm = b2*(q_in/2.0)*(12.0/dlnpdr)
        kyetg = 1000.0   ! does not impact SAT2
        cky=3.0
        sqcky=SQRT(cky)
        kycut = b0*kymax_out
        cz1 = 1.1*czf
      endif
!      write(*,*)"Bt0 = ",Bt0_out," B0 = ",B_geo0_out," gradr0 = ",grad_r0_out
!      write(*,*)"G1 = ",sat_geo1_out,"  G2 = ",sat_geo2_out,"  sat_geo0 = ",sat_geo0_out
!      write(*,*)"d1 = ",d1,"  d2 = ",d2
! coefficients for spectral shift model for ExB shear
      ax=0.0
      ay=0.0
      exp_ax = 1
      if(alpha_quench_in.eq.0.0)then
      !spectral shift model parameters
        ax = 1.15
        ay = 0.56
        exp_ax = 4
       if(sat_rule_in.eq.2)then
         ax = 1.55
         ay = 1.0
         exp_ax = 2
       endif
      endif
!    write(*,*)" ax= ",ax," ay= ",ay," kycut = ",kycut
!  apply the spectral shift temporal suppression factor
!
     if(first_pass .eqv. .FALSE.)then  ! second pass for spectral shift model
        do i=1,nky
          kx = spectral_shift_out(i)
          if(sat_rule_in.eq.2)then
            ky0=ky_spectrum(i)
            if(ky0.lt.kycut)then
              kx_width = kycut/grad_r0_out
            else
              kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
            endif
            kx = kx*ky0/kx_width
          endif
          gamma_net(i) = eigenvalue_spectrum_out(1,i,1)/(1.0 + (ax*kx)**exp_ax)
!         write(*,*)ky0,kx,(ax*kx)**4
!         write(*,*)i,"gamma_net = ",gamma_net(i)
        enddo
!     write(*,*)"gammamax_out = ",vzf_out*kymax_out, gamma_net(jmax_out),"  vexb_shear = ",vexb_shear_s!!
! find the maximum of gamma_net/ky
        CALL get_zonal_mixing(nky,ky_spectrum,gamma_net,vzf_mix,kymax_mix,jmax_mix)
        vzf_out = vzf_mix
        kymax_out = kymax_mix
        jmax_out = jmax_mix
      endif   ! second pass complete
! compute multi-scale phi-intensity spectrum field_spectrum(2,,) = phi_bar_out
      ! note that the field_spectrum(1,,) = v_bar_out = 1.0 for sat_rule_in = 1
      gammamax1= vzf_out*kymax_out
      kymax1 = kymax_out
      jmax1 = jmax_out
      vzf1 = vzf_out
!      vzf2 = gammamax2/kymax2
!      dvzf = MAX(vzf2-vzf1,0.0)
!      dgamma= dvzf*kymax1
!      gamma1=gammamax1
!      write(*,*)"dvzf = ",dvzf," vzf1 = ",vzf1," vzf2 = ",vzf2
!      write(*,*)"gammamax1 = ",gammamax1
!      write(*,*)"kymax1 = ",kymax1,"  kycut = ",kycut
      do j=1,nky
! include zonal flow effects on growthrate model
!          gamma=0.0
          gamma0 = gamma_net(j)
          ky0=ky_spectrum(j)
          if(sat_rule_in.eq.1)then
            if(ky0.lt.kymax1)then
              gamma(j) = MAX(gamma0  - cz1*(kymax1 - ky0)*vzf1,0.0)
            else
              gamma(j) = cz2*gammamax1  +  Max(gamma0 - cz2*vzf1*ky0,0.0)          
            endif
          elseif(sat_rule_in.eq.2)then
! note that the gamma mixing is for gamma_model not G*gamma_model
            if(ky0.lt.kymax1)then
              gamma(j) = gamma0
            else
              gamma(j) = gammamax1 +Max(gamma0 - cz1*vzf1*ky0,0.0)
            endif
          endif
          gamma_kymix(j) = gamma(j)     ! used if USE_MIX=F 
      enddo
    if(USE_MIX)then
      !mix over ky > kymax with integration weight = sqcky*ky0**2/(ky0**2 + cky*(ky-ky0)**2)
      do j=jmax1+2,nky
        gamma_ave = 0.0
        ky0 = ky_spectrum(j)
        mixnorm = ky0*(ATAN(sqcky*(ky_spectrum(nky)/ky0-1.0))-  &
                  ATAN(sqcky*(ky_spectrum(jmax1+1)/ky0-1.0)))
        do i=jmax1+1,nky-1
          ky1 = ky_spectrum(i)
          ky2 = ky_spectrum(i+1)
          mix1 = ky0*(ATAN(sqcky*(ky2/ky0-1.0))- ATAN(sqcky*(ky1/ky0-1.0)))
          delta = (gamma(i+1)-gamma(i))/(ky2-ky1)
          mix2 = ky0*mix1 + (ky0*ky0/(2.0*sqcky))*(LOG(cky*(ky2-ky0)**2+ky0**2)- &
                 LOG(cky*(ky1-ky0)**2+ky0**2))
          gamma_ave = gamma_ave + (gamma(i)-ky1*delta)*mix1 + delta*mix2
        enddo  
        gamma_kymix(j) = gamma_ave/mixnorm
!        write(*,*)j,ky0,gamma(j),gamma_kymix(j)
      enddo  
    endif      
! intensity model
       do j=1,nky
        gamma0 = eigenvalue_spectrum_out(1,j,1)
        ky0 = ky_spectrum(j)
        kx = spectral_shift_out(j)
        if(sat_rule_in.eq.1)then
          sat_geo_factor = SAT_geo0_out
          kx_width = ky0
        endif
        if(sat_rule_in.eq.2)then
          if(ky0.lt. kycut)then
            kx_width = kycut/grad_r0_out
            sat_geo_factor = SAT_geo0_out*d1*SAT_geo1_out
           else
             kx_width = kycut/grad_r0_out + b1*(ky0 - kycut)*Gq
             sat_geo_factor = SAT_geo0_out*(d1*SAT_geo1_out*kycut +  &
                              (ky0 - kycut)*d2*SAT_geo2_out)/ky0
           endif
           kx = kx*ky0/kx_width
    !       write(*,*)"sat2 ",B_geo0_out,grad_r0_out,SAT_geo1_out,SAT_geo2_out,sat_geo_factor
        endif
        do i=1,nmodes_in
          gammaeff = 0.0
          if(gamma0.gt.small)gammaeff = &
               gamma_kymix(j)*(eigenvalue_spectrum_out(1,j,i)/gamma0)**expsub
          if(ky0.gt.kyetg)gammaeff = gammaeff*SQRT(ky0/kyetg)
          field_spectrum_out(2,j,i) = measure*cnorm*((gammaeff/(kx_width*ky0))/(1.0+ay*kx**2))**2
          if(units_in.ne.'GYRO')field_spectrum_out(2,j,i) = sat_geo_factor*field_spectrum_out(2,j,i)
        enddo
     enddo
     ! recompute the intensity and flux spectra
      do j=1,nky
         do i=1,nmodes_in
            phinorm=field_spectrum_out(2,j,i) 
            field_spectrum_out(1,j,i) = phinorm
            field_spectrum_out(3,j,i) = QL_field_spectrum_out(3,j,i)*phinorm
            field_spectrum_out(4,j,i) = QL_field_spectrum_out(4,j,i)*phinorm
            do is=1,ns
               intensity_spectrum_out(1,is,j,i) = QL_intensity_spectrum_out(1,is,j,i)*phinorm
               intensity_spectrum_out(2,is,j,i) = QL_intensity_spectrum_out(2,is,j,i)*phinorm
               intensity_spectrum_out(3,is,j,i) = QL_intensity_spectrum_out(3,is,j,i)*phinorm
               intensity_spectrum_out(4,is,j,i) = QL_intensity_spectrum_out(4,is,j,i)*phinorm
               do k=1,3
                  flux_spectrum_out(1,is,k,j,i) = QL_flux_spectrum_out(1,is,k,j,i)*phinorm
                  flux_spectrum_out(2,is,k,j,i) = QL_flux_spectrum_out(2,is,k,j,i)*phinorm
                  flux_spectrum_out(3,is,k,j,i) = QL_flux_spectrum_out(3,is,k,j,i)*phinorm
                  flux_spectrum_out(4,is,k,j,i) = QL_flux_spectrum_out(4,is,k,j,i)*phinorm
                  flux_spectrum_out(5,is,k,j,i) = QL_flux_spectrum_out(5,is,k,j,i)*phinorm
              enddo
           enddo
        enddo
      enddo    

      END SUBROUTINE get_multiscale_spectrum
!
!-------------------------------------------------------------------------------
!

 SUBROUTINE get_zonal_mixing(nmix,ky_mix,gamma_mix,vzf_mix,kymax_mix,jmax_mix)
!
!  finds the maximum of gamma/ky spectrum vzf_out and kymax_out
!
    USE tglf_dimensions
    USE tglf_global
    USE tglf_species
    IMPLICIT NONE
    INTEGER :: nmix,i,k,j,j1,j2,jmax1,jmax2
    REAL :: test,testmax1,testmax2
    REAL :: kx, kyhigh, kycut
    REAL :: gammamax1,kymax1,gammamax2,kymax2,ky0,ky1,ky2
    REAL :: f0,f1,f2,a,b,c,x0,x02,dky,xmax
    REAL :: vzf1, vzf2
    REAL :: kymax_mix, vzf_mix
    INTEGER :: jmax_mix
    REAL, DIMENSION(nkym) :: gamma_mix, ky_mix
!  initialize output of subroutine
    vzf_mix = 0.0
    kymax_mix = 0.0
!
! find the maximum of gamma_mix/ky_mix
!
      gammamax1= gamma_mix(1)
      kymax1 = ky_mix(1)
      testmax1 = gammamax1/kymax1
      testmax2 = 0.0
      jmax1=1
      jmax2=0
      kycut=0.8*ABS(zs(2))/SQRT(taus(2)*mass(2))
      kyhigh=0.15*ABS(zs(1))/SQRT(taus(1)*mass(1))
!      write(*,*)" kycut = ",kycut," kyhigh = ",kyhigh
      j1=1
      j2=1
      ! find the low and high ky peaks of gamma/ky
      do j=2,nky
         ky0 = ky_mix(j)
         if(ky0 .lt. kycut)j1=j1+1
         if(ky0 .lt. kyhigh)j2=j2+1
!         write(*,*)"j=",j,"ky = ",ky0," gamma_net = ",gamma_net(j)
         test = gamma_mix(j)/ky0
         if(ky0 .lt. kycut)then
            if(test .gt. testmax1)then
              testmax1=test
              kymax1 = ky0
              jmax1=j
            endif
         endif
         if(ky0 .gt. kycut)then
           if(test .gt. testmax2)then
             testmax2 = test
             kymax2 = ky0
             jmax2=j
           endif
         endif
      enddo
!      write(*,*)"testmax2 = ",testmax2,"  testmax1*SQRT(kymax1/kymax2) = ",testmax1*SQRT(kymax1/kymax2)
!      write(*,*)"kymax1 = ",kymax1,"  kymax2 = ",kymax2
!      write(*,*)"testmax1 = ",testmax1,"  testmax2 = ",testmax2
!      write(*,*)"j1 = ",j1," j2 = ",j2
      ! handle exceptions
      if(j1.eq.nky)then  ! the maximum ky in the ky-spectrum is less than kycut
         j1=nky-1        ! note that j2=nky in this case
      endif
      if(jmax2.eq.0)jmax2=j2    ! there was no high-k peak set kymax2 to kyhigh or the highest ky
!      if(jmax1.lt.nky)then
!        test=gamma_net(jmax1+1)/ky_spectrum(jmax1+1)
!        if(testmax1.le.test)then
          ! there is no low-k peak
!          write(*,*)" gammamax1/kymax1 = ",gammamax1/kymax1, test
!          jmax1=jmax2
!        endif
!      endif
      gammamax2 = gamma_mix(jmax2)
      kymax2 = ky_mix(jmax2)
      gammamax1 = gamma_mix(jmax1)
      kymax1 = ky_mix(jmax1)
      vzf1 = gammamax1/kymax1
      vzf2 = gammamax2/kymax2
 !     write(*,*)" jmax1 = ",jmax1," jmax2= ",jmax2
 !     write(*,*)" g/k 1 = ",gammamax1/kymax1," g/k 2 = ",gammamax2/kymax2
      !interpolate to find a more accurate low-k maximum gamma/ky
      ! this is cut of at j1 since a maximum may not exist in the low-k range
      if(jmax1.gt.1.and.jmax1.lt.j1)then
!      write(*,*)"refining low-k maximum"
         f0 =  gamma_mix(jmax1-1)/ky_mix(jmax1-1)
         f1 =  gamma_mix(jmax1)/ky_mix(jmax1)
         f2 =  gamma_mix(jmax1+1)/ky_mix(jmax1+1)
         dky = (ky_mix(jmax1+1)-ky_mix(jmax1-1))
         x0 = (ky_mix(jmax1)-ky_mix(jmax1-1))/dky
         a = f0
         x02 = x0*x0
         b = (f1 - f0*(1-x02)-f2*x02)/(x0-x02)
         c = f2 - f0 - b
         xmax = -b/(2.0*c)
         if(xmax .ge. 1.0)then
           kymax1 = ky_mix(jmax1+1)
           gammamax1 = f2*kymax1
         elseif(xmax.lt.0.0)then
           kymax1 = ky_mix(jmax1-1)
           gammamax1 = f0*kymax1
         else
           kymax1 = ky_mix(jmax1-1)+dky*xmax
           gammamax1 = (a+b*xmax+c*xmax*xmax)*kymax1
         endif
      endif
      vzf_mix = gammamax1/kymax1
      kymax_mix = kymax1
      jmax_mix = jmax1
!      write(*,*)"get_zonal_mxing"
!      write(*,*)"gammamax1 = ",gammamax1," kymax1 = ",kymax1

 END SUBROUTINE get_zonal_mixing
