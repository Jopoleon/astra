      subroutine spider_run(ncoils, vcoils, equil_in, equil_out, params)

      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters 
      use keys, only: key_0st, key_prs

      implicit none

      integer, intent(in) :: ncoils
      real*8, intent(in) ::  vcoils(ncoils)
      type(type_equilibrium), intent(in)  :: equil_in
      type(type_parameters) , intent(in)  :: params
      type(type_equilibrium), intent(out) :: equil_out

      !locals
      integer nstep, key_dmf, k_grid, k_auto, k_fixfree, key_start
      integer kpr, kname, npointzz
      real*8 time, dt, dpsdt
      integer i, neql
      integer itmastra
      real*8 um, psi_boundary, press0, btor

      npointzz = size(equil_in%profiles_1d%psi)
      nstep     = params%nstep
      time      = params%time
      dt        = params%dt
      key_dmf   = params%key_dmf
      k_grid    = params%k_grid
      k_auto    = params%k_auto
      k_fixfree = params%k_fixfree
      dpsdt     = params%dpsdt
      key_start = params%key_start

      kpr = params%kpr
      call kpr_calc(kpr)
      kname = params%kname 
 
      call put_name(params%prename, kname)
! set key_plc for plasma current
      call put_key_plc(params%key_plc)
! set key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
! set the circuit equation accuracy
      call put_enels(params%enels)
      call put_ndmf(params%n_dmf)
      key_0st = params%key_0stp
      key_prs = params%key_pres

      call put_tim(dt, time)

      itmastra = 0
      if(associated(equil_in%profiles_1d%psi)) then
         itmastra = 1
      endif

      if(itmastra == 1) then
         call itm2spider(equil_in, params)
      endif

      call spider(nstep, time, dt, key_dmf, k_grid, k_auto, k_fixfree, 
     &            dpsdt, key_start, vcoils)  

      call cur_avg
      call wr_spik

! equil_out
      if(k_fixfree == 0) then
         psi_boundary = 0.d0
      else if(k_grid == 0) then
         call get_umup(um, psi_boundary)
      else if(k_grid == 1) then
         call f_get_umup(um, psi_boundary)
      endif
      if(associated(equil_in%profiles_1d%pressure)) then
         press0 = equil_in%profiles_1d%pressure(npointzz)
      endif
      btor = equil_in%global_param%toroid_field%b0

      call get_eq(equil_out, psi_boundary, press0, btor)

      call put_eq(equil_out, params)     

      return
      end subroutine spider_run

!-------------------------------------------------------------
      subroutine itm2spider(equil_in, params)

      use sp_parameters, only: nrp, ntp, n_ursp, twopi
      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters
      use durs_d_modul
      use ppf_modul
      use bnd_modul
      use keys, only: key_0st, key_prs

      implicit none

      type(type_equilibrium) equil_in
      type(type_parameters) params

      integer na1, nbnd, neql, nteta, nstep, key_ini
      real*8 ipl, rtor, btor

      integer, parameter :: nbtabp=1000
      integer i, ib
      real*8 ps0, psb, pscale, fscale
      character*40 eqdfn

      if(.not.associated(equil_in%profiles_1d%psi)) then
          write(*, *) 'itm2astra: INVALID input equilibrium CPO; return'
          pause
          return
      endif
 
! set the key_plc for plasma current
      call put_key_plc(params%key_plc)
! set the key_out for currents and inductive voaltages update in circuit equation
      call put_key_out(params%key_out)
! set the circuit equation accuracy
      call put_enels(params%enels)
        
! input dimension by the size of allocated psi
      na1 = size(equil_in%profiles_1d%psi)
      nbnd = equil_in%eqgeometry%boundary%npoints
      neql    = params%neql
      nteta   = params%nteta
      nstep   = params%nstep
      key_ini = params%key_ini
      key_0st = params%key_0stp
      key_prs = params%key_pres
      i_eqdsk = params%i_eqdsk
      eqdfn   = params%eqdfn

! check dimensions
      if(nbnd > nbtabp) then
         write(*, *) 'spider:nbnd exeeds maximum value', nbtabp
         write(*, *) 'nbnd=', nbnd   
         write(*, *) 'program is interrupted'
         stop      
      endif    
      if(na1 > n_ursp) then
         write(*, *) 'spider:na1 exeeds maximum value', n_ursp  
         write(*, *) 'na1=', na1   
         write(*, *) 'program is interrupted'
         stop      
      endif  
      if(neql > nrp .or. nteta > ntp) then
         write(*, *) 'spider: neql or nteta exeeds maximum value'  
         write(*, *) 'neql=', neql, 'nrp=', nrp   
         write(*, *) 'nteta=', nteta, 'ntp=', ntp   
         write(*, *) 'program is interrupted'
         stop      
      endif
 
! toroidal field, radius, plasma current
      rtor = equil_in%global_param%toroid_field%r0
      btor = equil_in%global_param%toroid_field%b0
      ipl  = equil_in%global_param%i_plasma*1e-6

! p' and ff' profiles
      nutab = na1     
      if(allocated(pstab)) deallocate( pstab, pptab, fptab )
      allocate( pstab(na1), pptab(na1), fptab(na1) )
 
      ps0 = equil_in%profiles_1d%psi(1)
      psb = equil_in%profiles_1d%psi(na1)
      pstab(1) = 0.d0
      do i=2, nutab
         pstab(i) = (equil_in%profiles_1d%psi(i) - ps0)/(psb - ps0)
         if(pstab(i) <= pstab(i-1))   then
            write(*, *) 'psi is nonmonotonic! ', pstab(i), i
         endif
      enddo   
      do i=1, nutab
         pptab(i) = equil_in%profiles_1d%pprime(i)
         fptab(i) = equil_in%profiles_1d%ffprime(i)
      enddo
! Normalization
      pscale = 2 * TWOPI * 1.d-7 * TWOPI
      fscale = TWOPI
      do i=1, nutab
         pptab(i) = -pptab(i)*pscale
         fptab(i) = -fptab(i)*fscale
      enddo

! Boundary
      nbtab = nbnd      
      if(allocated(rbtab)) deallocate( rbtab, zbtab )
      allocate( rbtab(nbtab), zbtab(nbtab) )

      do ib=1, nbtab
         rbtab(ib) = equil_in%eqgeometry%boundary%r(ib)
         zbtab(ib) = equil_in%eqgeometry%boundary%z(ib)
      enddo
   
      if(associated(equil_in%profiles_1d%rho_tor)) then
         call S_profro_p(  equil_in, params)
         call S_profro_sbc(equil_in, params)
         if(associated(equil_in%profiles_1d%pressure)) then
            call S_proffi_pres(equil_in, params)
         endif
      endif

      n_tht  = nteta
      n_psi  = neql
      igdf   = 2
      nurs   = -399
      i_betp = 0
      keyctr = 0 
      i_betp = 0

! fixed boundary equilibrium accuracy

      epsro = params%epsro
      tokf = ipl
      b0 = btor
      r0 = rtor
      rax = 0.5d0*( MAXVAL(rbtab) + MINVAL(rbtab) )
      zax = 0.5d0*( MAXVAL(zbtab) + MINVAL(zbtab) )
      call put_Ipl(tokf)

      if(nstep == 0 .and. key_ini == 0) then
          if(i_betp == 0) then
             betplx = 0.d0
          endif
          psax = ps0
          if(keyctr /= 2) then
             psax = 1.d0
          endif
          if(nurs < 0) then
             alf0 = 0.d0
             alf1 = 0.d0
             alf2 = 0.d0
             bet0 = 0.d0
             bet1 = 0.d0
             bet2 = 0.d0
          endif
      endif

! read eqdsk instead
      if(nstep == 0 .and. key_ini == 0) then
         i_eqdsk = 1
         call tab_efit(tokf, psax, eqdfn, rax, zax, b0, r0)
         keyctr  = 0 
         key_0st = 0
      endif

      if(keyctr == 0) call taburs(0, 1.d0, nurs)

      call aspid_flag(1)
  
      return
      end subroutine itm2spider

!-------------------------------------------------------------
      subroutine put_key_plc(k)

      use keys, only: key_plc
      implicit none

      integer, intent(in) :: k

      key_plc = k

      return
      end subroutine put_key_plc

!--------------------------------------------------------------
      subroutine put_key_out(k)

      use keys, only: key_out

      implicit none

      integer, intent(in) :: k

      key_out = k

      return
      end subroutine put_key_out

!--------------------------------------------------------------
      subroutine put_enels(e)

      implicit none

      real*8, intent(in) :: e
      common /com_enels/ ENELS
      real*8 ENELS

      ENELS = e

      return
      end subroutine put_enels

!--------------------------------------------------------------
      subroutine put_ndmf(k_dmf)

      implicit none

      integer, intent(in) :: k_dmf
      common /com_ndmf/ n_dmf
      integer n_dmf

      n_dmf = k_dmf

      return
      end subroutine put_ndmf

!---------------------------------------------------------------
      subroutine spline_prof(nspl, iplas, rhow, funw, rosn, yout)
! March 2021 git

      use sp_parameters, only: n_ursp

      implicit none

      integer, parameter :: n_ursp4=n_ursp+4, n_ursp6=n_ursp4*6

      integer, intent(in) :: nspl, iplas
      real*8,  intent(in) , dimension(n_ursp) :: funw, rhow
      real*8,  intent(in) , dimension(iplas) :: rosn
      real*8,  intent(out), dimension(iplas) :: yout

      integer :: i, IFAIL
      real*8, dimension(n_ursp4) :: rrk, cck
      real*8 :: zrho, cwk(4), wrk(n_ursp6)

      CALL E01BAF(nspl, rhow, funw, RRK, CCK, 
     &            nspl+4, WRK, 6*nspl+16, IFAIL)
      do i=2, iplas-1
         zrho = rosn(i)
         CALL E02BCF(nspl+4, RRK, CCK, zrho, 0, CWk, IFAIL)
         yout(i) = cwk(1)
      enddo
      yout(1)     = funw(1)
      yout(iplas) = funw(nspl)

      return
      end subroutine spline_prof

!--------------------------------------------------------------
      subroutine S_profro_p(equil_in, parameters_spider)
! April 2013

      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters
      use sp_parameters

      implicit none

      type(type_equilibrium), intent(in) :: equil_in
      type(type_parameters) , intent(in) :: parameters_spider

      integer :: i, nspl
      real*8 :: roc, pscale
      real*8, dimension(nrp) :: rosn
      real*8, dimension(n_ursp)  :: ppr, rhow

      include 'compol.inc'

      iplas = parameters_spider%neql
      nspl = size(equil_in%profiles_1d%psi)

      do i=1, iplas
         rosn(i) = (i - 1.d0)/(iplas - 1.d0)
      enddo

      roc = equil_in%profiles_1d%rho_tor(nspl)

      do i=1, nspl
         rhow(i) = equil_in%profiles_1d%rho_tor(i)/roc
         ppr(i)  = equil_in%profiles_1d%pprime(i)
      enddo 

      call spline_prof(nspl, iplas, rhow, ppr, rosn, dpdpsi(1:iplas))

      pscale = 2 * TWOPI * 1.d-7 * TWOPI

      do i=1, iplas
         dpdpsi(i) = -dpdpsi(i)*pscale
      enddo

      return
      end subroutine S_profro_p

!---------------------------------------------------------------
      subroutine S_profro_sbc(equil_in, parameters_spider)
! April 2013

      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters
      use sp_parameters

      implicit none

      integer :: i, nspl
      real*8 :: roc, btor
      real*8, dimension(n_ursp)  :: funw, rhow
      real*8, dimension(nrp) :: rosn, C_sig, T_el, C_bts, C_driv

      type(type_equilibrium) :: equil_in
      type(type_parameters) :: parameters_spider
      
      include 'compol.inc'

      common /com_sigcd/ C_sig, T_el, C_bts, C_driv

      iplas = parameters_spider%neql
      nspl = size(equil_in%profiles_1d%psi)
  
      do i=1, iplas
         rosn(i) = (i - 1.d0)/(iplas - 1.d0)
      enddo

      roc = equil_in%profiles_1d%rho_tor(nspl)
      btor = equil_in%global_param%toroid_field%b0

      do i=1, nspl
         rhow(i) = equil_in%profiles_1d%rho_tor(i)/roc
      enddo

      do i=1, nspl
         funw(i) = equil_in%profiles_1d%sigmapar%value(i)
      enddo
      call spline_prof(nspl, iplas, rhow, funw, rosn, C_sig(1:iplas))

      do i=1, nspl
         funw(i) = equil_in%profiles_1d%jni%value(i)
      enddo 
      call spline_prof(nspl, iplas, rhow, funw, rosn, C_bts(1:iplas))

      do i=1, nspl
         funw(i) = equil_in%profiles_1d%te%value(i)
      enddo 
      call spline_prof(nspl, iplas, rhow, funw, rosn, T_el(1:iplas))

      do i=1, iplas
         C_bts(i)  = C_bts(i)*btor
         C_driv(i) = 0.d0
         T_el(i)   = T_el(i)*1.d3
      enddo

      return
      end subroutine S_profro_sbc

!----------------------------------------------------------------
      subroutine S_proffi_pres(equil_in, parameters_spider)

      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters
      use sp_parameters

      implicit none

      type(type_equilibrium) :: equil_in
      type(type_parameters) :: parameters_spider

      integer :: i, na1, nspl, knin, knout, kopt, mdamat
      real*8 :: rh0, rh1, rh2, rh3, px0, px1, px2, px3, 
     &   prescale, ptaus, pbclft, pbcrgt, romin
      real*8, dimension(nrp) :: rokn, dpdro, dpdfi
      real*8, dimension(0:nrp) :: P_rho

      real*8, allocatable, dimension(:) :: PXIN, PYIN, PXOUT, PYOUT, 
     &  PYOUTP, PYOUTPP, PY2, PWORK, PYINNEW, p05_prim, rho05, d2pf
      real*8, allocatable, dimension(:, :) :: PAMAT

      include 'compol.inc'
 
      common /com_pres/ P_rho, dpdro, dPdFi, romin

      iplas = parameters_spider%neql
      na1 = size(equil_in%profiles_1d%psi)
      nspl = na1 + 1

      allocate(rho05(na1+1))
      allocate(p05_prim(na1+1))
      allocate(d2pf(iplas))
     
      do i=2, na1
         rho05(i) = 0.5d0*( 
     &      equil_in%profiles_1d%rho_tor(i-1)**2 +
     &      equil_in%profiles_1d%rho_tor(i)**2 ) /
     &      equil_in%profiles_1d%rho_tor(na1)**2
      enddo
      rho05(1)     = 0.d0
      rho05(na1+1) = 1.d0
    
      do i=2, na1
         p05_prim(i) = (
     &      equil_in%profiles_1d%pressure(i) -
     &      equil_in%profiles_1d%pressure(i-1) ) /
     &     (equil_in%profiles_1d%rho_tor(i)**2 -
     &      equil_in%profiles_1d%rho_tor(i-1)**2) *
     &      equil_in%profiles_1d%rho_tor(na1)**2
      enddo
      prescale = 1.d-6
      do i=2, na1
         p05_prim(i) = prescale*p05_prim(i)
      enddo

! Extrapolation on magn. axis
 
      rh0 = 0.d0
      rh1 = rho05(2)
      rh2 = rho05(3)
      rh3 = rho05(4)

      px1 = p05_prim(2)
      px2 = p05_prim(3)
      px3 = p05_prim(4)

      call EXTRP2(rh0, px0, rh1, rh2, rh3, px1, px2, px3)

      p05_prim(1) = px0

! Esxtrapolation to boundry

      rh0 = 1.d0
      rh1 = rho05(na1)     
      rh2 = rho05(na1-1)     
      rh3 = rho05(na1-2)     

      px1 = p05_prim(na1)  
      px2 = p05_prim(na1-1)
      px3 = p05_prim(na1-2)

      call EXTRP2(rh0, px0, rh1, rh2, rh3, px1, px2, px3) 
 
      p05_prim(na1+1) = px0

      do i=1, iplas
         rokn(i) = ((i - 1.d0)/(iplas - 1.d0))**2
      enddo
   
      rokn(iplas) = 1.d0
      rokn(1)     = 0.d0

      KNIN = nspl
      KNOUT = iplas
      KOPT = 1
      PTAUS = 1d-6
      allocate(PYINNEW(KNIN))
      allocate(PY2(KNIN))
      allocate(PWORK(KNIN))
      MDAMAT = 7
      PBCLFT = 2.d0
      PBCRGT = 2.d0
      allocate(PAMAT(MDAMAT, KNIN))
      CALL INTRPTAU(rho05, p05_prim, PYINNEW, PY2, KNIN, rokn, dPdFi, 
     &    d2pf, PYOUTPP, KNOUT, KOPT, PTAUS, PWORK, PAMAT, MDAMAT, 
     &    PBCLFT, PBCRGT)

      return
      end subroutine S_proffi_pres
     
!-----------------------------------------------------------------
      subroutine pres_d_psi

      use sp_parameters
 
      implicit none

      integer, parameter :: nrsp=nrp+1, nrsp4=nrsp+4, nrsp6=nrsp4*6

      integer :: i, ifail, n3spl
      real*8 :: romin, rh0, rh1, rh2, rh3, qx0, qx1, qx2, qx3,
     &   zrho, qkn
      real*8, dimension(nrsp) :: rhocn, funw
      real*8, dimension(nrp) :: rhokn, dpdro, dpdfi
      real*8, dimension(0:nrp) :: P_rho
      real*8, dimension(nrsp4) :: rrk, cck
      real*8, dimension(nrsp6) :: wrk
      real*8 :: cwk(4)

      include 'compol.inc'
      common /com_pres/ P_rho, dpdro, dPdFi, romin

      do i=1, iplas-1
         rhocn(i+1) = (i-0.5d0)/(iplas-1.d0)
      enddo
   
      rhocn(iplas+1) = 1.d0
      rhocn(1)       = 0.d0
    
      do i=1, iplas
         rhokn(i) = (i-1.d0)/(iplas-1.d0)
      enddo

! Extrapolation q to magn. axis
 
      rh0 = 0.d0
      rh1 = 0.5d0     
      rh2 = 1.5d0     
      rh3 = 2.5d0     

      qx1 = q(1)
      qx2 = q(2)
      qx3 = q(3)

      call EXTRP2(rh0, qx0, rh1, rh2, rh3, qx1, qx2, qx3)

      dpdpsi(1) = -amu0*qx0*dpdfi(1)/flucfm
      funw(1) = qx0

      do i=1, iplas
         funw(i+1) = q(i)
      enddo

      n3spl = iplas + 1
      CALL E01BAF(n3spl, rhocn, funw, RRK, CCK, 
     &            n3spl+4, WRK, 6*n3spl+16, IFAIL)

      do i=2, iplas
         zrho = rhokn(i)
         CALL E02BCF(n3spl+4, RRK, CCK, zrho, 0, CWk, IFAIL)
         qkn = CWk(1)
         dpdpsi(i) = -amu0*qkn*dpdfi(i)/flucfm
      enddo

      return
      end subroutine pres_d_psi

!-------------------------------------------------------------
      subroutine get_eq(equil_out, psi_boundary, press0, btor)

      use imas_ids, only: type_equilibrium
      use sp_parameters

      implicit none

      real*8, intent(in) :: psi_boundary, press0, btor
      type(type_equilibrium) equil_out

      integer :: i, j
      real*8 :: pscale, fscale, um, up, dpsi
      real*8, dimension(nrp) :: C_sig, T_el, C_bts, C_driv, 
     &   Bj_av, curfi_av
      real*8, allocatable :: wrk(:)
      real*8, pointer :: rho_p(:)

      include 'compol.inc'

      common /com_sigcd/ C_sig, T_el, C_bts, C_driv
      common /com_jb/ BJ_av, curfi_av
     
      pscale = 2*TWOPI*1.d-7 * TWOPI
      fscale = TWOPI

      if(.NOT. associated(equil_out%profiles_1d%psi)) then
         allocate(equil_out%profiles_1d%psi(iplas))
         allocate(equil_out%profiles_1d%pressure(iplas))
         allocate(equil_out%profiles_1d%phi(iplas))
         allocate(equil_out%profiles_1d%pprime(iplas))
         allocate(equil_out%profiles_1d%ffprime(iplas))
         allocate(equil_out%profiles_1d%F_dia(iplas))
         allocate(equil_out%profiles_1d%q(iplas))
      endif
    
      if(.NOT. associated(equil_out%coord_sys%position%r)) then
         allocate(equil_out%coord_sys%position%r(iplas, nt1))
         allocate(equil_out%coord_sys%position%z(iplas, nt1))    
         allocate(equil_out%coord_sys%position%teta2d(nt1))    
         allocate(equil_out%coord_sys%position%rmin(iplas, nt1))    
         allocate(equil_out%coord_sys%position%psirz(iplas, nt1))    
      endif

      do i=1, iplas

! psim, psip for fixed boundary
         equil_out%profiles_1d%psi(i) = 
     &         ((psim-psip)*psia(i) + psi_boundary)*TWOPI
     
         equil_out%profiles_1d%phi(i) = flx_fi(i)
         equil_out%profiles_1d%pprime(i) = -dpdpsi(i)/pscale
         equil_out%profiles_1d%ffprime(i) = -dfdpsi(i)/fscale
         do j=1, nt1
            equil_out%coord_sys%position%r(i, j) = r(i, j+1)
            equil_out%coord_sys%position%z(i, j) = z(i, j+1)
            equil_out%coord_sys%position%teta2d(j) = teta(j+1)
            equil_out%coord_sys%position%rmin(i, j) = ro(i, j+1)
            equil_out%coord_sys%position%psirz(i, j) = -psi(i, j+1)
         enddo
      enddo

      equil_out%global_param%toroid_field%r0 = r0ax
      equil_out%global_param%toroid_field%b0 = b0ax
      equil_out%global_param%i_plasma = tokp*1e6
   
      if(.NOT. associated(equil_out%eqgeometry%boundary%r)) then
         allocate(equil_out%eqgeometry%boundary%r(nt1-1))
         allocate(equil_out%eqgeometry%boundary%z(nt1-1))
      endif
      equil_out%eqgeometry%boundary%npoints = nt1-1    
      equil_out%eqgeometry%boundary%r = 
     &   equil_out%coord_sys%position%r(iplas, 1:nt1-1)
      equil_out%eqgeometry%boundary%z = 
     &   equil_out%coord_sys%position%z(iplas, 1:nt1-1)

 !pressure     

      allocate(wrk(0:iplas+1))
    
      wrk(iplas) = 0.d0
      equil_out%profiles_1d%pressure(iplas) =
     &   wrk(iplas)/2.d0/TWOPI*1.d7
      do i=iplas-1, 1, -1 
         dpsi = (psim-psip)*(psia(i+1)-psia(i))
         wrk(i) = wrk(i+1)-0.5d0*(dpdpsi(i+1) + dpdpsi(i))*dpsi
         equil_out%profiles_1d%pressure(i) = 
     &   wrk(i)/2.d0/TWOPI*1.d7 + press0
      enddo
      deallocate(wrk)

      if(.NOT. associated(equil_out%profiles_1d%rho_tor)) then
         allocate( equil_out%profiles_1d%rho_tor(iplas) )
         allocate( equil_out%profiles_1d%jparallel(iplas) )
         allocate( equil_out%profiles_1d%sigmapar%value(iplas) )
         allocate( equil_out%profiles_1d%jni%value(iplas) )
         allocate( equil_out%profiles_1d%te%value(iplas) )
      endif
      do i=1, iplas
         equil_out%profiles_1d%rho_tor(i) = 
     &        sqrt( flx_fi(i)/(0.5*TWOPI*b0ax) )
         equil_out%profiles_1d%jparallel(i) = BJ_av(i)/b0ax/(.2d0*TWOPI)
         equil_out%profiles_1d%sigmapar%value(i) = C_sig(i)
         equil_out%profiles_1d%jni%value(i) = C_bts(i)/b0ax
         equil_out%profiles_1d%te%value(i) = T_el(i)*1d-3
      enddo

      rho_p=>equil_out%profiles_1d%rho_tor

      allocate(wrk(iplas))
    
      call cell2node(q, wrk, rho_p, iplas)
      do i=1, iplas
         equil_out%profiles_1d%q(i) = wrk(i)/TWOPI
      enddo
      call cell2node(f, wrk, rho_p, iplas)
      do i=1, iplas
         equil_out%profiles_1d%F_dia(i) = -wrk(i)
      enddo
   
      deallocate(wrk)

      return
      end subroutine get_eq

!-------------------------------------------------------
      subroutine cell2node(arr_c, arr_n, rho, iplas)

      implicit none

      integer, intent(in) :: iplas
      real*8, intent(in)  :: arr_c(*), rho(*)
      real*8, intent(out) :: arr_n(*)

      integer :: i, n3spl, IFAIL
      real*8 :: x0, x1, x2, x3, u0, u1, u2, u3, zspl
      real*8, allocatable, dimension(:) :: x, x05, RRK, CCK, WRK, CWK

      allocate(RRK(iplas+4), CCK(iplas+4), WRK((iplas+4)*6))
      allocate(CWK(4))

      allocate( x(iplas) )
      allocate( x05(iplas-1) )

      do i=1, iplas
         x(i) = rho(i)/rho(iplas)
      enddo

      do i=1, iplas-1
         x05(i) = 0.5d0*(x(i) + x(i+1))
      enddo

      x0 = x(1)
      x1 = x05(1)
      x2 = x05(2)
      x3 = x05(3)

      u1 = arr_c(1)
      u2 = arr_c(2)
      u3 = arr_c(3)

      call EXTRP2(x0, u0, X1, X2, X3, u1, u2, u3)
 
      arr_n(1) = u0

      x0 = x(iplas)
      x1 = x05(iplas-1)
      x2 = x05(iplas-2)
      x3 = x05(iplas-3)
 
      u1 = arr_c(iplas-1)
      u2 = arr_c(iplas-2)
      u3 = arr_c(iplas-3)
 
      call EXTRP2(x0, u0, X1, X2, X3, u1, u2, u3)
    
      arr_n(iplas) = u0
 
      n3spl = iplas-1

      CALL E01BAF(n3spl, x05, arr_c, RRK, CCK, 
     &            n3spl+4, WRK, 6*n3spl+16, IFAIL)

      do i=2, iplas-1
         zspl = x(i)
         CALL E02BCF(n3spl+4, RRK, CCK, zspl, 0, CWk, IFAIL)
         arr_n(i) = CWk(1)
      enddo

      deallocate(RRK, CCK, WRK)
      deallocate(CWK)
      deallocate( x )
      deallocate( x05 )

      return
      end subroutine cell2node

!-------------------------------------------------------------
      subroutine get_umup(umt, upt)

      use comblc, only: um, up

      implicit none

      real*8, intent(out) :: umt, upt

      umt = um
      upt = up

      return
      end subroutine get_umup

!-------------------------------------------------------------
      subroutine f_get_umup(umt, upt)

      use sp_parameters

      implicit none

      real*8, intent(out) :: umt, upt

      include 'compol.inc'

      umt = psim
      upt = psip

      return
      end subroutine f_get_umup

!-------------------------------------------------------------
      subroutine vol2d_eq(equil_out, i1, i2, fcell2d, fvav)

      use imas_ids, only: type_equilibrium

      implicit none

      integer, intent(in) :: i1, i2
      type(type_equilibrium) equil_out
      real*8 :: fcell2d(i1, i2)
      real*8 :: fvav(i1+1)
    
      real*8, allocatable :: fvavs(:), rhowr(:), rhowrh(:)
      integer iplas, nt1
      real*8 r1, r2, r3, r4, r0, z1, z2, z3, z4, sss
      real*8 dvoli, fvavi
      real*8, allocatable :: RRK(:), CCK(:), WRK(:)
      real*8 CWK(4)
      real*8 zspl
      integer n3spl, IFAIL

      integer i, j, jerr
      integer ip0, ip1, ip2, ip3
      real*8 x0, x1, x2, x3
      real*8 funsq

      iplas = size(equil_out%coord_sys%position%r, 1)
      nt1 = size(equil_out%coord_sys%position%r, 2)

      allocate(fvavs(iplas+1), rhowr(iplas), rhowrh(iplas+1))

      do i=1, iplas-1
        dvoli = 0
        fvavi = 0
        do j=1, nt1-1
          r1 = equil_out%coord_sys%position%r(i, j)
          r2 = equil_out%coord_sys%position%r(i+1, j)
          r3 = equil_out%coord_sys%position%r(i+1, j+1)
          r4 = equil_out%coord_sys%position%r(i, j+1)
          r0 = (r1 + r2 + r3 + r4)*0.25d0
          z1 = equil_out%coord_sys%position%z(i, j)
          z2 = equil_out%coord_sys%position%z(i+1, j)
          z3 = equil_out%coord_sys%position%z(i+1, j+1)
          z4 = equil_out%coord_sys%position%z(i, j+1)
          sss = funsq(r1, r2, r3, r4, z1, z2, z3, z4)
          dvoli = dvoli + r0*sss
          fvavi = fvavi + fcell2d(i, j)*r0*sss
        enddo
        fvavs(i+1) = fvavi/dvoli
      enddo
      do i=1, iplas
        rhowr(i) = 
     &   sqrt(equil_out%profiles_1d%phi(i)
     &   /equil_out%profiles_1d%phi(iplas))
      enddo
      do i=2, iplas
        rhowrh(i) = (rhowr(i) + rhowr(i-1))*0.5d0
      enddo
      rhowrh(iplas+1) = rhowr(iplas)
      rhowrh(1)       = 0.d0
! Extrapolation to the axis
      ip0 = 1
      ip1 = 2
      ip2 = 3
      ip3 = 4
      x3 = rhowrh(ip3)
      x2 = rhowrh(ip2)
      x1 = rhowrh(ip1)
      x0 = rhowrh(ip0)
      call b_extrp(x0, x1, x2, x3, ip1, ip2, ip3, ip0, fvavs, jerr)
! Extrapolation to the boundary
      ip0 = iplas + 1
      ip1 = iplas-2
      ip2 = iplas-1
      ip3 = iplas
      x1 = rhowrh(ip3)
      x2 = rhowrh(ip2)
      x3 = rhowrh(ip1)
      x0 = rhowrh(ip0)
      call b_extrp(x0, x1, x2, x3, ip3, ip2, ip1, ip0, fvavs, jerr)  
! Cubic splines
      n3spl = iplas + 1
      allocate(RRK(n3spl+4), CCK(n3spl+4), WRK(6*n3spl+16))
      CALL E01BAF(n3spl, rhowrh, fvavs, RRK, CCK, 
     &            n3spl+4, WRK, 6*n3spl+16, IFAIL)
      do i=1, iplas
        zspl = rhowr(i)
        CALL E02BCF(n3spl+4, RRK, CCK, zspl, 0, CWK, IFAIL)
        fvav(i) = CWK(1)
      enddo
      fvav(iplas) = fvavs(iplas+1)

      return
      end subroutine vol2d_eq

!-------------------------------------------------------------
      subroutine put_eq(equil_out, parameters_spider)

      use imas_ids, only: type_equilibrium
      use parameters_a2spider, only: type_parameters

      use sp_parameters

      implicit none

      type(type_equilibrium) equil_out
      type(type_parameters) parameters_spider

      integer, parameter :: nfour_maxx=5, nrpl=nrp+1

      integer :: i, j, i1, i2, ip0, ip1, ip2, ip3, ji, jerr, 
     &   jindx, n3spl
      real*8 :: dsqi, Zcurr, totcurr, Rcurr, dumm1, dumm2, dum1,
     &   psplexs, dpsplexs, qedge, dum_area, dum_perim, ggreen, 
     &   ddllt, psplex, pressvoli, theta_poloidal, x0, x1, x2, x3,
     &   r0, r1, r2, r3, r4, z0, z1, z2, z3, z4, u1, u2, u3, u4,
     &   btor, rtor, voli, dvoli, cur0, cur1, cur2, cur3, cur4,
     &   bmax, bp2, bp2av, bmx, bmn, bavr, bavr2, bavrd2, bnor,
     &   bm_pol, bp_pol, bpol2v, bpol2, b_fi2, b_tot, b_tot2,
     &   svol, b_ij, sss, dvdr, dvdz, gr2, 
     &   avr_v, av_r2, av_gr1, av_grr2, avr2,
     &   yrr, yrmax, yrmin, yzmax, yzmin, yrzmax, yrzmin
      real*8, dimension(3) :: xxxx1, yyyy1, pppp1
      real*8, dimension(nrp) :: upls

! 1d allocatables
      real*8, allocatable, dimension(:) :: rhowr, rhowrh, rhos,
     &   g11, g22, g33, g2int, g41, gradro, droda, dpsidv, rbp_b2, 
     &   bplfs, dvdpsi, surf_s, surf_x, g11s, g22s, g33s, g2ints, 
     &   g41s, gradros, drodas, dpsidvs, bplfss, dvdpsis, sa, volum, 
     &   dvdro, b_maxt, b_mint, b_db02, b_db0, b_0db2, bmaxt, bmint, 
     &   bdb02, bdb0, b0db2, yFOFB, areats, areat, perims, perim,
     &   ya, yra, yri, yshif, yshiv, yelon, ytria_u, ytria_l,
     &   traps, dvoliz

! 2d allocatables
      real*8, allocatable, dimension(:, :) :: 
     &   acosB2a , asinB2a , acosBlnBa , asinBlnBa, 
     &   acosB2as, asinB2as, acosBlnBas, asinBlnBas, 
     &   arr2, Btot, fcell2d1, fcell2d2

      real*8, allocatable, dimension(:, :, :) :: fcell2d21, fcell2d22,
     &   fcell2d23, fcell2d24

      include 'compol.inc'

      common /com_volt/ upls


!2d in cells
      allocate(equil_out%coord_sys%gradvcell(iplas-1, nt1-1))
      allocate(equil_out%coord_sys%bpcell(iplas-1, nt1-1))
      allocate(equil_out%coord_sys%bcell(iplas-1, nt1-1))
      allocate(equil_out%coord_sys%rcell(iplas-1, nt1-1))
! for new computations
      allocate(fcell2d1(iplas-1, nt1-1))
      allocate(fcell2d2(iplas-1, nt1-1))
      allocate(fcell2d21(iplas-1, nt1-1, nfour_maxx))
      allocate(fcell2d22(iplas-1, nt1-1, nfour_maxx))
      allocate(fcell2d23(iplas-1, nt1-1, nfour_maxx))
      allocate(fcell2d24(iplas-1, nt1-1, nfour_maxx))

      allocate( g11(iplas), g22(iplas), g33(iplas) )
      allocate( g41(iplas) )
      allocate( rbp_b2(iplas), bplfs(iplas))
      allocate( g2int(iplas) , dpsidv(iplas))
      allocate( gradro(iplas), droda(iplas) )
      allocate( volum(iplas) )
      allocate( dvdpsi(iplas) )
      allocate( dvdro(iplas) )
      allocate( surf_x(iplas) )
      allocate( surf_s(iplas+1) )
      allocate( areat(iplas) )
      allocate( areats(iplas+1) )
      allocate( perim(iplas) )
      allocate( perims(iplas+1) )

      allocate( acosB2a(iplas, nfour_maxx) )
      allocate( asinB2a(iplas, nfour_maxx) )
      allocate( acosBlnBa(iplas, nfour_maxx) )
      allocate( asinBlnBa(iplas, nfour_maxx) )

      allocate( acosB2as(iplas+1, nfour_maxx) )
      allocate( asinB2as(iplas+1, nfour_maxx) )
      allocate( acosBlnBas(iplas+1, nfour_maxx) )
      allocate( asinBlnBas(iplas+1, nfour_maxx) )

      allocate( g11s(iplas+1), g22s(iplas+1), g33s(iplas+1) )
      allocate( g41s(iplas+1) )
      allocate( g2ints(iplas+1), dpsidvs(iplas+1) )
      allocate( gradros(iplas+1), drodas(iplas+1) )
      allocate( bplfss(iplas+1))
      allocate( dvoliz(iplas+1))
      allocate( dvdpsis(iplas+1) )
      allocate( sa(iplas+1) )

      allocate( arr2(iplas, nt) )
    
      allocate( rhowr(iplas), rhowrh(iplas+1) )
      allocate( rhos(iplas) )
      allocate ( Btot(iplas, nt) )

      allocate (b_maxt(iplas+1), b_mint(iplas+1), 
     &          b_db02(iplas+1), b_db0(iplas+1), b_0db2(iplas+1))
   
      allocate (bmaxt(iplas), bmint(iplas), 
     &          bdb02(iplas), bdb0(iplas), b0db2(iplas))
     
      allocate (traps(iplas+1))
      allocate (yFOFB(iplas))
    
      btor = b0ax
      rtor = r0ax

      volum(1) = 0.d0
      pressvoli = 0.d0
      do i=2, iplas
         voli = 0.d0
         do j=2, nt1
             voli = voli + vol(i-1, j)
         enddo
         volum(i) = voli + volum(i-1)
         pressvoli = pressvoli +
     &      equil_out%profiles_1d%pressure(i-1)*voli
      enddo

      do i=1, iplas
         rhowr(i) = sqrt(flx_fi(i)/flx_fi(iplas))
         rhos(i)  = sqrt(flx_fi(i)/(pi*btor))
      enddo

!note that SPIDER natural grid is equispaced rhos=sqrt(Phi/(pi*btor)) !

      do i=2, iplas
         rhowrh(i) = (rhowr(i) + rhowr(i-1))*0.5d0
      enddo
      rhowrh(iplas+1) = rhowr(iplas)
      rhowrh(1) = 0.d0

!  A=GREENI(R, Z, RP, ZP)
      psplexs = 0.
      psplex = 0.
      dpsplexs = 0.
      Zcurr = 0.
      Rcurr = 0.
      totcurr = 0.

      do i=1, iplas-1
         dsqi = 0.d0
         dvoli = volum(i+1)-volum(i)
         u1 = volum(i)
         u2 = volum(i+1)
         u3 = volum(i+1)
         u4 = volum(i)
         dvdro(i+1) = twopi*dvoli/(rhos(i+1)-rhos(i))

         do j=2, nt1
            cur1 = cur(i, j)
            cur2 = cur(i+1, j)
            cur3 = cur(i+1, j+1)
            cur4 = cur(i, j+1)
            cur0 = (cur1 + cur2 + cur3 + cur4)*0.25d0

            r1 = r(i, j)
            r2 = r(i+1, j)
            r3 = r(i+1, j+1)
            r4 = r(i, j+1)
            r0 = (r1 + r2 + r3 + r4)*0.25d0

            z1 = z(i, j)
            z2 = z(i+1, j)
            z3 = z(i+1, j+1)
            z4 = z(i, j+1)
            z0 = (z1 + z2 + z3 + z4)*0.25d0
    
            Zcurr  = Zcurr + z0*cur0*s(i, j)
            Rcurr  = Rcurr + r0*cur0*s(i, j)
            totcurr = totcurr + cur0*s(i, j)

            sss = s(i, j)
            dvdr =  0.5d0*((u1 + u2)*(z2-z1) + (u2 + u3)*(z3-z2) + 
     &                     (u3 + u4)*(z4-z3) + (u4 + u1)*(z1-z4))/sss

            dvdz = -0.5d0*((u1 + u2)*(r2-r1) + (u2 + u3)*(r3-r2) + 
     &                     (u3 + u4)*(r4-r3) + (u4 + u1)*(r1-r4))/sss

            gr2 = dvdr**2 + dvdz**2
            arr2(i, j) = gr2
            equil_out%coord_sys%gradvcell(i, j-1) = sqrt(gr2)*twopi
         enddo
      enddo

      Zcurr = Zcurr/totcurr
      Rcurr = Rcurr/totcurr

      call avr2_c(arr2, iplas, nt, g11s(2))
      areats(1) = 0.
      perims(1) = 0.
      g2ints = 0.

      do i=1, iplas-1

         dsqi = 0.d0
         dvoli = 0.
         dvoliz(i) = dvoli
       
         av_r2 = 0.d0
         avr2 = 0.d0
         av_gr1 = 0.d0
         av_grr2 = 0.d0
         dum_area = 0.
         dum_perim = 0.

         do j=2, nt1
            dsqi = dsqi + sr(i+1, j)
            dvoli = dvoli + vol(i, j)

            r1 = r(i, j)
            r2 = r(i+1, j)
            r3 = r(i+1, j+1)
            r4 = r(i, j+1)
            r0 = (r1 + r2 + r3 + r4)*0.25d0

            gr2 = arr2(i, j)      
            av_gr1  = av_gr1  + vol(i, j)*dsqrt(gr2)
            av_grr2 = av_grr2 + vol(i, j)*gr2/r0**2
            av_r2 = av_r2 + vol(i, j)/r0**2
            avr2 = avr2 + r0**2.*vol(i, j)

            equil_out%coord_sys%rcell(i, j-1) = r0
            dum_area = dum_area + s(i, j)   
            dum_perim = dum_perim + dlt(i+1, j)   
         enddo

         sa(i+1) = dsqi*twopi
         g22s(i+1) = av_grr2/dvoli
         g2ints(i+1) = g2ints(i) + av_grr2
         g33s(i+1) = av_r2/dvoli
         g41s(i+1) = avr2/dvoli
         gradros(i+1) = av_gr1/dvoli
         drodas(i+1) = dvoli/(rhowr(i+1)-rhowr(i))
         surf_s(i+1) = sa(i+1)
         areats(i+1) = areats(i) + dum_area
         perims(i+1) = dum_perim
      enddo

!  magnetic field

      bp2av = 0.d0
      bp2 = 0.d0
      sa(1) = 0.

      do i=1, iplas-1 

         dsqi = 0.d0
          dvoli = 0.d0
          bmx = 0.d0
          bmn = 1.d9*btor
          bavr = 0.d0
          bavr2 = 0.d0
          bavrd2 = 0.d0
          theta_poloidal = 0.d0

          do j=2, nt1

             dum1 = dvoli
             dvoli = dvoli + vol(i, j)
             r1 = r(i, j)
             r2 = r(i+1, j)
             r3 = r(i+1, j+1)
             r4 = r(i, j+1)
             r0 = (r1 + r2 + r3 + r4)*0.25d0

             bm_pol = psim*(psia(i+1)-psia(i))/st(i, j)
             bp_pol = psim*(psia(i+1)-psia(i))/st(i, j+1)

             if(i /= 1) then
                bpol2v = bm_pol**2*( vol1(i, j)/sin1(i, j)   +
     &                               vol2(i, j)/sin2(i, j) ) + 
     &                   bp_pol**2*( vol3(i, j)/sin3(i, j)   +
     &                               vol4(i, j)/sin4(i, j) )
             else
                bpol2v = bm_pol**2*( vol2(i, j)/sin2(i, j) ) + 
     &                   bp_pol**2*( vol3(i, j)/sin3(i, j) ) 
            endif

            bpol2 = bpol2v /vol(i, j)
            b_fi2 = (f(i)/r0)**2
            equil_out%coord_sys%bpcell(i, j-1) = sqrt(bpol2)
            b_tot2 = bpol2 + b_fi2
            b_tot = dsqrt(b_tot2)
            equil_out%coord_sys%bcell(i, j-1) = b_tot
            fcell2d1(i, j-1) = bpol2/b_tot2*r0**2.0
            dsqi = dsqi + dlt(i+1, j)/sqrt(bpol2)
            Btot(i, j) = b_tot  ! aai 06/02/10
            bmx = dmax1(b_tot, bmx)
            bmn = dmin1(b_tot, bmn)
            bavr2 = bavr2 + b_tot2*vol(i, j)
            bavr  = bavr  + b_tot *vol(i, j)
            bavrd2 = bavrd2 + vol(i, j)/(b_tot2 + 1.d-8)
            bp2av = bp2av + bpol2v
            bp2   = bp2   + bpol2v
         enddo
         dvoliz(i) = dvoli   
         dvdpsis(i+1) = dsqi
         dpsidvs(i+1) = 1./dvdpsis(i+1)
         b_maxt(i+1) = bmx
         b_mint(i+1) = bmn
         b_db02(i+1) = bavr2/(dvoli*btor**2)
         b_0db2(i+1) = (bavrd2*btor**2)/(dvoli)
         b_db0(i+1) = bavr/(dvoli*btor)

      enddo 

      equil_out%global_param%li3 = 4.d0*pi*bp2/
     &   (rtor*(.4d0*pi*tokp)**2)
      equil_out%global_param%betpol = 
     &   (pressvoli/bp2)*2.*(2.*TWOPI*1.E-7)
      equil_out%global_param%wkin = TWOPI*pressvoli
      equil_out%global_param%bpkin = TWOPI*bp2/(2.*(2.*TWOPI*1.E-7))

      do i=1, iplas-1
         Bmax = 0.d0
         Svol = 0.d0   
         do j=2, nt-1
            B_ij = Btot(i, j)
            Bmax = dmax1(Bmax, B_ij)
            Svol = Svol + vol(i, j)
         enddo

         avr_v = 0.d0
         do j=2, nt-1
            Bnor = Btot(i, j)/Bmax
            avr_v = avr_v +  (b0ax/Btot(i, j))**2 *
     &         ( 1.d0-dsqrt(1.d0-Bnor)*(1.d0 + .5d0*Bnor) )*vol(i, j)
         enddo
    
         traps(i+1) = avr_v/Svol
     
      enddo

! stuff for psi as b.c.

      qedge = equil_out%profiles_1d%q(iplas)
      if (parameters_spider%k_grid == 0) then
         call psib_pla(psplexs)    ! spider G integral
         psplexs = psplexs/(btor*rhos(iplas)**2)*qedge  ! spider G integral
      endif

      if (parameters_spider%k_grid == 1) then
         call psib_pla(psplexs)    ! spider G integral
         psplexs = psplexs/(btor*rhos(iplas)**2)*qedge  ! spider G integral
      endif

! cosinus sinus fourier modes of theta_gen

      ip0 = max(2, minloc(abs(teta(1:nt1)), 1))
      do i=1, iplas-1 

         dsqi = 0.d0
         dvoli = 0.d0
         theta_poloidal = 0.d0

         do j=ip0, nt1
            dvoli = dvoli + equil_out%coord_sys%bcell(i, j-1)*vol(i, j)
            theta_poloidal = TWOPI*dvoli/dvoliz(i)*1./(btor*b_db0(i+1))

            do jindx=1, nfour_maxx
               fcell2d21(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)**2.*
     &            cos(jindx*theta_poloidal)
               fcell2d22(i, j-1, jindx) =
     &            equil_out%coord_sys%bcell(i, j-1)*
     &            log(equil_out%coord_sys%bcell(i, j-1))*
     &            cos(jindx*theta_poloidal)
               fcell2d23(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)**2.*
     &            sin(jindx*theta_poloidal)
               fcell2d24(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)*
     &            log(equil_out%coord_sys%bcell(i, j-1))*
     &            sin(jindx*theta_poloidal)
            enddo
         enddo

         do j=2, ip0-1
            dvoli = dvoli + equil_out%coord_sys%bcell(i, j-1)*vol(i, j)
            theta_poloidal = TWOPI*dvoli/dvoliz(i)*1./(btor*b_db0(i+1))

            do jindx=1, nfour_maxx
               fcell2d21(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)**2.*
     &            cos(jindx*theta_poloidal)
               fcell2d22(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)*
     &            log(equil_out%coord_sys%bcell(i, j-1))*
     &            cos(jindx*theta_poloidal)
               fcell2d23(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)**2.*
     &            sin(jindx*theta_poloidal)
               fcell2d24(i, j-1, jindx) = 
     &            equil_out%coord_sys%bcell(i, j-1)*
     &            log(equil_out%coord_sys%bcell(i, j-1))*
     &            sin(jindx*theta_poloidal)
            enddo
         enddo      
      enddo 

! COmpute <R^2 Bp^2 / B^2>
      call vol2d_eq(equil_out, iplas-1, nt1-1, fcell2d1, rbp_b2)

      do jindx=1, nfour_maxx
         call vol2d_eq(equil_out, iplas-1, nt1-1, 
     &      fcell2d21(:, :, jindx), acosB2a(:, jindx))
         call vol2d_eq(equil_out, iplas-1, nt1-1, 
     &      fcell2d22(:, :, jindx), acosBlnBa(:, jindx))
         call vol2d_eq(equil_out, iplas-1, nt1-1, 
     &      fcell2d23(:, :, jindx), asinB2a(:, jindx))
         call vol2d_eq(equil_out, iplas-1, nt1-1, 
     &      fcell2d24(:, :, jindx), asinBlnBa(:, jindx))
      enddo

!Compute bp at low field side
      dumm1 = 100.
      do jindx=2, nt1
         dumm2 = abs(z(iplas, jindx)-z(1, 2))
         if (dumm2 < dumm1 .and. r(iplas, jindx) > r(1, 2)) then
            j = jindx
            dumm1 = dumm2
         endif
      enddo

      j = minloc(abs(teta(2:nt1)), 1)
      do i=1, iplas-1 
         bplfss(i+1) = equil_out%coord_sys%bpcell(i, j)
      enddo
! Extrapolation to the axis

      ip0 = 1
      ip1 = 2
      ip2 = 3
      ip3 = 4

      x3 = rhowrh(ip3)
      x2 = rhowrh(ip2)
      x1 = rhowrh(ip1)
      x0 = rhowrh(ip0)
    
      g11s(1) = 0.d0
      g22s(1) = 0.d0
      g2ints(1) = 0.d0
      g33s(1) = 1.d0/rm**2
      g41s(1) = rm**2

      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, dpsidvs, jerr)
      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, gradros, jerr)

      drodas(1) = 0.d0
      areats(1) = 0.d0
      perims(1) = 0.d0
      bplfss(1) = 0.d0
      surf_s(1) = 0.d0

      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, b_db02, jerr)
      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, b_0db2, jerr)
      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, b_db0, jerr)
      call b_extrp(x0, X1, X2, X3, ip1, ip2, ip3, ip0, traps, jerr)

      b_maxt(1) = b_db0(1)*btor
      b_mint(1) = b_maxt(1)

! Extrapolation to the boundary

      ip0 = iplas+1
      ip1 = iplas-2
      ip2 = iplas-1
      ip3 = iplas

      x1 = rhowrh(ip3)
      x2 = rhowrh(ip2)
      x3 = rhowrh(ip1)
      x0 = rhowrh(ip0)

      call EXTRAP_EFSP(rhowrh(1:iplas), 
     &      g22s(1:iplas), x0, 
     &      iplas, g22s(ip0), 1, iplas)   

      call EXTRAP_EFSP(rhowrh(1:iplas), 
     &      dpsidvs(1:iplas), x0, 
     &      iplas, dpsidvs(ip0), 2, iplas)   

      psplex = psplexs

      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, g11s   , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, g2ints , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, g33s   , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, g41s   , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, gradros, jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, drodas , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, areats , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, perims , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, bplfss , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, surf_s, jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, b_db02 , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, b_0db2 , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, b_db0  , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, traps  , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, b_maxt , jerr)
      call b_extrp(x0, X1, X2, X3, ip3, ip2, ip1, ip0, b_mint , jerr)

! Extrapolation to the boundary

      n3spl = iplas+1

      call spline_prof(n3spl,iplas,rhowrh,g11s   ,rhowr,g11(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,g22s   ,rhowr,g22(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,g33s   ,rhowr,g33(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,g41s   ,rhowr,g41(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,dpsidvs,rhowr,dpsidv(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,g2ints ,rhowr,g2int(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,gradros,rhowr,gradro(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,drodas ,rhowr,droda(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,bplfss ,rhowr,bplfs(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,b_db02 ,rhowr,bdb02(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,b_0db2 ,rhowr,b0db2(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,b_db0  ,rhowr,bdb0(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,traps  ,rhowr,yFOFB(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,surf_s ,rhowr,surf_x(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,areats ,rhowr,areat(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,perims ,rhowr,perim(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,b_mint ,rhowr,bmint(1:iplas))
      call spline_prof(n3spl,iplas,rhowrh,b_maxt ,rhowr,bmaxt(1:iplas))

      if(.NOT. associated(equil_out%profiles_1d%g1)) then
         allocate(equil_out%profiles_1d%gm1(iplas))
         allocate(equil_out%profiles_1d%gm4(iplas))
         allocate(equil_out%profiles_1d%gm5(iplas))
         allocate(equil_out%profiles_1d%gm41(iplas))
         allocate(equil_out%profiles_1d%rbp_b2(iplas))
         allocate(equil_out%profiles_1d%bplfs(iplas))
         allocate(equil_out%profiles_1d%g1(iplas))
         allocate(equil_out%profiles_1d%g2(iplas))
         allocate(equil_out%profiles_1d%g2int(iplas))
         allocate(equil_out%profiles_1d%fofb(iplas))
         allocate(equil_out%profiles_1d%areat(iplas))
         allocate(equil_out%profiles_1d%perim(iplas))
         allocate(equil_out%profiles_1d%ggradro(iplas))
         allocate(equil_out%profiles_1d%gdroda(iplas))
         allocate(equil_out%profiles_1d%bmaxt(iplas))
         allocate(equil_out%profiles_1d%bmint(iplas))
         allocate(equil_out%profiles_1d%bdb0(iplas))
         allocate(equil_out%profiles_1d%dPSIdV(iplas))
         allocate(equil_out%profiles_1d%surface(iplas))

         allocate(equil_out%profiles_1d%acosB2a(iplas, nfour_maxx))
         allocate(equil_out%profiles_1d%asinB2a(iplas, nfour_maxx))
         allocate(equil_out%profiles_1d%acosBlnBa(iplas, nfour_maxx))
         allocate(equil_out%profiles_1d%asinBlnBa(iplas, nfour_maxx))
      endif

      equil_out%global_param%Zcurr = Zcurr
      equil_out%global_param%Rcurr = Rcurr
      equil_out%global_param%Vloop = upls(iplas)

      equil_out%profiles_1d%dPSIdV = dpsidv
      equil_out%profiles_1d%gm1 = g33
      equil_out%profiles_1d%gm4 = bdb02 * btor**2
      equil_out%profiles_1d%gm41 = g41 
      equil_out%profiles_1d%gm5 = b0db2 / btor**2
      equil_out%profiles_1d%surface = surf_x
      equil_out%profiles_1d%rbp_b2 = rbp_b2
      equil_out%profiles_1d%rbp_b2(1) = 0.
      equil_out%profiles_1d%bplfs = bplfs
      equil_out%profiles_1d%acosB2a = acosB2a
      equil_out%profiles_1d%asinB2a = asinB2a
      equil_out%profiles_1d%acosBlnBa = acosBlnBa
      equil_out%profiles_1d%asinBlnBa = asinBlnBa
      equil_out%profiles_1d%g1 = g11*twopi**2
      equil_out%profiles_1d%g2 = g22*twopi**2
      equil_out%profiles_1d%g2int = g2int*twopi**2
      equil_out%profiles_1d%fofb = yFOFB
      equil_out%profiles_1d%areat = areat
      equil_out%profiles_1d%perim = perim
      equil_out%profiles_1d%ggradro = gradro*twopi
      equil_out%profiles_1d%gdroda = droda
      equil_out%profiles_1d%bmaxt = bmaxt
      equil_out%profiles_1d%bmint = bmint
      equil_out%profiles_1d%bdb0 = bdb0 * btor
 
!geometry
! calculation of 3M from SPIDER output geometry
      allocate( ya(iplas), yra(iplas), yri(iplas), yshif(iplas), 
     &   yshiv(iplas), yelon(iplas), ytria_u(iplas), ytria_l(iplas) )

      do ji=2, iplas
         yrmax = -99999.d0
         yrmin = 99999.d0
         yzmax = -99999.d0
         yzmin = 99999.d0

         i = minloc(z(ji, 2:nt1), 1)  
         i = i+1
         i1 = i-1
         i2 = i+1
         if (i == 1) i1 = nt1
         if (i == nt1) i2 = 2
         xxxx1(1) = r(ji, i1)
         xxxx1(2) = r(ji, i)
         xxxx1(3) = r(ji, i2)
         yyyy1(1) = z(ji, i1)
         yyyy1(2) = z(ji, i)
         yyyy1(3) = z(ji, i2)
         call polyfitcc_spid(xxxx1, yyyy1, pppp1)
         yrzmin = -pppp1(2)/(2.*pppp1(1))
         yzmin = pppp1(1)*yrzmin**2.0+pppp1(2)*yrzmin+pppp1(3)

         i = maxloc(z(ji, 2:nt1), 1)  
         i = i+1
         i1 = i-1
         i2 = i+1
         if (i == 1) i1 = nt1
         if (i == nt1) i2 = 2
         xxxx1(1) = r(ji, i1)
         xxxx1(2) = r(ji, i)
         xxxx1(3) = r(ji, i2)
         yyyy1(1) = z(ji, i1)
         yyyy1(2) = z(ji, i)
         yyyy1(3) = z(ji, i2)
         call polyfitcc_spid(xxxx1, yyyy1, pppp1)
         yrzmax = -pppp1(2)/(2.*pppp1(1))
         yzmax = pppp1(1)*yrzmax**2.0+pppp1(2)*yrzmax+pppp1(3)
         i = minloc(r(ji, 2:nt1), 1)  
         yrmin = r(ji, i+1)
         i = maxloc(r(ji, 2:nt1), 1)  
         yrmax = r(ji, i+1)
         yrr = .5d0*(yrmax+yrmin)
         ya(ji) = .5*(yrmax-yrmin)
         yra(ji)  =  yrmax
         yri(ji)  =  yrmin
         yshif(ji) = yrr-rtor
         yshiv(ji) = .5d0*(yzmin+yzmax)
         yelon(ji) = (yzmax-yzmin)/(yrmax-yrmin)
         ytria_u(ji) = (yrr-yrzmax)/ya(ji)
         ytria_l(ji) = (yrr-yrzmin)/ya(ji)
      enddo
      ya(1) = 0.d0
      yra(1) = r(1, 2)
      yri(1) = r(1, 2)
      yelon(1) = yelon(2)
      ytria_u(1) = ytria_u(2)
      ytria_l(1) = ytria_l(2)
      yshif(1) = r(1, 2)-rtor
      yshiv(1) = z(1, 2)

      if(.NOT. associated(equil_out%profiles_1d%volume)) then
         allocate(equil_out%profiles_1d%volume(iplas))
         allocate(equil_out%profiles_1d%r_inboard(iplas))
         allocate(equil_out%profiles_1d%r_outboard(iplas))
         allocate(equil_out%profiles_1d%elongation(iplas))
         allocate(equil_out%profiles_1d%tria_upper(iplas))
         allocate(equil_out%profiles_1d%tria_lower(iplas))
         allocate(equil_out%profiles_1d%shif(iplas))
         allocate(equil_out%profiles_1d%shiv(iplas))
      endif
      equil_out%profiles_1d%volume = volum*twopi
      equil_out%profiles_1d%r_inboard = yri
      equil_out%profiles_1d%r_outboard = yra
      equil_out%profiles_1d%elongation = yelon
      equil_out%profiles_1d%tria_upper = ytria_u
      equil_out%profiles_1d%tria_lower = ytria_l        
      equil_out%profiles_1d%shif = yshif
      equil_out%profiles_1d%shiv = yshiv 

!get rectangular psi
      call get2d_eq_ef(equil_out)

!get psplex
      equil_out%global_param%psplex = psplex

!Deallocation    
      deallocate( g11, g22, g33 )
      deallocate( g41 )
      deallocate( rbp_b2, bplfs)
      deallocate( g2int , dpsidv)
      deallocate( gradro, droda )
      deallocate( volum )
      deallocate( dvdpsi )
      deallocate( dvdro )
      deallocate( surf_s )
      deallocate( surf_x )
      deallocate( areat )
      deallocate( areats )
      deallocate( perim )
      deallocate( perims )

      deallocate( acosB2a )
      deallocate( asinB2a )
      deallocate( acosBlnBa )
      deallocate( asinBlnBa )

      deallocate( acosB2as )
      deallocate( asinB2as )
      deallocate( acosBlnBas )
      deallocate( asinBlnBas )

      deallocate( g11s, g22s, g33s )
      deallocate( g41s)
      deallocate( g2ints, dpsidvs )
      deallocate( gradros, drodas )
      deallocate( bplfss)
      deallocate( dvoliz)
      deallocate( dvdpsis )
      deallocate( sa )

      deallocate( arr2 )
    
      deallocate( rhowr, rhowrh )
      deallocate( rhos )
      deallocate ( Btot )

      deallocate (b_maxt, b_mint, b_db02, b_db0, b_0db2)
      deallocate (bmaxt, bmint, bdb02, bdb0, b0db2)
     
      deallocate (traps)
      deallocate (yFOFB)
      deallocate(ya, yra, yri, yshif, yshiv, 
     &  yelon, ytria_u, ytria_l)

      return
      end subroutine put_eq

!----------------------------------------------------------
      subroutine get2d_eq_ef(equil_out)
!get rectangular part

      use imas_ids, only: type_equilibrium
      use sp_parameters, only: twopi
      use comblc, only: ni1, nj1, r, z, u, um, up

      implicit none

      integer :: i, j

      type(type_equilibrium) :: equil_out

      if(.NOT. associated(equil_out%eqgeometry%rectgrid%r2d)) then
         allocate(equil_out%eqgeometry%rectgrid%r2d(ni1))
         allocate(equil_out%eqgeometry%rectgrid%z2d(nj1))
         allocate(equil_out%eqgeometry%rectgrid%psirz2d(ni1, nj1))
      endif

      equil_out%eqgeometry%rectgrid%npointsr = ni1
      equil_out%eqgeometry%rectgrid%npointsz = nj1

      do j=1, nj1
         do i=1, ni1
            equil_out%eqgeometry%rectgrid%r2d(i) = r(i)
            equil_out%eqgeometry%rectgrid%z2d(j) = z(j)
            equil_out%eqgeometry%rectgrid%psirz2d(i, j) = u(i, j)
         enddo
      enddo

      equil_out%global_param%psibound = -TWOPI*up
      equil_out%global_param%psiaxis  = -TWOPI*um

      return
      end subroutine get2d_eq_ef

!-----------------------------------------------------
      subroutine polyfitcc_spid(x, y, P)

      implicit none

      double precision x(3), y(3), P(3)
      double precision y21, y32, x21, x32, h21, h32

      y32 = y(3)-y(2)     
      y21 = y(2)-y(1)     
      x32 = x(3)-x(2)     
      x21 = x(2)-x(1)     
      h32 = x(3)+x(2)     
      h21 = x(2)+x(1)     

      P(1) = (x21*y32-x32*y21)/(x21*x32*(h32-h21))
      P(2) = y21/x21-P(1)*h21
      P(3) = y(3)-P(1)*x(3)**2.-P(2)*x(3)

      return
      end subroutine polyfitcc_spid

!------------------------------------------------------
      subroutine polyfitcc_spid_1(x, y, P)

      implicit none

      double precision x(2), y(2), P(2)
      double precision y21, x21

      y21 = y(2)-y(1)     
      x21 = x(2)-x(1)     

      P(1) = y21/x21
      P(2) = y(1)-x(1)*y21/x21

      return
      end subroutine polyfitcc_spid_1

!--------------------------------------------------------
      subroutine polyfitcc_spid_0(x, y, P)

      implicit none

      double precision x(1), y(1), P(1)

      P(1) = y(1)

      return
      end subroutine polyfitcc_spid_0

!--------------------------------------------------------
      subroutine EXTRAP_EFSP(x_input, y_input, x_extrap, 
     &  j_extrap, y_extrap, ex_order, nagrid)
! Extrapolations: assume that x is of r-type, i.e. interpolation in 0
! has zero odd derivatives

      implicit none

      integer k1, k2, k3, ex_order, j_extrap, nagrid
      integer jsign
      double precision x_input(nagrid), x_extrap
      double precision y_input(nagrid), y_extrap
      double precision P(3)
   
      if (j_extrap == nagrid)  jsign = -1
      if (j_extrap == 1)  jsign = 1
     
! Constant interpolation
      if (ex_order == 0) then
         k1 = j_extrap
         call polyfitcc_spid_0(x_input(k1), y_input(k1), P(1))
         y_extrap = P(1)
      endif

! Linear interpolation
      if (ex_order == 1) then
         if (jsign < 0) then
            k1 = j_extrap+jsign*1
            k2 = j_extrap
            call polyfitcc_spid_1(x_input(k1:k2),y_input(k1:k2),P(1:2))
            y_extrap = P(1)*x_extrap+P(2)
         endif
         if (jsign > 0) then
            k1 = j_extrap
            call polyfitcc_spid_0(x_input(k1), y_input(k1), P(1))
            y_extrap = P(1)
         endif
      endif

! Quadratic interpolation
      if (ex_order == 2) then
         if (jsign < 0) then
            k1 = j_extrap+jsign*2
            k2 = j_extrap+jsign*1
            k3 = j_extrap
            call polyfitcc_spid(x_input(k1:k3), y_input(k1:k3), P)
            y_extrap = P(1)*x_extrap**2.0+P(2)*x_extrap+P(3)
         endif
         if (jsign > 0) then
            k1 = j_extrap
            k2 = j_extrap+jsign*1
            k3 = j_extrap+jsign*2
            call polyfitcc_spid(x_input(k1:k3), y_input(k1:k3), P)
            y_extrap = P(1)*x_extrap**2.0+P(2)*x_extrap+P(3)
         endif
   
      endif

      return
      end subroutine EXTRAP_EFSP
