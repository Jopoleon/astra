      subroutine eqb( alf0, alf1, alf2, bet0, bet1, bet2, 
     &   alw0, alw1, alw2, betplx, i_betp,
     &   keyctr, nstep, platok, rax, zax, b0cen, r0cen, psax, igdf,
     &   n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, psi_bnd, psi0_bnd)

      use sp_parameters, only: nursp
      use iopath, only: path
      use keys, only: kpr, kstep
      use compol, only: nr, nr1, nr2, nt, nt1, nt2, iplas, iplas1, ngav,
     & nitbeg, nitdel, iter, itin, tokp, fvac,
     & psia, psiax, psibon, psi_eav, psibon0, psim, psip, psipla,
     &     r0ax, b0ax, tok, dfdpsi, dpdpsi, cnor, erru, r, z, rm, zm

      implicit none

      integer, parameter :: nursp4=nursp+4,nursp6=nursp4*6

      integer, intent(in) :: i_betp, keyctr, nstep, igdf, n_tht, n_psi,
     &   nurs, i_eqdsk, i_bsh
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2, 
     &   alw0, alw1, alw2, betplx, b0cen, r0cen,
     &   epsro, psi0_bnd
      real*8, intent(inout) :: psi_bnd, platok, rax, zax, psax

      integer ::i, igrdf, imov, max_max_iterb, nbsh, itrmax, nitmax,
     &   NiMax
      real*8 :: rm0, zm0, rc0, zc0, asp0, el_up, el_lw, tr_up, tr_lw,
     &   alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, alw0p, alw1p, alw2p,
     &   betpol, bettot, dummy_nbsh, erro, errpsi, errod, qax, q_giv0
      real*8 :: tabp, tabf
      character(len=80) :: fname

      common /combsh/ rm0, zm0, rc0, zc0, asp0, el_up, el_lw, 
     &                tr_up, tr_lw, nbsh
      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, bet2f,
     &                alw0p, alw1p, alw2p
      common /maxiterib/ max_max_iterb

      NiMax = 500

      if(keyctr.ge.100) return

      nt    = n_tht
      iplas = n_psi
      ngav  = keyctr
      tok   = platok
      nbsh  = i_bsh
      b0ax  = b0cen
      r0ax  = r0cen
      psiax = psax
      if(ngav.eq.-1) then
         psibon=psi_bnd
      else
         psi_eav=psi_bnd
      endif

      psibon0=psi0_bnd
      kstep = nstep

! input initial data

      rm0   = rax
      zm0   = zax 
         
      alf0p = alf0 
      alf1p = alf1 
      alf2p = alf2 

      bet0f = bet0
      bet1f = bet1
      bet2f = bet2

      alw0p = alw0 
      alw1p = alw1 
      alw2p = alw2 

      if(i_eqdsk.eq.1) then
         nbsh=1
      endif

      nr  = iplas
      nr1 = nr-1
      nt1 = nt-1
      nr2 = nr-2
      nt2 = nt-2

      iplas1 = iplas-1
      itrmax = 50
      nitmax = 5
      nitdel = 10
      nitbeg = 5
      igrdf = igdf

      if(nstep.eq.0) then
         if(nbsh.eq.0) then
            write(fname,'(a,a)') TRIM(path),'/inpol.dat'
            open(1,file=fname)
               read(1,*) dummy_nbsh
               read(1,*) rc0
               read(1,*) zc0
               read(1,*) asp0
               read(1,*) el_up
               read(1,*) el_lw
               read(1,*) tr_up
               read(1,*) tr_lw
            close(1)    ! for "inpol.dat" reading
         endif

         igrdf = igdf
         if(igdf.eq.2) igrdf=1

         if(nbsh.eq.0.or.nbsh.eq.1) then  
            call grid_0(igrdf, nstep)
            call taburs(0, 1.d0, nurs)
         elseif(nbsh.eq.-1) then
            call grid_b(igrdf,nstep)
         elseif(nbsh.eq.-2) then
            call grid_p0(igrdf,nstep)
            call taburs(0,1.d0,nurs)
         endif

         do i=1,iplas
            dpdpsi(i)=tabp(psia(i))
            dfdpsi(i)=tabf(psia(i))
         enddo

         if(ngav.eq.0) call presol(i_betp,betplx,betpol)

      else 
         if(nbsh.eq.0.or.nbsh.eq.1) then  
            call grid_1(igrdf,nstep)
         elseif(nbsh.eq.-1) then
            call grid_b1(igrdf,nstep)
         elseif(nbsh.eq.-2) then
            call grid_p1(igrdf,nstep)
         endif

      endif   !nstep.eq.0

      if(ngav.lt.0) then
         call bongri
      endif   

      if(nbsh.eq.-1) then
         call psib_ext(psi_eav)
      endif   

      cnor=1.d0

      iter  = 0
      itin  = 0
      imov  = 1

      do
         iter  = iter+1
         itin  = itin+1
         igrdf = igdf

         call metrix
         call matcof
         call matpla
         call rightg

	 if(i_betp.eq.1) then
	    if(iter.gt.4) call skbetp(betplx,betpol)
         endif
	 if(ngav.le.0 ) then
	    call qst_b
         endif

	 call grdef(igrdf)
         call solint(imov)
         call remesh(erro,errpsi,imov)

         errod = 0.5d0*erro/(dabs(z(iplas,2)-zm)+dabs(r(iplas,2)-rm))
         erru=erro
         if(errod.lt.epsro) EXIT
         if(itin.ge.NiMax) then 
            write(*,*) 'SPIDER: MAX NUMBER OF ITERATIONS IS EXCEEDED'
            write(*,*) 'ERROD=',errod
            write(*,*) 'ITER=',itin
            call err_catch_a !EFable ereror catchin in astra 
	    EXIT
         endif             

	 if (kpr.eq.-2.and.iter.ge.max_max_iterb) then
            errod=1.e-12
            EXIT
         endif

      enddo

      if(ngav.le.0) then
         call qst_b
         psiax=psim
      endif

      platok = tokp
      rax    = rm
      zax    = zm
      psax   = psim

      if(ngav.le.0) then
         psax   = psim
         psipla = psim-psip
      endif
      if(ngav.eq.-3) then
         psi_bnd=psi_eav
      endif

      if(ngav.gt.0) qax = q_giv0
      if(ngav.lt.0) call retab_L
      if(ngav.gt.0) call retab_p
      call bt_pol(betpol)
      call bt_tot(bettot)

      if(kpr.eq.1) then
         write(*,*) '*******************************************'
         write(*,*) 'Iterations have been converged for epsro =',epsro
         write(*,*) 'Achieved accuracy ................ errod =',errod
         write(*,*) 'Number of iterations ............. iter  =',iter
         write(*,*) 'Magnetic axis coordinates ........ rm    =',rm
         write(*,*) '                                   zm    =',zm
         write(*,*) 'Plasma current ................... tokp  =',tokp
         write(*,*) '                                   cnor  =',cnor
         write(*,*) 'Magnetic axis PSI value .......... psim  =',psim
         write(*,*) 'Poloidal beta value ............. betpol =',betpol
         write(*,*) 'Total    beta value ............. bettot =',bettot
         write(*,*) 'Plas.boun.poloidal current value... fvac =',fvac
         write(*,*) '*******************************************'
      endif

      call wrb
!      call grid_b_ef2(igdf,nstep) ! git useless subroutine

      return
      end subroutine eqb

!----------------------------------------------------------------
      subroutine eqb_contour( alf0, alf1, alf2,  
     &    bet0, bet1, bet2, alw0, alw1, alw2, betplx, i_betp,
     &    keyctr, nstep, platok, rax, zax, b0cen, r0cen, psax, igdf,
     &    n_tht, n_psi, epsro, nurs, i_eqdsk, i_bsh, psi_bnd, psi0_bnd)

      use sp_parameters, only: nursp
      use keys, only: kstep
      use compol, only: nr, nr1, nr2, nt, nt1, nt2, iplas, iplas1,
     & ngav, nitdel, nitbeg, iter, itin,
     & r0ax, b0ax, psiax,
     &     psia, psiax, psibon, psibon0, psi_eav, psim, psip, psipla,
     & dpdpsi, dfdpsi, cnor, tok, tokp, rm, zm

      implicit none

      integer, parameter :: nursp4=nursp+4,nursp6=nursp4*6

      integer, intent(in) :: i_betp, keyctr, nstep, igdf, n_tht, n_psi,
     &   nurs, i_eqdsk, i_bsh
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2, 
     &   alw0, alw1, alw2, betplx, b0cen, r0cen,
     &   epsro, psi_bnd, psi0_bnd
      real*8, intent(out) :: platok, rax, zax, psax

      integer ::i, igrdf, imov, max_max_iterb, nbsh, itrmax, nitmax
      real*8 :: rm0, zm0, rc0, zc0, asp0, el_up, el_lw, tr_up, tr_lw,
     &   alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, alw0p, alw1p, alw2p,
     &   betpol, bettot
      real*8 :: tabp, tabf

      common /combsh/ rm0, zm0, rc0, zc0, asp0, el_up, el_lw, 
     &                tr_up, tr_lw, nbsh
      common /comabw/ alf0p, alf1p, alf2p, bet0f, bet1f, bet2f,
     &                alw0p, alw1p, alw2p
      common /maxiterib/ max_max_iterb

      nt    = n_tht
      iplas = n_psi
      ngav  = keyctr
      tok   = platok
      nbsh  = i_bsh
      b0ax  = b0cen
      r0ax  = r0cen
      psiax = psax
      if(ngav.eq.-1) then
         psibon=psi_bnd
      else
         psi_eav=psi_bnd
      endif

      psibon0 = psi0_bnd
      kstep = nstep

      rm0   = rax
      zm0   = zax 
         
      alf0p = alf0 
      alf1p = alf1 
      alf2p = alf2 

      bet0f = bet0
      bet1f = bet1
      bet2f = bet2

      alw0p = alw0 
      alw1p = alw1 
      alw2p = alw2 

      nr     = iplas
      nr1    = nr-1
      nt1    = nt-1
      nr2    = nr-2
      nt2    = nt-2

      iplas1 = iplas-1

      itrmax = 50
      nitmax = 5
      nitdel = 10
      nitbeg = 5
      igrdf = igdf
      if(igdf.eq.2) igrdf=1

      call grid_b_ef(igrdf,nstep)

      do i=1,iplas
         dpdpsi(i)=tabp(psia(i))
         dfdpsi(i)=tabf(psia(i))
      enddo

      cnor=1.d0

      iter  = 0
      itin  = 0
      imov  = 0
      iter  = iter+1
      itin  = itin+1

      call metrix
      call rightg
      call grdef(igrdf)
      call qst_b

      psiax  = psim
      platok = tokp
      rax    = rm
      zax    = zm
      psax   = psim

      psax   = psim
      psipla = psim-psip
      call bt_pol(betpol)
      call bt_tot(bettot)
      call wrb
!      call grid_b_ef2(igrdf,nstep) ! git useless subroutine

      return
      end subroutine eqb_contour

!----------------------------------------------------------------
      subroutine retab_L

      use sp_parameters, only: nursp
      use compol, only: iplas, psia, dfdpsi, dpdpsi

      implicit none

      integer :: i, nurs, it, ic
      real*8 :: wes, zpsi, dpsi, furs_n, purs_n
      real*8, dimension(nursp) :: psit, purs, furs, wurs,
     &   ppp, fff, www, pstab

      common/comurs/ psit, purs, furs, wurs, nurs
      common/comppp/ ppp, fff, www

      wes=0.5d0
      wes=1.0d0

      do i=1,iplas
         pstab(i)=1.d0-psia(i)
      enddo

      do it=2,nurs-1
         zpsi=1.d0-psit(it)
         do i=1,iplas-1
            if(zpsi.gt.pstab(i) .AND. zpsi.le.pstab(i+1)) then
               ic=i
               EXIT
            endif
         enddo

         dpsi=pstab(ic+1)-pstab(ic)
         furs_n=( dfdpsi(ic)*(pstab(ic+1)-zpsi)+
     &            dfdpsi(ic+1)*(zpsi-pstab(ic)) )/dpsi
         furs(it)=wes*furs_n+(1.d0-wes)*furs(it)

         purs_n=( dpdpsi(ic)*(pstab(ic+1)-zpsi)+
     &            dpdpsi(ic+1)*(zpsi-pstab(ic)) )/dpsi
         purs(it)=wes*purs_n+(1.d0-wes)*purs(it)
      enddo

      furs(1)=wes*dfdpsi(iplas)+(1.d0-wes)*furs(1)
      furs(nurs)=wes*dfdpsi(iplas)+(1.d0-wes)*furs(nurs)

      purs(1)=wes*dpdpsi(iplas)+(1.d0-wes)*purs(1)
      purs(nurs)=wes*dpdpsi(iplas)+(1.d0-wes)*purs(nurs)

      dpsi=1.d0/(nurs-1.d0)
      ppp(1)=0.d0
      fff(1)=0.d0

      do i=2,nurs
         ppp(i)=ppp(i-1) +(purs(i-1)+purs(i))*dpsi*0.5d0
         fff(i)=fff(i-1) +(furs(i-1)+furs(i))*dpsi*0.5d0
      enddo

      return
      end subroutine retab_L

!----------------------------------------------------------------
      subroutine retab_p

      use sp_parameters, only: nursp
      use compol, only: iplas, psia, dpdpsi

      implicit none

      integer :: i, nurs, it, ic
      real*8 :: zpsi, dpsi
      real*8, dimension(nursp) :: psit, purs, furs, wurs,
     &   ppp, fff, www, pstab

      common/comurs/ psit, purs, furs, wurs, nurs
      common/comppp/ ppp, fff, www

      write(*,*) 'retab_p:linear recalculation p table'

      do i=1,iplas
         pstab(i)=1.d0-psia(i)
      enddo

      do it=2,nurs-1
         zpsi=1.d0-psit(it)
         do i=1,iplas-1
            if(zpsi.gt.pstab(i) .AND. zpsi.le.pstab(i+1)) then
               ic=i
               EXIT
            endif
         enddo

         dpsi=pstab(ic+1)-pstab(ic)
         purs(it)=( dpdpsi(ic)*(pstab(ic+1)-zpsi)+
     &              dpdpsi(ic+1)*(zpsi-pstab(ic)) )/dpsi
      enddo

      purs(1)=dpdpsi(iplas)
      purs(nurs)=dpdpsi(1)

      dpsi=1.d0/(nurs-1.d0)
      ppp(1)=0.d0

      do i=2,nurs
         ppp(i)=ppp(i-1) +(purs(i-1)+purs(i))*dpsi*0.5d0
      enddo

      return
      end subroutine retab_p

!----------------------------------------------------------------
      subroutine bt_tot(bettot)

      use compol, only: nt1, psin, iplas1, b0ax, psim, psip, cnor,
     & vol, vol1, vol2, vol3, vol4

      implicit none

      real*8, intent(out) :: bettot

      integer :: i, j
      real*8 :: volcen, volpl, psn, zpres, pintg, volk, paverg
      real*8, external :: funppp

      volcen=0.d0

      do j=2,nt1
         volcen=volcen+vol(1,j)*0.5d0
      enddo

      volpl=volcen
      psn=psin(1,2)
      zpres=funppp(psn)
      pintg=zpres*volcen

      do i=2,iplas1
         do j=2,nt1
            psn=psin(i,j)
            zpres=funppp(psn)
            volk=vol1(i,j)+vol2(i-1,j)+vol3(i-1,j-1)+vol4(i,j-1)
            pintg=pintg+zpres*volk
            volpl=volpl+volk
         enddo
      enddo

      paverg=pintg/volpl
      bettot=2.d0*paverg/(b0ax*b0ax)
      bettot=bettot*(psim-psip)*cnor

      return
      end subroutine bt_tot

!----------------------------------------------------------------
      subroutine presol(i_betp, betplx, betpol)

      use keys, only: kpr
      use compol, only: iter, itin, z, zm, psi, psim, psip, psipla,
     & iplas, iplas1, nt1, cnor, tokp, rm, fvac, f

      implicit none

      integer, intent(in) :: i_betp
      real*8, intent(in) :: betplx
      real*8, intent(out) :: betpol

      integer :: i, j, imov, imax, jmax
      real*8 :: erro, errod, errpsi, psimax, platok, rax, zax, psax

      imov=0
      iter=0
      itin=0

      call metrix
      call matcof
      call matpla

      do
         iter = iter+1
         itin = itin+1
         if(kpr.eq.1) then
            write(*,*)'iter=',iter,itin
         endif
         call rightg

	 if(i_betp.eq.1) then
	    if(iter.gt.4) call skbetp(betplx,betpol)
         endif

         call solint(imov)
         call remesh(erro,errpsi,imov)

         if(kpr.eq.1) then
            write(*,*)'presol: errpsi',errpsi
         endif

         errod = 0.5d0*erro/dabs(z(iplas,2)-zm)
         if(errpsi.lt.1.d-4 .or. iter.ge.50) EXIT

      enddo

      psimax=psi(1,2)
      imax=1
      jmax=2

      do i=2,iplas1
         do j=2,nt1
            if(psi(i,j).gt.psimax) then
               psimax=psi(i,j)
               imax=i
               jmax=j
            endif
         enddo
      enddo

      if(kpr.eq.1) then
         write(*,*) 'presol: '
         write(*,*) 'number of iterations ',iter
         write(*,*) 'accuracy ',errpsi
      endif

      call prgrid(imax,jmax,erro,1,errpsi)

      platok = tokp
      rax    = rm
      zax    = zm
      psax   = psim

      if(kpr.eq.1) then
         write(*,*) 'rax,zax',rax,zax
         write(*,*) 'plas.current',platok,cnor
         write(*,*) 'psax',psax
      endif

      psax   = psim
      psipla = psim-psip
      fvac   = f(iplas)

      call bt_pol(betpol)
      call wrb

      return
      end subroutine presol

!----------------------------------------------------------------
      subroutine put_tim(dt, time)

      use tim, only: dtim, ctim

      implicit none

      real*8, intent(in) :: dt, time

      dtim=dt
      ctim=time

      return
      end subroutine put_tim

!----------------------------------------------------------------
      subroutine get_tim(dt, time)

      use tim, only: dtim, ctim

      implicit none

      real*8, intent(out) :: dt, time

      dt=dtim
      time=ctim

      return
      end subroutine get_tim
