      subroutine auto(rk, zk, tk, nk, nstep, ngrid,
     &                rk1, zk1, rk2, zk2,
     &                rk3, zk3, rk4, zk4,
     &                ntipe, necon, wecon )

      use iopath, only: path
      use comblc, only: nctrl, ni, nj, icont, nblm, rblm, zblm,
     &   rmin, rmax, zmin, zmax, rm0, zm0, 
     &   alp, qcen, psi_bon, 
     &   rx1, rx2, rx10, rx20, zx1, zx2, zx10, zx20

      implicit none

      integer, intent(in) :: nk, nstep, ngrid
      integer, intent(in), dimension(*) :: ntipe, necon
      real*8, intent(in), dimension(*) :: tk, wecon,
     &   rk, rk1, rk2, rk3, rk4, zk, zk1, zk2, zk3, zk4

      integer :: i, n_ctrl
      real*8 :: psi_bnd ! Where does the routine get psi_bnd from? git
      character(len=80) :: fname

! input initial data

      write(fname,'(a,a)') TRIM(path), '/data.dat'
      open(1,file=fname,form='formatted')
         read(1,*) icont
         read(1,*) ni
         read(1,*) nj
         read(1,*) rmin
         read(1,*) rmax
         read(1,*) zmin
         read(1,*) zmax
         read(1,*) alp
         read(1,*) nctrl
         read(1,*) qcen
         read(1,*) rx10
         read(1,*) zx10
         read(1,*) rx20
         read(1,*) zx20
         read(1,*) rm0
         read(1,*) zm0
      close(1)

      write(fname,'(a,a)') TRIM(path), '/data_d.wr'
      open(1,file=fname,form='formatted')
         write(1,*) icont, ni, nj, nctrl
         write(1,*) rmin, rmax, zmin, zmax, alp, qcen,
     &              rx10, zx10, rx20, zx20, rm0, zm0
      close(1)

      n_ctrl=nctrl
      psi_bon=psi_bnd

      rx1=rx10
      zx1=zx10
      rx2=rx20
      zx2=zx20

      write(fname,'(a,a)') TRIM(path), '/limpnt.dat'
      open(1,file=fname,form='formatted')
         read(1,*) nblm
         do i=1,nblm
            read(1,*) rblm(i),zblm(i)
         enddo
      close(1)

      write(fname,'(a,a)') TRIM(path), '/limpnt_d.wr'
      open(1,file=fname,form='formatted')
         write(1,*) nblm
         do i=1,nblm
            write(1,*) rblm(i),zblm(i)
         enddo
      close(1)

!index

      call glbind
      call grid_spid

      if(icont.eq.0) then
         call exfmat(rk, zk, tk, nk, rk1, zk1, rk2, zk2,
     &               rk3, zk3, rk4, zk4, ntipe, necon, wecon)
         call cfr_mat(rk, zk, tk, nk, NECON, WECON)
      endif

      return
      end subroutine auto

!----------------------------------------------------------------
      subroutine noauto

      use iopath, only: path
      use comblc, only: nctrl, ni, nj, icont, nblm, rblm, zblm,
     &   rmin, rmax, zmin, zmax, rm0, zm0, 
     &   alp, qcen, psi_bon, 
     &   rx1, rx2, rx10, rx20, zx1, zx2, zx10, zx20

      implicit none

      integer :: i, n_ctrl
      real*8 :: psi_bnd ! Used but not defined ! git
      character(len=80) :: fname

! input initial data

      write(fname,'(a,a)') TRIM(path), '/data_d.wr'
      open(1,file=fname,form='formatted')
         read(1,*) icont, ni, nj, nctrl
         read(1,*) rmin, rmax, zmin, zmax, alp, qcen,
     &             rx10, zx10, rx20, zx20, rm0, zm0
      close(1)

      n_ctrl=nctrl
      psi_bon=psi_bnd

      rx1=rx10
      zx1=zx10
      rx2=rx20
      zx2=zx20

      write(fname,'(a,a)') TRIM(path), '/limpnt_d.wr'
      open(1,file=fname,form='formatted')
         read(1,*) nblm
         do i=1,nblm
            read(1,*) rblm(i),zblm(i)
         enddo
      close(1)

      return
      end subroutine noauto

!----------------------------------------------------------------
      subroutine eq_0( pcequi, psitok, ncequi, nstep, ngrid, 
     &                 alf0, alf1, alf2, bet0, bet1, bet2, 
     &                 betplx, ngav1, 
     &                 ftok, tokout, psicen, nursb, 
     &                 psi_bnd, alp_b, rax, zax, n_ctrl, b_0, r_0 )

      use iopath, only: path
      use comblc, only: ni, nj, nnstpp, nctrl, icont, iter, iterbf, 
     &   itin, nitl, nitin, nrun, ix1, ix2, jx1, jx2, 
     &   alp, clr, clz, eps, tok, psi_bon, ucen, b0ax, r0ax, rm0, zm0

      implicit none

      integer, intent(in) :: ncequi, nstep, ngrid, ngav1, nursb
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2,
     &   betplx, ftok, psicen, psi_bnd, rax, zax, b_0, r_0
      real*8, intent(in), dimension(*) :: pcequi, psitok
      integer, intent(out) :: n_ctrl
      real*8, intent(out) :: tokout, alp_b

      integer :: isol, ien, nflag, nwr
      real*8 :: alf0p, alf1p, alf2p, bet0f, bet1f, bet2f, 
     &   eps0, epsin, epscrz, sigm, omg, ceps, coin
      character(len=80) :: fname

      alf0p=alf0
      alf1p=alf1
      alf2p=alf2
      bet0f=bet0
      bet1f=bet1
      bet2f=bet2
      b0ax=b_0
      r0ax=r_0
      rm0=rax 
      zm0=zax

      nnstpp=nstep
      n_ctrl=nctrl
      psi_bon=psi_bnd

      if(nstep.gt.0) icont=1
      tok  = ftok
      ucen = psicen

      write(fname,'(a,a)') TRIM(path), '/itpr.dat'
      open(1,file=fname,form='formatted')
         read(1,*) eps0
         read(1,*) eps
         read(1,*) epsin
         read(1,*) epscrz
         read(1,*) iterbf
         read(1,*) Nitl
         read(1,*) Nitin
         read(1,*) nrun
         read(1,*) nwr
         read(1,*) omg
         read(1,*) sigm
         read(1,*) ceps
      close(1)

! index
      call glbind

      ix1=1
      jx1=1
      ix2=1
      jx2=1

      call grid_spid
      call geom

      if(icont.eq.0) then
         call bndmat
         call wrdbnd
         icont=1
      else
         call rddbnd
      endif

      call extfil(pcequi,ncequi)

      coin=1.d0
      if(nstep.eq.0) then
         ien=1
         call taburs(ien,coin,nursb)
      endif

      isol  = 0
      nflag = 0

      if(ngrid.eq.0) then
         Ni=(Ni+1)/2
         Nj=(Nj+1)/2
         alp=alp-0.005d0

         call glbind
         call botlev(nstep)
         call geom
      endif

      call matrix
      call zero_ax(isol)

      clr=0.d0
      clz=0.d0

      if(nflag.eq.0) call psiful
      tokout = tok
      alp_b = alp
      call wrd

      iter=0
      itin=0

      return
      end subroutine eq_0

!----------------------------------------------------------------
      subroutine eq( pcequi, psicon, ncequi, nstep, ngrid,
     &               alf0, alf1, alf2, bet0, bet1, bet2,
     &               betpol, betplx,
     &               ngav1, tokout, pscout, 
     &               nursb, psi_bnd, alp_b, rax, zax )

      use iopath, only: path
      use keys, only: kpr
      use comblc, only: ni, nj, ix1, ix2, jx1, jx2, iter, itin, iterbf,
     &   imax, jmax, nitl, nitin, nrun, 
     &   alp, clr, clz, cnor, eps, erru, f_cur, g, qcen, tok,
     &   ui, um, up, r0ax, b0ax,
     &   rm, zm, rl, zl, dr, dz,
     &   rx1, rx2, rx10, rx20, zx0

      implicit none

      integer, intent(in) :: ncequi, nstep, ngav1, nursb
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2,
     &   betplx, psi_bnd, rax, zax
      real*8, intent(in) , dimension(*) :: pcequi
      real*8, intent(out), dimension(*) :: psicon
      integer, intent(inout) :: ngrid
      real*8, intent(out) :: betpol, tokout, pscout, alp_b

      integer :: ich, iflag, il, jl, isol, ien, nflag, nwr, itl,
     &   icelm, jcelm, nstp, istep
      real*8 :: eps0, epsin, epscrz, eps0l, epsl, epsinl, sigm, omg,
     &   coin, zcoin, rl0, zl0, r_ax, z_ax,
     &   ceps, cepsl, clr0, clr1, clz0, clz1, crz,
     &   dll, dllim, ddrr, ddzz, ddrl, ddzl, dcrdr, dcrdz, 
     &   dczdr, dczdz, det, delrl, delzl
      real*8 :: omega, sigma, alf0n
      character(len=80) :: fname

      common /comomg/ omega, sigma

      if(nstep.eq.0) then
         rm=rax		     
         zm=zax
      endif

      write(fname,'(a,a)') TRIM(path), '/itpr.dat'
      open(1,file=fname,form='formatted')
         read(1,*) eps0
         read(1,*) eps
         read(1,*) epsin
         read(1,*) epscrz
         read(1,*) iterbf
         read(1,*) Nitl
         read(1,*) Nitin
         read(1,*) nrun
         read(1,*) nwr
         read(1,*) omg
         read(1,*) sigm
         read(1,*) ceps
      close(1)

      eps0L=  eps0/3.d0
      epsL =  eps/3.d0
      cepsL = ceps/3.d0
      epsinL= epsin/3.d0

      if(nstep.gt.0) then
          iterbf=1
          goto 7890
      endif

! index
      call glbind

      ix1=1
      jx1=1
      ix2=1
      jx2=1

      call grid_spid
      call geom
      call rddbnd

 7890    continue

      call extfil(pcequi,ncequi)

      coin=1.d0
      if(nstep.eq.0) then
         ien=1
         call taburs(ien,coin,nursb)
      endif

      isol  = 0
      nflag = 0

      if(ngrid.eq.0) then
         Ni=(Ni+1)/2
         Nj=(Nj+1)/2
         alp=alp-0.005d0

         call glbind
         call botlev(nstep)
         call geom
      endif

 777  continue
      call matrix

      if(nstep.eq.0 .AND. nflag.eq.0) then
         call zero_ax(isol)
         clr=0.d0
         clz=0.d0
      endif

      if(nflag.eq.0) call psiful

      erru=1.d0
      omega=1.d0
      sigma=1.d0
      ich=0
      iflag=0

      call wrd

      itl=-1
      iter=0
      itin=0

 1000 continue

      iter=iter+1
      itin=itin+1

      if(iter.eq.iterbf .or. nstep.gt.0) then
         rl=rm
         zl=zm
         call shab(il,jl,icelm,jcelm)
      endif

      if(itin.le.3) then
         rx1=rx10
         rx2=rx20
      endif

      if(ich.eq.0) then

         if(erru.lt.eps .OR. itin.ge.Nitin) then
            itin=0
            omega=1.d0
            sigma=1.d0
            itl=itl+1
            ich=ich+1
            iflag = 1

            rl0=rl
            zl0=zl
            clr0=clr
            clz0=clz

            ddzl=dz(jmax)*ceps
            zl=zl+ddzl
            call shab(il,jl,icelm,jcelm)
         endif

      elseif(ich.eq.1) then

         if(erru.lt.epsin .OR. itin.ge.Nitin) then

            itin=0
            omega=1.d0
            sigma=1.d0
            ich=ich+1

            clr1=clr 
            clz1=clz 
            ddrl=dr(imax)*ceps
            rl=rl+ddrl
            zl=zl0
            call shab(il,jl,icelm,jcelm)
         endif      
                          
      elseif(ich.eq.2) then

         if(erru.lt.epsin .OR. itin.ge.Nitin) then
            itin=0 
            omega=1.d0
            sigma=1.d0

            dcrdr=(clr-clr0)/ddrl
            dczdr=(clz-clz0)/ddrl
            dcrdz=(clr1-clr0)/ddzl
            dczdz=(clz1-clz0)/ddzl
            det=dcrdr*dczdz-dczdr*dcrdz
                          
            delrl= (clz0*dcrdz-clr0*dczdz)/det
            delzl= (clr0*dczdr-clz0*dcrdr)/det                      
            dll=dsqrt(delrl**2 + delzl**2)
            dllim=0.25d0*dr(imax)
                          
            if(dll .gt. dllim) then
               nstp=dll/dllim
               ddrr=delrl/nstp
               ddzz=delzl/nstp
				        
               if(nstp.gt.10) nstp=10

               do istep=1,nstp
                  rl=rl0+ ddrr*istep
                  zl=zl0+ ddzz*istep
                  call shab(il,jl,icelm,jcelm)
                  rx1=rx10
                  rx2=rx20
                  call right0(il,jl,icelm,jcelm,ngav1)
                  call solve(isol,g)
                  call bound
                  call solve(isol,ui)
                  call psiful
               enddo

            else

               rl=rl0+ delrl
               zl=zl0+ delzl
               call shab(il,jl,icelm,jcelm)
            endif

            ich=0
            rx1=rx10
            rx2=rx20

         endif
      endif

      if(itl.ge.Nitl) then
         write(*,*) 'limit number of ext. iterations is exeeded'
         write(*,*) 'itl = ',itl,' Nitl = ',Nitl
         goto 1111
      endif

      if( ngav1.eq.4 .OR. ngav1.eq.5 ) then
         if(itin.ge.5) then
            call gridpl
            call qst(qcen,cnor,b0ax,r0ax)
         endif
      endif

      call right0(il,jl,icelm,jcelm,ngav1)
      call solve(isol,g)
      call bound
      call solve(isol,ui)
      call psiful

      if(ngav1.eq.1 .OR. ngav1.eq.3 .OR. ngav1.eq.5) then
         if(itin.ge.3) then
            call btpol( betpol )
            zcoin =  betplx / betpol
            zcoin =  (1.d0/zcoin-1.d0)*tok/(cnor*f_cur)+1.d0
            coin =  1.d0/zcoin
            coin =  (coin-1.d0)*0.33d0+1.d0
            call taburs( 1,coin,nursb)
         endif
      endif

      crz=dabs(clr*rm*rm/(um-up))+dabs(clz*(zm-zx0)/(um-up))

      if(erru.le.epsin .AND. crz.lt.epscrz .AND. itl.ge.0
     &   .AND. ich.eq.0) go to 1111

      omega = omg
      sigma = sigm

      if(erru.le.eps.and.kpr.eq.1)then
         write(*,'("iter erru rm zm ", i4,6(1pe12.5))'),
     &      iter,erru,rm,zm

         write(*,'("ich  clr clz ", i4,6(1pe12.5))'),
     &      ich,clr,clz
      endif

      if(iter.le.nrun) go to 1000

 1111  continue

      tokout = tok
      pscout = um
      r_ax = rm
      z_ax = zm
      alp_b = alp

      if(ngrid.eq.1) then
         call loop
         call wrd
         call psi_fil( psicon,ncequi )
      endif

      call btpol(betpol)

      if(ngrid.eq.0) then

         alp=alp+0.00499999d0

         ni=2*ni-1
         nj=2*nj-1

         imax=2*imax-1
         jmax=2*jmax-1

         ix1=2*ix1-1
         jx1=2*jx1-1

         ix2=2*ix2-1
         jx2=2*jx2-1

         call glbind
         call grid_spid
         call geom
         call uplev
         call extfil(pcequi,ncequi)
         call rddbnd
         call psiful
         call inter2(imax,jmax)
         call inter2(ix1,jx1)
         call inter2(ix2,jx2) 

         iterbf=1
         Nitl = 15
         ngrid=1
         nflag=1
         isol=0

         eps    = 1.0d-7
         eps0    = 2.0d-7
         epsin  = 1.0d-7
         epscrz = 5.0d-7
         ceps=0.02d0
         omega=1.d0
         sigma=1.d0

         goto 777

      endif

      call tabnor(cnor)
      alf0n=alf0
      write(*,*) '..................................................'
      write(*,*) 'Convergence of free bound. equilibrium iterations:'
      write(*,*) '  iter = ',iter,' nrun = ',nrun
      write(*,*) '  itl  = ',itl ,' Nitl = ',Nitl
      write(*,*) '  erru = ',erru
      write(*,*) '  crz  = ',crz

      return
      end subroutine eq

!----------------------------------------------------------------
      subroutine eq_ax( pcequi, psicon, ncequi, nstep, ngrid,
     &                  alf0, alf1, alf2, bet0, bet1, bet2,
     &                  betpol, betplx, zli3,
     &                  ngav1,
     &                  ftok, tokout, psicen, pscout,
     &                  EREVE0, ERPS,
     &                  psi_bnd, alp_b, rax, zax, isymm )

      use comblc, only: iter, itin, iterbf,
     &   nitl, nitin, nrun, icont, nnstpp,
     &   alp, alpnew, clr, clz, cnor, erru, f_cur, g, qcen, tok,
     &   ucen, ui, um, r0ax, b0ax, psi_bon,
     &   rm, zm, rl, zl, rx0, zx0

      implicit none

      integer, intent(in) :: ncequi, nstep, ngrid, ngav1, isymm
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2,
     &   betplx, ftok, psicen, psi_bnd, EREVE0
      real*8, intent(in) , dimension(*) :: pcequi
      real*8, intent(out), dimension(*) :: psicon
      real*8, intent(out) :: betpol, zli3, tokout, 
     &   pscout, alp_b, rax, zax, ERPS

      integer :: isol, nwr, il, jl, icelm, jcelm, nursb
      real*8 :: eps0, coin, zcoin
      real*8 :: omega, sigma, alf0n
      double precision :: alpef, alpnewef, rx0ef, zx0ef

      common /comomg/ omega, sigma
      common /comsav/ alf0n
      common /xpointvertpos/ alpef, alpnewef, rx0ef, zx0ef

      if(nstep.gt.0) icont=1
      nnstpp=nstep

      tok  = ftok
      ucen = psicen
      psi_bon = psi_bnd

      eps0   = EREVE0

      Nitl   = 9
      Nitin  = 61
      nrun   = 1000
      nwr    = 1000
      iterbf = 99999999

      call ext_fil(pcequi, ncequi)

      isol = 1

      clr=0.d0
      clz=0.d0

      call psiful

      erru=1.d0
      omega=1.d0
      sigma=1.d0

      iter=iter+1
      itin=itin+1

      rl=rm
      zl=zm

      call shab(il,jl,icelm,jcelm)

      if( ngav1.eq.4 .OR. ngav1.eq.5 ) then
         call gridpl
         call qst(qcen,cnor,b0ax,r0ax)
      endif

      call right0(il, jl, icelm, jcelm, ngav1)
      call solve(isol, g)  !gbound = 0
      call bound
      call solve(isol, ui) !uibound = dgdn integral

      if(isymm.eq.1) then ! symmetric case         
         call updown
      endif

      call psiful

      if( ngav1.eq.1 .OR. ngav1.eq.3 .OR. ngav1.eq.5 ) then
         call btpol( betpol )
         zcoin =  betplx / betpol
         zcoin =  (1.d0/zcoin-1.d0)*tok/(cnor*f_cur)+1.d0
         coin =  1.d0/zcoin
         coin =  (coin-1.d0)*0.33d0+1.d0
         call taburs( 1,coin,nursb) ! nursb used but not defined! git
      endif

      call psi_fil( psicon, ncequi )

      tokout = tok
      pscout = um
      rax = rm
      zax = zm
      alp_b = alp
      if(erru.le.eps0) then
         call loop
      endif

      ERPS=ERRU
      alpef=alp
      alpnewef=alpnew
      rx0ef=rx0
      zx0ef=zx0

      call tabnor(cnor)
      alf0n=alf0

      return
      end subroutine eq_ax

!----------------------------------------------------------------
      subroutine eq_ax0( pcequi, psicon, ncequi, nstep, ngrid,
     &                   alf0, alf1, alf2, bet0, bet1, bet2,
     &                   betplx,
     &                   ngav1,
     &                   ftok, psicen,
     &                   EREVE0,
     &                   psi_bnd)

      use comblc, only: nitl, nitin, icont, iter, iterbf, itin, nrun, 
     &   nnstpp, clr, clz, cnor, erru, qcen, psi_bon, b0ax, r0ax, 
     &   tok, ucen, rl, rm, zl, zm

      implicit none

      integer, intent(in) :: ncequi, nstep, ngrid, ngav1
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2,
     &   betplx, ftok, psicen, psi_bnd, EREVE0
      real*8, intent(in), dimension(*) :: pcequi, psicon

      integer :: il, jl, icelm, jcelm
      real*8 :: omega, sigma, alf0n
      double precision :: alpef, alpnewef, rx0ef, zx0ef

      common /comomg/ omega, sigma
      common /comsav/ alf0n
      common /xpointvertpos/ alpef, alpnewef, rx0ef, zx0ef

      if(nstep.gt.0) icont=1
      nnstpp=nstep

      tok  = ftok
      ucen = psicen
      psi_bon = psi_bnd

      Nitl   = 9
      Nitin  = 61
      nrun   = 1000
      iterbf = 99999999

      call ext_fil(pcequi,ncequi)

      clr=0.d0
      clz=0.d0

      call psiful

      erru=1.d0
      omega=1.d0
      sigma=1.d0

      iter=iter+1
      itin=itin+1

      rl=rm
      zl=zm

      call shab(il, jl, icelm, jcelm)

      if( ngav1.eq.4 .OR. ngav1.eq.5 ) then
         call gridpl
         call qst(qcen, cnor, b0ax, r0ax)
      endif

      call right0(il, jl, icelm, jcelm, ngav1)

      return
      end subroutine eq_ax0

!----------------------------------------------------------------
      subroutine eq_ax2( pcequi, psicon, ncequi, nstep, ngrid,
     &                   alf0, alf1, alf2, bet0, bet1, bet2,
     &                   betpol, betplx,
     &                   ngav1,
     &                   ftok, tokout, psicen, pscout,
     &                   EREVE0, ERPS,
     &                   psi_bnd, alp_b, rax, zax, isymm )

      use comblc, only: nitl, nitin, icont, iter, iterbf, itin, nrun, 
     &   nnstpp, alp, alp_b, alpnew, clr, clz, cnor, erru, f_cur, g, 
     &   psi_bon, tok, ucen, ui, um, 
     &   rx0, zx0, rl, rm, zl, zm

      implicit none

      integer, intent(in) :: ncequi, nstep, ngrid, ngav1, isymm
      real*8, intent(in) :: alf0, alf1, alf2, bet0, bet1, bet2,
     &   betplx, ftok, psicen, psi_bnd, EREVE0
      real*8, intent(in) , dimension(*) :: pcequi
      real*8, intent(out), dimension(*) :: psicon
      real*8, intent(out) :: betpol, tokout, 
     &   pscout, alp_b, rax, zax, ERPS

      integer :: isol, nwr, nursb
      real*8 :: eps0, coin, zcoin
      real*8 :: omega, sigma, alf0n
      double precision :: alpef, alpnewef, rx0ef, zx0ef

      common /comomg/ omega, sigma
      common /comsav/ alf0n
      common /xpointvertpos/ alpef, alpnewef, rx0ef, zx0ef

      if(nstep.gt.0) icont=1
      nnstpp=nstep

      tok  = ftok
      ucen = psicen
      psi_bon = psi_bnd

      eps0   = EREVE0

      Nitl   = 9
      Nitin  = 61
      nrun   = 1000
      nwr    = 1000
      iterbf = 99999999
      isol  = 1

      clr=0.d0
      clz=0.d0
      erru=1.d0
      omega=1.d0
      sigma=1.d0

      iter=iter+1
      itin=itin+1

      rl=rm
      zl=zm

      call solve(isol, g) !gbound = 0
      call bound
      call solve(isol, ui) !uibound = dgdn integral

      if(isymm.eq.1) then ! Symmetric case          
         call updown
      endif

      call psiful

      if( ngav1.eq.1 .OR. ngav1.eq.3 .OR. ngav1.eq.5 ) then
         call btpol( betpol )
         zcoin =  betplx / betpol
         zcoin =  (1.d0/zcoin-1.d0)*tok/(cnor*f_cur)+1.d0
         coin =  1.d0/zcoin
         coin =  (coin-1.d0)*0.33d0+1.d0
         call taburs(1, coin, nursb)
      endif

      call psi_fil(psicon, ncequi)

      tokout = tok
      pscout = um
      rax = rm
      zax = zm
      alp_b = alp

      if(erru.le.eps0) then
         call loop
      endif

      ERPS=ERRU
      alpef=alp
      alpnewef=alpnew
      rx0ef=rx0
      zx0ef=zx0

      call tabnor(cnor)

      alf0n=alf0

      return
      end subroutine eq_ax2

!----------------------------------------------------------------
      subroutine shab(il, jl, icelm, jcelm)

      use comblc, only: ni1, nj1, r, rl, z, zl, rmax

      implicit none

      integer, intent(out) :: il, jl, icelm, jcelm

      integer :: i, j, k, l, ic, jc, ilp, jlp
      real*8 :: ri, zj, sdmn, ddl

      do i=1,ni1
         ic=i
         if( rl.ge.r(i) .AND. rl.lt.r(i+1) ) EXIT
      enddo

      do j=1,nj1
         jc=j
         if( zl.ge.z(j) .AND. zl.lt.z(j+1) ) EXIT
      enddo

      icelm=ic
      jcelm=jc
      sdmn=rmax

      do k=0,1
         ilp=ic+k
         ri=r(ilp)
         do l=0,1
            jlp=jc+l
            zj=z(jlp)
            ddl=dsqrt( (rl-ri)**2 + (zl-zj)**2 )
            if(ddl.lt.sdmn) then
               sdmn=ddl
               il=ilp
               jl=jlp
            endif
         enddo
      enddo

      return
      end subroutine shab

!----------------------------------------------------------------
      subroutine glbind

      use comblc, only: ni, ni1, ni2, nj, nj1, nj2, nbnd

      implicit none

      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2

      nbnd=2*(ni1+nj1)

      return
      end subroutine glbind

!----------------------------------------------------------------
      subroutine botlev(nstep)

      use comblc, only: ni, nj, nbnd, r, z, ue, ui, binadg

      implicit none

      integer, intent(in) :: nstep

      integer :: i, j, ib, ibs

      do i=1,ni
         r(i)=r(2*i-1)
      enddo
      do j=1,nj
         z(j)=z(2*j-1)
      enddo

      do i=1,ni
         do j=1,nj
            ue(i,j)=ue(2*i-1,2*j-1)
            if(nstep.eq.0) CYCLE
            ui(i,j)=ui(2*i-1,2*j-1)
         enddo
      enddo

      do ib=1,nbnd
         do ibs=1,nbnd
            binadg(ibs,ib)=binadg(2*ibs-1,2*ib-1)+binadg(2*ibs,2*ib-1)
         enddo
      enddo

      return
      end subroutine botlev

!----------------------------------------------------------------
      subroutine uplev

      use comblc, only: ni, ni1, nj, nj1, ui, dr, dz

      implicit none

      integer :: i, j, nnil, nnjl

      nNil=(Ni+1)/2
      nNjl=(Nj+1)/2

      do i=nnil,1,-1
         do j=nnjl,1,-1
            ui(2*i-1,2*j-1)=ui(i,j)
         enddo
      enddo

      do i=1,ni,2
         do j=2,nj1,2
            ui(i,j)=(ui(i,j-1)*dz(j-1)+ui(i,j+1)*dz(j))/(dz(j-1)+dz(j))
         enddo
      enddo

      do i=2,ni1,2
         do j=1,nj
            ui(i,j)=(ui(i-1,j)*dr(i-1)+ui(i+1,j)*dr(i))/(dr(i-1)+dr(i))
         enddo
      enddo

      return
      end subroutine uplev

!----------------------------------------------------------------
      subroutine inter2(ip, jp)

      use comblc, only: ni2, nj2, r, z, u

      implicit none

      integer, parameter :: nshp=10

      integer, intent(in) :: ip, jp

      integer :: i, j, k, l, nsh
      real*8 :: rpi, zpj, upij, ri, zj
      real*8, dimension(5) :: dp
      real*8, dimension(nshp) :: xs, ys, fun

      if(ip.le.2 .OR. ip.gt.ni2 .OR. jp.le.2 .OR. jp.gt.nj2) return

      nsh=1
      xs(nsh)=r(ip)
      ys(nsh)=z(jp)
      fun(nsh)=u(ip,jp)

      do k=-2,2,2
         i= ip+k
         do l=-2,2,2
            j= jp+l
            if(i.eq.ip .AND. j.eq.jp) CYCLE
            nsh=nsh+1
            xs(nsh)=r(i)
            ys(nsh)=z(j)
            fun(nsh)=u(i,j)
         enddo
      enddo

      call deriv5(xs, ys, fun, nsh, 5, dp)

      rpi=r(ip)
      zpj=z(jp)
      upij=u(ip,jp)

      do k=-2,2
         i= ip+k
         ri=r(i)
         do l=-2,2
            j= jp+l
            zj=z(j)
            u(i,j)=upij + dp(1)*(ri-rpi) + dp(2)*(zj-zpj) +
     &                0.5*dp(3)*(ri-rpi)*(ri-rpi) +
     &                    dp(4)*(ri-rpi)*(zj-zpj) +
     &                0.5*dp(5)*(zj-zpj)*(zj-zpj)
         enddo
      enddo

      return
      end subroutine inter2

!----------------------------------------------------------------
      subroutine psiful

      use comblc, only: ni, nj, u, ue, ui, um, up, clr, clz, erru, r, z

      implicit none

      integer :: i, j
      real*8 :: uold, delu, del_um

      erru=0.d0

      do i=1,ni
         do j=1,nj
            uold=u(i,j)
            u(i,j)=ui(i,j)+ue(i,j)+clz*z(j)+clr*r(i)*r(i)
            delu=dabs(u(i,j)-uold)
            del_um=dabs(um-up)+1.d-9
            erru=dmax1(delu/del_um,erru)
         enddo
      enddo

      return
      end subroutine psiful

!----------------------------------------------------------------
      subroutine updown

      use comblc, only: ni, nj, ui

      implicit none

      integer :: i, j, nj12
      real*8 :: udown, uup, usym

      nj12= (nj+1)/2

      do j=1,nj12
         do i=1,ni
            udown=ui(i,j)
            uup=ui(i,nj-j+1)
	    usym=0.5d0*(uup+udown)
            ui(i,j)=usym
            ui(i,nj-j+1)=usym 
         enddo
      enddo

      return
      end subroutine updown

!----------------------------------------------------------------
      subroutine get_par(psi_bnd)

      use comblc, only: up

      implicit none

      real*8, intent(out) :: psi_bnd

      psi_bnd=up

      return
      end subroutine get_par
