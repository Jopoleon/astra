      subroutine f_rightg

! right-hand side for problem L(g)=J

      use sp_parameters, only: neqp, amu0
      use keys, only: kstep, kastr, key_plc, kpr, key_prs
      use compol, only: nr, nt, nt1, neq, ngav, iplas, cur, right, erru,
     &    dpdpsi, dfdpsi, cnor, tok, tokp, tokff, tokpp, r, rm, psin,
     & sq1, sq2, sq3, sq4

      implicit none

      integer :: i, j, il
      real*8 :: cur_mu, curp, curf, sqk, r0, psn, curcen, sqcen
      integer, external :: numlin
      real*8, external :: tabf, tabp

      save cur_mu
        
      tokp = 0.d0

      do il=1,neqp
          right(il) = 0.d0
      enddo

      do i=1, Nr
          do j=1, Nt
              cur(i, j) = 0.d0
          enddo
      enddo

      if (kstep == 0 .OR. (ngav == 0 .AnD. kastr == 0)) then
! central point
          curp = tabp(1.d0)
          curf = tabf(1.d0)
          dpdpsi(1) = curp
          dfdpsi(1) = curf
          curcen = rm*curp + curf/rm
          sqcen = 0.d0
          do j=2, nt1
              sqcen = sqcen + sq1(1, j) + sq4(1, j)
          enddo

          tokp  = curcen*sqcen
          tokff = curf*sqcen/rm
          tokpp = curp*sqcen*rm

          il = numlin(1, 1, nt)
          right(il) = curcen*sqcen

! regular points

          do i=2, Iplas
              do j=2, Nt1
                  if (i /= iplas) then
                      sqk = sq1(i, j) + sq2(i-1, j) +
     &                      sq3(i-1, j-1) + sq4(i, j-1)
                      r0 = r(i, j)
                      psn = psin(i, j)
                  else
                      sqk = sq2(i-1, j) + sq3(i-1, j-1)
                      r0 = (3.d0*r(i, j) + r(i-1, j))*0.25d0
                      psn = (3.d0*psin(i, j) + psin(i-1, j))*0.25d0
                  endif

                  curp = tabp(psn)
                  curf = tabf(psn)
                  cur(i, j) = r0*curp + curf/r0
                  tokp  = tokp  + cur(i, j)*sqk
                  tokff = tokff + curf*sqk/r0
                  tokpp = tokpp + curp*sqk*r0
                  il = numlin(i, j, nt)
                  right(il) = cur(i, j)*sqk
              enddo
              dpdpsi(i) = curp
              dfdpsi(i) = curf
          enddo

          cnor = amu0*tok/tokp

          if (key_plc == 0) cnor = 1.d0
          cur_mu = amu0*tok

          if (kpr == 1) then
              write(*, *) 'cnor, tok, tokp'
              write(*, *) cnor, tok, tokp/amu0
          endif

          do j=1, Nt
              cur(1, j) = curcen*cnor
          enddo

          do i=2, iplas
              do j=2, Nt1
                  cur(i, j) = cur(i, j)*cnor
              enddo
          enddo

          do i=1, nr
              cur(i, 1 ) = cur(i, nt1)
              cur(i, nt) = cur(i, 2)
          enddo

          do il=1, neq
              right(il) = right(il)*cnor
          enddo

          do i=1, iplas
              dpdpsi(i) = dpdpsi(i)*cnor
              dfdpsi(i) = dfdpsi(i)*cnor
          enddo

          tokp = tok

      else
          if (kastr == 1 .AND. key_prs == 1 .AND. erru < 5.d-3) then
              call pres_d_psi
          endif
          call f_procof(1, cur_mu)

! central point

          curcen = rm*dpdpsi(1) + dfdpsi(1)/rm
          sqcen = 0.d0

          do j=2, nt1
             sqcen = sqcen + sq1(1, j) + sq4(1, j)
          enddo

          tokp = curcen*sqcen
          il = numlin(1, 1, nt)
          right(il) = curcen*sqcen

          do j=1, Nt
              cur(1, j) = curcen
          enddo

! regular points

          do i=2, Iplas
              do j=2, Nt1
                  if (i /= iplas) then
                      sqk = sq1(i, j) + sq2(i-1, j) +
     &                      sq3(i-1, j-1) + sq4(i, j-1)
                      r0 = r(i, j)
                  else
                      sqk = sq2(i-1, j) + sq3(i-1, j-1)
                      r0 = (3.d0*r(i, j) + r(i-1, j))*0.25d0
                  endif
                  curp = dpdpsi(i)
                  curf = dfdpsi(i)
                  cur(i, j) = r0*curp + curf/r0
                  tokp = tokp + cur(i, j)*sqk
                  il = numlin(i, j, nt)
                  right(il) = cur(i, j)*sqk
              enddo
          enddo

          if (ngav == 2) call toksur
          cnor = cur_mu/tokp
          if (key_plc == 0) cnor = 1.d0

          do i=1, nr
              cur(i, 1 ) = cur(i, nt1)
              cur(i, nt) = cur(i, 2)
          enddo

          tokp = tokp/amu0
 
          if (ngav == 0 ) then
              tokp=cur_mu/amu0
              do i=1, iplas
                  do j=1, Nt
                      cur(i, j)=cur(i, j)*cnor
                  enddo
              enddo
              do il=1, neq
                  right(il)=right(il)*cnor
              enddo
          endif

          if (kpr == 1) then
              write(*, *) 'cnor cur_mu', cnor, cur_mu
          endif

      endif

      return
      end subroutine f_rightg

!----------------------------------------------------------------
      subroutine f_rightp

      use sp_parameters, only: nrp, ntp
      use compol_add, only: dgdn, binadg, psii, g
      use compol, only: nr, nr1, nt, nt1, dlt, right

      implicit none

      integer :: i, j, jb, il, numlin
      real*8 :: a1, a2, a3, a7, a8, a9, g1, g2, g3, dltk, dgdnl, psb
      real*8, dimension(ntp) :: psib
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      common/comaaa/ a12, a23, a34, a14, a13, a24

      i = nr

      do j=2, nt1
          a1 = a13(i-1, j-1)
          a2 = a34(i-1, j-1) + a12(i-1, j)
          a3 = a24(i-1, j)
          g1 = g(i-1, j-1)
          g2 = g(i-1, j)
          g3 = g(i-1, j+1)
          dltk = (dlt(i, j-1) + dlt(i, j))*0.5d0
          dgdnl = a1*g1 + a2*g2 + a3*g3
          dgdn(j) = dgdnl/dltk
      enddo

      dgdn(1 ) = dgdn(nt1)
      dgdn(nt) = dgdn(2)

      do j=2, nt1
          psb = 0.d0
          do jb=2, nt1
              psb = psb + binadg(jb, j)*(dgdn(jb) + dgdn(jb+1))*0.5d0
          enddo
          psib(j) = psb
          psii(i, j) = -psb
      enddo

      psib(1 ) = psib(nt1)
      psib(nt) = psib(2)
      psii(i, 1 ) = -psib(nt1)
      psii(i, nt) = -psib(2)

      i = nr1
      do j=2, nt1
          a7 = a24(i, j-1)
          a8 = a34(i, j-1) + a12(i, j)
          a9 = a13(i, j)
          il = numlin(i, j, nt)
          right(il) = (a7*psib(j-1) + a8*psib(j) + a9*psib(j+1))
      enddo

      return
      end subroutine f_rightp

!---------------------------------------------------------------------
      subroutine rigext

      use sp_parameters, only: neqp, nrp, ntp
      use compol_add, only: psie
      use compol, only: nt, nt1, iplas, right

      implicit none

      integer :: i, j, il
      real*8 :: a7, a8, a9
      real*8, dimension(nrp) :: psib
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      integer, external :: numlin

      common/comaaa/ a12, a23, a34, a14, a13, a24

      do il=1, neqp
          right(il) = 0.d0
      enddo

      do j=2, nt1
          psib(j) = psie(iplas, j)
      enddo

      psib(1 ) = psib(nt1)
      psib(nt) = psib(2)

      i = iplas - 1
      do j=2, nt1
          a7 = a24(i, j-1)
          a8 = a34(i, j-1) + a12(i, j)
          a9 = a13(i, j)
          il = numlin(i, j, nt)
          right(il) = -(a7*psib(j-1) + a8*psib(j) + a9*psib(j+1))
      enddo

      return
      end

!----------------------------------------------------------------
      subroutine toksur

      use sp_parameters, only: nrp, ntp
      use compol_add, only: psii, psie, fpv
      use compol, only: nt, nt1, iplas, f, fvac, cur, curs, 
     & toksfi, right, dlt, r, sq2, sq3

      implicit none

      integer :: i, j, il
      real*8 :: a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     &   ps1, ps2, ps3, ps4, ps5, ps6, ps7, ps8, ps9, 
     &   DpiDni, DpiDne, DpeDni, DpsDni, DpsDne, sqk, dlt0, r0, fvv
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      integer, external :: numlin

      common/comaaa/ a12, a23, a34, a14, a13, a24

      i = iplas  
      toksfi = 0.d0
                     
      do j=2, nt1
          il = numlin(i, j, nt)
          a1 = a13(i-1, j-1)
          a2 = a34(i-1, j-1) + a12(i-1, j)
          a3 = a24(i-1, j)
          a4 = a23(i-1, j-1)
          a6 = a23(i-1, j)
          a5 = -(a1 + a2 + a3 + a4 + a6)

          ps1 = psii(i-1, j-1)
          ps2 = psii(i-1, j)
          ps3 = psii(i-1, j+1)
          ps4 = psii(i, j-1)
          ps5 = psii(i, j)
          ps6 = psii(i, j+1)

          DpiDni = a1*ps1 + a2*ps2 + a3*ps3 + a4*ps4 + a5*ps5 +
     &             a6*ps6 - right(il)

          ps1 = psie(i-1, j-1)
          ps2 = psie(i-1, j)
          ps3 = psie(i-1, j+1)
          ps4 = psie(i, j-1)
          ps5 = psie(i, j)
          ps6 = psie(i, j+1)

          DpeDni = a1*ps1 + a2*ps2 + a3*ps3 + a4*ps4 + a5*ps5 + a6*ps6

          a4 = a14(i, j-1)
          a6 = a14(i, j)
          a7 = a24(i, j-1)
          a8 = a34(i, j-1) + a12(i, j)
          a9 = a13(i, j)
          a5 = -(a4 + a6 + a7 + a8 + a9)

          ps4 = psii(i, j-1)
          ps5 = psii(i, j)
          ps6 = psii(i, j+1)
          ps7 = psii(i+1, j-1)
          ps8 = psii(i+1, j)
          ps9 = psii(i+1, j+1)

          DpiDne = -(a4*ps4 + a5*ps5 + a6*ps6 + a7*ps7 + a8*ps8 +
     &               a9*ps9)
          DpsDni = DpiDni + DpeDni
          DpsDne = DpiDne + DpeDni

          sqk = sq2(i-1, j) + sq3(i-1, j-1)
          dlt0 = (dlt(i, j-1) + dlt(i, j))*0.5d0
          r0 = r(i, j)
          curs(j) = (dlt0/r0**2)*fpv/(DpsDni + DpsDne)
          right(il) = right(il) + curs(j)*dlt0
          cur(i, j) = curs(j)*dlt0/sqk + cur(i, j)
          toksfi = toksfi + curs(j)*dlt0
      enddo

      curs(1 ) = curs(nt1)
      curs(nt) = curs(2)

      fvv = SQRT(f(iplas)**2 + fpv)
      write(6, *) 'Fp, Fvac, fv', f(iplas), Fvac, fvv

      return
      end subroutine toksur

!----------------------------------------------------------------
      subroutine f_psib_ext(psex_av)

      use sp_parameters, only: ntp
      use compol_add, only: psie
      use compol, only: nt, nt1, iplas

      implicit none

      real*8, intent(out) :: psex_av

      integer :: i, j
      real*8, dimension(ntp) :: psex_b
      real*8 :: avr_bnd

      i = iplas

      do j=1, nt
          psex_b(j) = psie(i, j)
      enddo

      psex_b(1 ) = psex_b(nt1)
      psex_b(nt) = psex_b(2)
      psex_av    = avr_bnd(psex_b)

      return
      end subroutine f_psib_ext

!----------------------------------------------------------------
      subroutine f_psib_pla(pspl_av)

      use sp_parameters, only: ntp
      use compol_add, only: psii
      use compol, only: nt, nt1, iplas

      implicit none

      real*8, intent(out) :: pspl_av

      integer :: i, j
      real*8, dimension(ntp) :: pspl_b
      real*8 :: avr_bnd

      i = iplas

      do j=2, nt1
          pspl_b(j) = psii(i, j)
      enddo

      pspl_b(1 ) = pspl_b(nt1)
      pspl_b(nt) = pspl_b(2)
      pspl_av    = avr_bnd(pspl_b)

      return
      end subroutine f_psib_pla

!----------------------------------------------------------------
      subroutine put_psib0(psi0_bnd)

      use compol, only: psibon0

      implicit none

      real*8, intent(in) :: psi0_bnd

      psibon0 = psi0_bnd
	            
      return
      end subroutine put_psib0
