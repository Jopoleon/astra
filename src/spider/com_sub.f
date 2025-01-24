      subroutine grdef(igdf)

      use sp_parameters, only: nrp
      use compol, only: iplas, iplas1, psia, dpsda, q

      implicit none

      integer, intent(in) :: igdf

      integer :: i
      real*8 :: znam, psnor, zn_gp
      real*8, dimension(nrp) :: h_fi, fia

      if(igdf.eq.0) then
         do i=1,iplas
            psia(i)=(iplas-i)/(iplas-1.d0)
         enddo
         do i=1,iplas1
	    dpsda(i)=-1.d0
         enddo
      elseif(igdf.eq.1) then
	 do i=1,iplas
            psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
         enddo
	 do i=1,iplas1
            dpsda(i)=(1-2*i)/(iplas-1.d0)
         enddo
      elseif(igdf.eq.2) then
         psia(1)=0.d0
         do i=2,iplas
            znam=q(i-1)*(iplas-1)**2
            psia(i)=psia(i-1)+(2*i-3)/znam
         enddo
	 psnor=psia(iplas)
         do i=1,iplas
            psia(i)=1.d0-psia(i)/psnor
         enddo
         do i=1,iplas1
            dpsda(i)=-(2*i-1)/( (iplas-1.d0)*q(i)*psnor )
         enddo
      elseif(igdf.eq.3)then
         zn_gp=0.95d0
         fia(1)=0.d0
         psia(1)=0.d0
         do i=1,iplas-1
            if( i.lt.iplas/2 ) then
               h_fi(i)=1.d0
            elseif( i.ge.iplas/2 ) then
               h_fi(i)=h_fi(i-1)*zn_gp
            endif
            fia(i+1)=fia(i)+h_fi(i)
            psia(i+1)=psia(i)+(fia(i+1)**2-fia(i)**2)/q(i)
         enddo
         psnor=psia(iplas)
         do i=1,iplas
            psia(i)=1.d0-psia(i)/psnor
         enddo

      endif

      return
      end subroutine grdef

!----------------------------------------------------------------
      subroutine axdef(rma, zma, psima, dpm)

      use sp_parameters, only: ntp
      use compol, only: r, z, rm, zm, psi, nt1

      implicit none

      integer, parameter :: nshp=ntp+1

      real*8, intent(out) :: rma, zma, psima
      real*8, intent(out), dimension(*) :: dpm

      integer :: j, nsh
      real*8 :: rm0, zm0, det, erroma, psim0
      real*8, dimension(nshp) :: xs, ys, fun

      rm0=rm
      zm0=zm
      psim0=psi(1,1)

      nsh=1

      xs(nsh)=rm0
      ys(nsh)=zm0 
      fun(nsh)=psim0

      do j=2,nt1
         nsh=nsh+1
         xs(nsh)=r(2,j)
         ys(nsh)=z(2,j)
         fun(nsh)=psi(2,j)
      enddo

      call deriv5(xs,ys,fun,nsh,5,dpm)

      DET = dpm(3)*dpm(5) - dpm(4)**2

      Rma = Xs(1) + ( dpm(2)*dpm(4) - dpm(1)*dpm(5) )/DET
      Zma = Ys(1) + ( dpm(1)*dpm(4) - dpm(2)*dpm(3) )/DET

      erroma=dsqrt((rma-rm0)**2+(zma-zm0)**2)
      erroma=erroma/(dabs(r(2,2)-rm0))

      psima = fun(1) + dpm(1)*(rma-rm0) + dpm(2)*(zma-zm0) +
     &           0.5d0*dpm(3)*(rma-rm0)*(rma-rm0)  +
     &                 dpm(4)*(rma-rm0)*(zma-zm0)  +
     &           0.5d0*dpm(5)*(zma-zm0)*(zma-zm0)

      return
      end subroutine axdef

!----------------------------------------------------------------
      subroutine avr2_c(arr2, nro, nteta, arr1)

      use compol, only: iplas, nt1

      implicit none

      integer, intent(in) :: nro, nteta
      real*8, intent(in), dimension(nro, nteta) :: arr2
      real*8, intent(out), dimension(nro) :: arr1
      real*8, external :: avr1_c

      integer :: i, j
      real*8, dimension(nteta) :: arrw
      
      do i=1,iplas-1
         do j=2,nt1
            arrw(j)=arr2(i, j)
         enddo
         arr1(i)=avr1_c(arrw, i)
      enddo

      return
      end subroutine avr2_c

!----------------------------------------------------------------
      real*8 function avr1_c(arr1_c, i)

      use compol, only: nt1, vol

      implicit none

      integer, intent(in) :: i
      real*8, intent(in), dimension(*) :: arr1_c

      integer :: j
      real*8 :: avrg, sum_vol
      
      avrg=0.d0
      sum_vol=0.d0
      do j=2,nt1
         avrg=avrg+arr1_c(j)*vol(i,j)
         sum_vol=sum_vol+vol(i,j)
      enddo
      avr1_c=avrg/sum_vol
      
      return
      end function avr1_c

!----------------------------------------------------------------
      real*8 function avr_bnd(arr)

      use compol, only: nt1, dlt, iplas

      implicit none

      real*8, intent(in), dimension(*) :: arr

      integer :: j
      real*8 :: avrg, sum_len

      avrg=0.d0
      sum_len=0.d0
      do j=2,nt1
         avrg=avrg+0.5d0*(arr(j+1)+arr(j))*dlt(iplas,j)
         sum_len=sum_len+dlt(iplas,j)
      enddo
      avr_bnd=avrg/sum_len

      return
      end

!----------------------------------------------------------------
      subroutine bongri

      use sp_parameters, only: ntp
      use compol, only: iplas, nt1, r, z

      implicit none

      integer :: j, jb
      real*8 :: rr, zz, r0, r1, z0, z1, Fint
      real*8, dimension(ntp) :: dg_dn
      real*8, dimension(ntp, ntp) :: bin_adg

      common/com_bgr/ bin_adg, dg_dn

      do j=1,nt1
         rr=r(iplas, j)
         zz=z(iplas, j)
         do jb=1,nt1
            r0=r(iplas, jb)
            z0=z(iplas, jb)
            r1=r(iplas, jb+1)
            z1=z(iplas, jb+1)
            call bint(rr, zz, r0, z0, r1, z1, Fint, 1)
            bin_adg(jb,j)=fint
         enddo 
      enddo 

      return
      end subroutine bongri

!----------------------------------------------------------------
      subroutine cof_bon(cps_bon, bps_bon, dps_bon)

      use sp_parameters, only: nrp, ntp
      use compol, only: nt1, iplas, nt, psi, sq2, sq3, dlt, cur, r

      implicit none

      real*8, intent(out) :: cps_bon, bps_bon, dps_bon

      integer :: i, j, l, jb
      real*8 :: a1, a2, a3, g1, g2, g3, sqk, dltk, sum_len, avrg, djg,
     &   psb, ajg, bjg, r0, dgdnl
      real*8, dimension(ntp) :: dg_dn, aj, aj_G, bj, bj_g, dj, dj_G
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      real*8, dimension(ntp, ntp) :: bin_adg

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/com_bgr/ bin_adg, dg_dn

      i=iplas

      do j=2,nt1
         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)

         g1=psi(i-1,j-1)-psi(i,j-1)
         g2=psi(i-1,j)-psi(i,j)
         g3=psi(i-1,j+1)-psi(i,j+1)

         sqk=sq2(i-1,j)+sq3(i-1,j-1)
         dltk=(dlt(i,j-1)+dlt(i,j))*0.5d0
         dgdnl=-(a1*g1+a2*g2+a3*g3)+ cur(i,j)*sqk  !
         dg_dn(j)=dgdnl/dltk
         aj(j)=-(a1+a2+a3)/dltk
         r0=(3.d0*r(i,j)+r(i-1,j))*0.25d0
         bj(j)=sqk/dltk/r0
         dj(j)=sqk*r0/dltk
      enddo

      dg_dn(1)=dg_dn(nt1)
      dg_dn(nt)=dg_dn(2)

      aj(1)=aj(nt1)
      aj(nt)=aj(2)
      bj(1)=bj(nt1)
      bj(nt)=bj(2)
      dj(1)=dj(nt1)
      dj(nt)=dj(2)

      do l=2,nt1
         psb=0.d0
         ajG=0.d0
         bjG=0.d0
         djG=0.d0
         do jb=2,nt1
            psb=psb+bin_adg(jb,l)*(dg_dn(jb)+dg_dn(jb+1))*0.5d0
            ajG=ajG+bin_adg(jb,l)*(aj(jb)+aj(jb+1))*0.5d0
            bjG=bjG+bin_adg(jb,l)*(bj(jb)+bj(jb+1))*0.5d0
            djG=djG+bin_adg(jb,l)*(dj(jb)+dj(jb+1))*0.5d0
         enddo
         aj_G(l)=ajG
         bj_G(l)=bjG
         dj_G(l)=djG
      enddo

      aj_G(1)=aj_G(nt1)
      aj_G(nt)=aj_G(2)
      bj_G(1)=bj_G(nt1)
      bj_G(nt)=bj_G(2)
      dj_G(1)=dj_G(nt1)
      dj_G(nt)=dj_G(2)

      sum_len=0.d0
      do j=2,nt1
         sum_len=sum_len+dlt(iplas,j)
      enddo

      avrg=0.d0
      do j=2,nt1
         avrg=avrg+0.5d0*(aj_G(j+1)+aj_G(j))*dlt(iplas,j)
      enddo
      cps_bon=avrg/sum_len

      avrg=0.d0
      do j=2,nt1
         avrg=avrg+0.5d0*(bj_G(j+1)+bj_G(j))*dlt(iplas,j)
      enddo
      bps_bon=avrg/sum_len

      avrg=0.d0
      do j=2,nt1
         avrg=avrg+0.5d0*(dj_G(j+1)+dj_G(j))*dlt(iplas,j)
      enddo
      dps_bon=avrg/sum_len

      return
      end subroutine cof_bon

!----------------------------------------------------------------
      subroutine get_flfi(flfi_m)

      use compol, only: iplas, q, psia, psim

      implicit none

      real*8, intent(out) :: flfi_m

      integer :: i
       
      flfi_m=0.d0

      do i=1,iplas-1
         flfi_m=flfi_m-q(i)*(psia(i+1)-psia(i))*psim
      enddo
	             
      return
      end subroutine get_flfi

!----------------------------------------------------------------
      subroutine metcof(alp22, alp33, delsc, delv, ac0, skcen, cur_I,
     *                  alp33k, delsk, delvk)

      use sp_parameters, only: nrp, ntp
      use compol, only: iplas, iplas1, nt1, r, s, vol, cur,
     &  psia, psipla, sq1, sq2, sq3, sq4

      implicit none

      real*8, intent(out) :: ac0, skcen
      real*8, intent(out), dimension(nrp) :: alp22, alp33, delsc, delv,
     &   alp33k, delsk, delvk, cur_i

      integer :: i, j
      real*8 :: zdelsc, zdelsk, zdelv, zavrc, zavrk, zamn, za0, zapl,
     &   acj, sqk, dpsids, samn, sa0, sapl, diff, dcur_i, dcurj_i,
     &   r0, r1, r2, r3, r4, a1, a2, a3, a4, a5, a6, a7, a8, a9
      real*8, dimension(nrp) :: psf, sqtor, curj_I, deviat
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/compsf/ psf, sqtor

      do i=1,iplas1
         zdelsc=0.d0
         zdelv=0.d0
         zavrc=0.d0
         do j=2,Nt1
            r1=r(i,j)
            r2=r(i+1,j)
            r3=r(i+1,j+1)
            r4=r(i,j+1)
            r0=(r1+r2+r3+r4)*0.25d0
            zdelv=zdelv+vol(i,j)
            zdelsc=zdelsc+s(i,j)
            zavrc=zavrc+s(i,j)/r0
         enddo
         alp33(i)=zdelsc/zavrc
         delsc(i)=zdelsc
         delv(i)=zdelv
      enddo

      sqtor(1)=0.d0

      do i=2,iplas
         sqtor(i)=sqtor(i-1)+delsc(i-1)
      enddo

      do i=2,Iplas
         zdelv=0.d0
         zdelsk=0.d0
         zavrk=0.d0
         do j=2,Nt1
            if(i.ne.iplas) then
               sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
               r0=r(i,j)
            else
               sqk=sq2(i-1,j)+sq3(i-1,j-1)
               r0=(3.d0*r(i,j)+r(i-1,j))*0.25d0
            endif
            zdelsk=zdelsk+sqk
            zdelv=zdelv +sqk*r0
            zavrk=zavrk+sqk/r0
         enddo
         alp33k(i)=zdelsk/zavrk
         delsk(i)=zdelsk
         delvk(i)=zdelv
      enddo

      do i=2,iplas1
         samn=0.d0
         sa0 =0.d0
         sapl=0.d0

         do j=2,nt1
            a1=a13(i-1,j-1)
            a2=a34(i-1,j-1)+a12(i-1,j)
            a3=a24(i-1,j)
            a4=a23(i-1,j-1)+a14(i,j-1)
            a6=a14(i,j)+a23(i-1,j)
            a7=a24(i,j-1)
            a8=a34(i,j-1)+a12(i,j)
            a9=a13(i,j)
            a5=-(a1+a2+a3+a4+a6+a7+a8+a9)

            zamn=a1+a2+a3
            za0 =a4+a5+a6
            zapl=a7+a8+a9

            samn=samn+zamn
            sa0 =sa0 +za0
            sapl=sapl+zapl
         enddo

         alp22(i-1)=samn*delsc(i-1)
         if(i.eq.iplas1) alp22(i)=sapl*delsc(i)
      enddo

      do i=1,iplas1
         dpsids=psipla*(psia(i+1)-psia(i))/delsc(i)
         cur_I(i)=dpsids*alp22(i)
      enddo

! equation in central node

      ac0=0.d0
      skcen=0.d0

      do j=2,nt1
         acj=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         ac0=ac0-acj
         skcen=skcen+sq1(1,j)+sq4(1,j-1)
      enddo

! deviation

      curj_I(1)=cur(1,2)*skcen
      deviat(1)=cur_I(1)-curj_I(1)

      do i=2,Iplas1
         dcurj_I=0.d0
         do j=2,Nt1
            if(i.ne.iplas) then
               sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
            else
               sqk=sq2(i-1,j)+sq3(i-1,j-1)
            endif
            dcurj_I=dcurj_I+cur(i,j)*sqk
         enddo
         curj_I(i)=curj_I(i-1)+dcurj_I
         dcur_I=cur_I(i)-cur_I(i-1)
         deviat(i)=(dcur_I-dcurj_I)/delsk(i)
      enddo

      diff=(alp22(iplas-1)-alp22(iplas-2))/delsk(iplas-1)       
      alp22(iplas)=alp22(iplas-1)+0.5d0*diff*delsk(iplas)

      return
      end subroutine metcof

!----------------------------------------------------------------
      subroutine qst_b

      use sp_parameters, only: ntp
      use compol, only: nt1, iplas, r, z, r0ax, b0ax, fvac, f, dfdpsi,
     & psin, psim, psia, psip, s, sq2, sq3, flx_fi, q, flucfm

      implicit none

      integer, parameter :: nshp=ntp+1

      integer :: i, j, iefit, nsh
      real*8 :: fcefit, rc_efit, Drr, Drz, Dzz, Dxx, Dyy, dpsi, ps14,
     &   tg2a, cos2a, sin2a, ffp, dflufi, q2pi, r0, r1, r2, r3, r4
      real*8, dimension(5) :: dp
      real*8, dimension(nshp) :: xs, ys, fun

      common/efites/ fcefit, rc_efit, iefit

      nsh=1

      xs(nsh)=r(1,1)
      ys(nsh)=z(1,1)
      fun(nsh)=psin(1,1)

      do j=2,nt1
         nsh=nsh+1
         xs(nsh)=r(2,j)
         ys(nsh)=z(2,j)
         fun(nsh)=psin(2,j)
      enddo

      call deriv5(xs,ys,fun,nsh,5,dp)

      Drr=dp(3)
      Drz=dp(4)
      Dzz=dp(5)

      tg2a=2.d0*drz/(drr-dzz)
      cos2a=1.d0/dsqrt(1.d0+tg2a**2)
      sin2a=cos2a*tg2a

      Dxx=0.5d0*(Drr+Dzz)+0.5d0*cos2a*(Drr-Dzz)+sin2a*Drz
      Dyy=0.5d0*(Drr+Dzz)+0.5d0*cos2a*(Dzz-Drr)-sin2a*Drz

      if(iefit.eq.1) r0ax=rc_efit

      fvac=b0ax*r0ax

      dpsi=(psia(iplas)-psia(iplas-1))
      ps14=-0.25d0*dpsi
      ffp=dfdpsi(iplas)
      f(iplas-1)=dsqrt(fvac**2-ffp*(psim-psip)*dpsi)

      do i=iplas-2,1,-1
         dpsi=0.5d0*(psia(i+2)-psia(i))*(psim-psip)
         ffp=dfdpsi(i+1)
         f(i)=dsqrt(f(i+1)**2-2.d0*ffp*dpsi)
      enddo

      f(iplas)=fvac
      flucfm=0.d0
      flx_fi(1)=0.d0

      do i=1,iplas
         dflufi=0.d0
         do j=2,nt1
            if(i.ne.iplas) then
               r1=r(i,j)
               r2=r(i+1,j)
               r3=r(i+1,j+1)
               r4=r(i,j+1)
               r0=(r1+r2+r3+r4)*0.25d0
               dflufi=dflufi+s(i,j)*f(i)/r0
            else
               r0=r(i,j)
               dflufi=dflufi+(sq3(i-1,j-1)+sq2(i-1,j))*f(i)/r0
            endif
         enddo

         if(i.ne.iplas) then
            q2pi=-dflufi/((psia(i+1)-psia(i))*(psim-psip))
            flucfm=flucfm+dflufi
            flx_fi(i+1)=flx_fi(i)+dflufi
         else
            q2pi=-2.d0*dflufi/((psia(i)-psia(i-1))*(psim-psip))
         endif

         q(i)=q2pi
      enddo

      q(iplas)=1.5d0*q(iplas-1)-0.5d0*q(iplas-2)

      return
      end subroutine qst_b

!----------------------------------------------------------------
      subroutine EXTRPC(X, F, X0, X1, X2, X3, F1, F2, F3)

      implicit none

      real*8, intent(in) :: X, X0, X1, X2, X3, F1, F2, F3
      real*8, intent(out) :: F
      real*8 :: x1m, x2m, x3m, d1f2m, d1f3m, d2f23

      X1M = 0.5*(X0+X1)
      X2M = 0.5*(X1+X2)
      X3M = 0.5*(X2+X3)

      D1F2M =(F2-F1)/(X2M-X1M)
      D1F3M =(F3-F2)/(X3M-X2M)
      D2F23 = (D1F3M-D1F2M)/(X3M-X1M)
      F = F1+(X-X1M)*(D1F2M + (X-X2M)*D2F23)

      return
      end subroutine EXTRPC

!----------------------------------------------------------------
      subroutine EXTRP2(X, F, X1, X2, X3, F1, F2, F3)

      implicit none

      real*8, intent(in) :: X, X1, X2, X3, F1, F2, F3
      real*8, intent(out) :: F
      real*8 :: dfdx, d2fdx

      dfdx=(f3-f1)/(x3-x1)
      d2fdx=( (f3-f2)/(x3-x2) - (f2-f1)/(x2-x1) )/(x3-x1)
      
      f= f2+(x-x2)*( dfdx + d2fdx*(x-x2) )
      
      return
      end subroutine EXTRP2

!----------------------------------------------------------------
      real*8 function qvadin(x, xm, x0, xp, fm, f0, fp)

      implicit none

      real*8, intent(in) :: x, xm, x0, xp, fm, f0, fp

      qvadin = f0 + (fp - fm + 
     &    ( (fp-f0)/(xp-x0)-(f0-fm)/(x0-xm) )*(x-x0) ) *(x-x0)/(xp-xm)

      return
      end function qvadin

!----------------------------------------------------------------
      subroutine bt_pol(betpol)

      use sp_parameters, only: pi, amu0
      use compol, only: nt1, iplas1, sq1, sq2, sq3, sq4, tokp, cnor,
     & psin, psim, psip 

      implicit none

      real*8, intent(out) :: betpol

      integer :: i, j
      real*8 :: sqcen, sqk, psn, zpres, pintg
      real*8, external :: funppp

      sqcen=0.d0
      do j=2,nt1
         sqcen = sqcen + sq1(1,j) + sq4(1,j)
      enddo
      psn=psin(1,2)
      zpres = funppp(psn)
      pintg = zpres*sqcen

      do i=2,iplas1
         do j=2,nt1
            psn=psin(i,j)
            zpres=funppp(psn)
            sqk=sq1(i,j)+sq2(i-1,j)+sq3(i-1,j-1)+sq4(i,j-1)
            pintg=pintg+zpres*sqk
         enddo
      enddo

      betpol=8.d0*pi*pintg/(tokp*tokp)
      betpol=betpol*(psim-psip)*cnor/amu0**2

      return
      end subroutine bt_pol

!----------------------------------------------------------------
      subroutine skbetp(betplx, betpol)

      use compol, only: ngav, tok, cnor, tokff, iplas, dpdpsi

      implicit none

      real*8, intent(in) :: betplx
      real*8, intent(out) :: betpol

      integer :: i, nurs
      real*8 :: zcoin, coin

      call bt_pol(betpol)

      if(ngav.eq.0) then
         zcoin =  betplx / betpol
         zcoin =  (1.d0/zcoin-1.d0)*tok/(cnor*tokff)+1.d0
         coin =  1.d0/zcoin
         coin =  (coin-1.d0)*0.33d0+1.d0
      else
         coin=betplx/betpol
         do i=1,iplas
            dpdpsi(i)=coin*dpdpsi(i)
         enddo
      endif

      call taburs(1,coin,nurs)

      return
      end subroutine skbetp

!----------------------------------------------------------------
      subroutine cur_avg

      use sp_parameters, only: nrp, amu0
      use keys, only: kastr
      use jb, only: Bj_av, curfi_av
      use tim, only: dtim
      use compol, only: iplas, iplas1, nt1, vol1, vol2, vol3, vol4,
     & psim, psia, r, rm, tokp, flx_fi, cur, f

      implicit none

      integer :: i, j
      real*8 :: rm0, ac0, ac0n, skcen, skcen0, volk, volk_i, bj_i,
     &   b2_av_m, xa, fun_sig, fun_sig_a, sebeb, fm, ebm, em, see, 
     &   dpsdt05, dfidt05
      real*8, dimension(nrp) :: psi0, fi0, f0, ri0, q0, dpsidt, 
     &   dfidt, B2_av, EB, EB_c, Epar_c, WDj,
     &   BfJf, alfa22, alfa33, ds, dsk, dv, dvk, cj, alp33k,
     &   sigma, z_nvzk

      common/savt0/ psi0, fi0, f0, ri0, q0, 
     &              dpsidt, dfidt, rm0, ac0n, skcen0
      common /com_b2/ B2_av
      common /com_eb/ EB, EB_c, Epar_c
      common/com_heat_Dj/ WDj

      BfJf(1)=cur(1,2)*f(1)/rm

      do i=2,iplas1
         volk_i=0.d0
         BJ_i=0.d0
         do j=2,nt1
            volk=vol1(i,j)+vol2(i-1,j)+vol3(i-1,j-1)+vol4(i,j-1)
            volk_i=volk_i+volk
            BJ_i=BJ_i+cur(i,j)*(f(i)+f(i-1))*0.5d0*volk/r(i,j)
         enddo
         BfJf(i)=BJ_i/volk_i
      enddo

      call metcof (alfa22, alfa33, ds, dv, ac0, skcen, cJ,
     *             alp33k, dsk, dvk)

      curfi_av(1)=cur(1,2)
      BJ_av(1)=BfJf(1)
      B2_av(1)=(f(1)/rm)**2
      cJ(iplas)=tokp*amu0
      do i=2,iplas
         BJ_av(i)=(cJ(i)*F(i-1)-cJ(i-1)*F(i))/dvk(i)
         curfi_av(i)=(cJ(i)-cJ(i-1))/dsk(i)
      enddo
      do i=1,iplas1
         B2_av(i)=( F(i)*(flx_fi(i+1)-flx_fi(i))-  
     &            cJ(i)*psim*(psia(i+1)-psia(i)) )/dv(i)
      enddo
        
      B2_av(1)=(f(1)/rm)**2
      i=iplas

      B2_av(i)=(3.d0*B2_av(i-1) - B2_av(i-2))/2.0d0 
      B2_av_m=(3.d0*B2_av(1) - B2_av(2))/2.0d0 

      do i=1,iplas
         xa=dfloat(i-1)/dfloat(iplas-1)	  	  
         if(kastr.eq.1) then
            sigma(i) =dabs( fun_sig(xa,i)*dtim /amu0)  
         else
            sigma(i) =dabs( fun_sig_a(xa,i)*dtim/amu0 ) 
         endif
      enddo

      do i=2,iplas-1
         EB(i)= dpsidt(i)*(flx_fi(i+1)-flx_fi(i-1))
     &          -dfidt(i)*(psia(i+1)-psia(i-1))*psim
         EB(i)=EB(i)/(dv(i)+dv(i-1))/dtim
         SEBEB=sigma(i)*EB(i)**2
         Wdj(i)=SEBEB/B2_av(i)
      enddo

      i=1
	   
      Fm=(3.d0*F(1)-F(2))/2.0d0
      EBm=dpsidt(1)*Fm/rm**2/dtim
      EB(1)= EBm
      Em=dpsidt(1)/rm/dtim
      SEE=sigma(1)*Em**2
      Wdj(1)=SEE

      i=iplas
      EB(i)= dpsidt(i)*(flx_fi(i)-flx_fi(i-1))
     &       -dfidt(i)*(psia(i)-psia(i-1))*psim
      EB(i)=EB(i)/(dv(i)+dv(i-1))/dtim
      SEBEB=sigma(i)*EB(i)**2
      Wdj(i)=SEBEB/B2_av(i)

! (E,B) in cells

      do i=1,iplas-1
         Dpsdt05=0.5d0*(dpsidt(i)+dpsidt(i+1))/dtim
         Dfidt05=0.5d0*(dfidt(i) + dfidt(i+1))/dtim
         EB_c(i)= Dpsdt05*(flx_fi(i+1)-flx_fi(i))
     &           -Dfidt05*(psia(i+1)-psia(i))*psim
         EB_c(i)=-0.5d0*(EB(i)+EB(i+1))
         Epar_c(i)=EB_c(i)/dsqrt(B2_av(i))
      enddo

      do i=1,iplas
         z_nvzk(i)=amu0*sigma(i)*EB(i)+BJ_av(i)
      enddo

      return
      end subroutine cur_avg

!----------------------------------------------------------------
      function numlin(i, j, nr, nt)

      integer, intent(in) :: i, j, nr, nt

      integer :: jj

! definition of element number using it indexes(i,j) in 2-dim. array
      if(i.ne.1) then
         if(j.eq.1) then
            jj=nt-1
         elseif(j.eq.nt) then
            jj=2
         else
            jj=j
         endif
         numlin=(i-2)*(Nt-2)+jj
      else
         numlin=1
      endif

      return
      end function numlin

!----------------------------------------------------------------
      subroutine matpla

!    IL  -- number of equation (matrix line)
!    IM  -- element number in one-dimensional array A (a(im))
! IA(IL) -- number of first nonzero element in line IL
! JA(IM) -- number of matrix column for element a(im)

      use sp_parameters, only: nrp, ntp
      use compol, only: ia, ja, a, nt, nt1, nt2, nr, iplas, neqpla,
     & itin, nitbeg, nitdel, app0, dapp 

      implicit none

      integer :: i, j, il, im, km
      real*8 :: a1, a2, a3, a4, a5, a6, a7, a8, a9
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      integer, external :: numlin

      common/comaaa/ a12, a23, a34, a14, a13, a24

! equation for central point

      il=1
      im=1

      ia(il)=1
      ja(im)=1
      a(im)=0.d0

      do j=2,nt1
         im=im+1
         ja(im)=numlin(2,j,nr,nt)
         a(im)=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         a(1)=a(1)-a(im)
      enddo

      do 100 i=2,iplas-1

         j=2
         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)
         a4=a23(i-1,j-1)+a14(i,j-1)
         a6=a14(i,j)+a23(i-1,j)
         a7=a24(i,j-1)
         a8=a34(i,j-1)+a12(i,j)
         a9=a13(i,j)
         a5=-(a1+a2+a3+a4+a6+a7+a8+a9)

         il=il+1
         ia(il)=im+1

         if(i.ne.2) then
!2!cof. to (i-1,j)
            im=im+1
            a(im)=a2
            ja(im)=numlin(i-1,j,nr,nt)

!3!cof. to (i-1,j+1)
            im=im+1
            a(im)=a3
            ja(im)=numlin(i-1,j+1,nr,nt)

!1!cof. to (i-1,j-1)
            im=im+1
            a(im)=a1
            ja(im)=numlin(i-1,j-1,nr,nt)

!5!cof. to (i,j)
            im=im+1
            a(im)=a5
            ja(im)=numlin(i,j,nr,nt)

!6!cof. to (i,j+1)
            im=im+1
            a(im)=a6
            ja(im)=numlin(i,j+1,nr,nt)

!4!cof. to (i,j-1)
            im=im+1
            a(im)=a4
            ja(im)=numlin(i,j-1,nr,nt)

            if(i.ne.iplas-1) then

!8!cof. to (i+1,j)
               im=im+1
               a(im)=a8
               ja(im)=numlin(i+1,j,nr,nt)

!9!cof. to (i+1,j+1)
               im=im+1
               a(im)=a9
               ja(im)=numlin(i+1,j+1,nr,nt)

!7!cof. to (i+1,j-1)
               im=im+1
               a(im)=a7
               ja(im)=numlin(i+1,j-1,nr,nt)

            endif

         else

!1!cof. to central point
            im=im+1
            a(im)=a1+a2+a3
            ja(im)=numlin(1,j,nr,nt)

!5!cof. to (2,j)
            im=im+1
            a(im)=a5
            ja(im)=numlin(2,j,nr,nt)

!6!cof. to (2,j+1)
            im=im+1
            a(im)=a6
            ja(im)=numlin(2,j+1,nr,nt)

!4!cof. to (2,j-1)
            im=im+1
            a(im)=a4
            ja(im)=numlin(2,j-1,nr,nt)

!8!cof. to (3,j)
            im=im+1
            a(im)=a8
            ja(im)=numlin(3,j,nr,nt)

!9!cof. to (3,j+1)
            im=im+1
            a(im)=a9
            ja(im)=numlin(3,j+1,nr,nt)

!7!cof. to (3,j-1)
            im=im+1
            a(im)=a7
            ja(im)=numlin(3,j-1,nr,nt)

         endif

         do 110 j=3,nt2
            a1=a13(i-1,j-1)
            a2=a34(i-1,j-1)+a12(i-1,j)
            a3=a24(i-1,j)
            a4=a23(i-1,j-1)+a14(i,j-1)
            a6=a14(i,j)+a23(i-1,j)
            a7=a24(i,j-1)
            a8=a34(i,j-1)+a12(i,j)
            a9=a13(i,j)
            a5=-(a1+a2+a3+a4+a6+a7+a8+a9)

            il=il+1
            ia(il)=im+1

            if(i.ne.2) then

!1!cof. to (i-1,j-1)
               im=im+1
               a(im)=a1
               ja(im)=numlin(i-1,j-1,nr,nt)

!2!cof. to (i-1,j)
               im=im+1
               a(im)=a2
               ja(im)=numlin(i-1,j,nr,nt)

!3!cof. to (i-1,j+1)
               im=im+1
               a(im)=a3
               ja(im)=numlin(i-1,j+1,nr,nt)

!3!cof. to (i,j-1)
               im=im+1
               a(im)=a4
               ja(im)=numlin(i,j-1,nr,nt)

!5!cof. to (i,j)
               im=im+1
               a(im)=a5
               ja(im)=numlin(i,j,nr,nt)

!6!cof. to (i,j+1)
               im=im+1
               a(im)=a6
               ja(im)=numlin(i,j+1,nr,nt)

               if(i.ne.iplas-1) then

!7!cof. to (i+1,j-1)
                  im=im+1
                  a(im)=a7
                  ja(im)=numlin(i+1,j-1,nr,nt)

!8!cof. to (i+1,j)
                  im=im+1
                  a(im)=a8
                  ja(im)=numlin(i+1,j,nr,nt)

!9!cof. to (i+1,j+1)
                  im=im+1
                  a(im)=a9
                  ja(im)=numlin(i+1,j+1,nr,nt)

               endif

            else

!1!cof. to central point
               im=im+1
               a(im)=a1+a2+a3
               ja(im)=numlin(1,j,nr,nt)

!2!cof. to (2,j-1)
               im=im+1
               a(im)=a4
               ja(im)=numlin(2,j-1,nr,nt)

!3!cof. to (2,j)
               im=im+1
               a(im)=a5
               ja(im)=numlin(2,j,nr,nt)

!4!cof. to (2,j+1)
               im=im+1
               a(im)=a6
               ja(im)=numlin(2,j+1,nr,nt)

!5!cof. to (3,j-1)
               im=im+1
               a(im)=a7
               ja(im)=numlin(3,j-1,nr,nt)

!6!cof. to (3,j)
               im=im+1
               a(im)=a8
               ja(im)=numlin(3,j,nr,nt)

!7!cof. to (3,j+1)
               im=im+1
               a(im)=a9
               ja(im)=numlin(3,j+1,nr,nt)

            endif

 110     continue

         j=nt1

         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)
         a4=a23(i-1,j-1)+a14(i,j-1)
         a6=a14(i,j)+a23(i-1,j)
         a7=a24(i,j-1)
         a8=a34(i,j-1)+a12(i,j)
         a9=a13(i,j)
         a5=-(a1+a2+a3+a4+a6+a7+a8+a9)

         il=il+1
         ia(il)=im+1

         if(i.ne.2) then

!3!cof. to (i-1,j+1)
            im=im+1
            a(im)=a3
            ja(im)=numlin(i-1,j+1,nr,nt)

!1!cof. to (i-1,j-1)
            im=im+1
            a(im)=a1
            ja(im)=numlin(i-1,j-1,nr,nt)

!2!cof. to (i-1,j)
            im=im+1
            a(im)=a2
            ja(im)=numlin(i-1,j,nr,nt)

!6!cof. to (i,j+1)
            im=im+1
            a(im)=a6
            ja(im)=numlin(i,j+1,nr,nt)

!4!cof. to (i,j-1)
            im=im+1
            a(im)=a4
            ja(im)=numlin(i,j-1,nr,nt)

!5!cof. to (i,j)
            im=im+1
            a(im)=a5
            ja(im)=numlin(i,j,nr,nt)

            if(i.ne.iplas-1) then

!9!cof. to (i+1,j+1)
               im=im+1
               a(im)=a9
               ja(im)=numlin(i+1,j+1,nr,nt)

!7!cof. to (i+1,j-1)
               im=im+1
               a(im)=a7
               ja(im)=numlin(i+1,j-1,nr,nt)

!8!cof. to (i+1,j)
               im=im+1
               a(im)=a8
               ja(im)=numlin(i+1,j,nr,nt)

            endif

         else

!1!cof. to central point
            im=im+1
            a(im)=a1+a2+a3
            ja(im)=numlin(1,j,nr,nt)

!6!cof. to (2,j+1)
            im=im+1
            a(im)=a6
            ja(im)=numlin(2,j+1,nr,nt)

!4!cof. to (2,j-1)
            im=im+1
            a(im)=a4
            ja(im)=numlin(2,j-1,nr,nt)

!5!cof. to (2,j)
            im=im+1
            a(im)=a5
            ja(im)=numlin(2,j,nr,nt)

!9!cof. to (3,j+1)
            im=im+1
            a(im)=a9
            ja(im)=numlin(3,j+1,nr,nt)

!7!cof. to (3,j-1)
            im=im+1
            a(im)=a7
            ja(im)=numlin(3,j-1,nr,nt)

!8!cof. to (3,j)
            im=im+1
            a(im)=a8
            ja(im)=numlin(3,j,nr,nt)

         endif

 100  continue

      neqpla=il

      il=il+1
      ia(il)=im+1

      if( (itin/nitdel*nitdel+nitbeg) .EQ. itin ) then
          do km=1,im 
             app0(km)=a(km)
          enddo
      elseif(itin.gt.nitbeg) then
         do km=1,im
            dapp(km)=a(km)-app0(km)
         enddo
      endif

      return
      end subroutine matpla
