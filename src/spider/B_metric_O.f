      subroutine metrix

      use compol, only: nt, iplas, iplas1, dlr, dlt, r, z, sr, st, nt1,
     &     cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4,
     & vol, vol1, vol2, vol3, vol4, sq1, sq2, sq3, sq4, s

      implicit none

      integer :: i, j
      real*8 :: r0, r1, r2, r3, r4, z0, z1, z2, z3, z4,
     &   r12, r14, r23, r34, z12, z14, z23, z34,
     &   sq_1, sq_2, sq_3, sq_4, dr12, dr14, dr23, dr34,
     &   dz12, dz14, dz23, dz34, dlm1, dlm2, dlm3, dlm4,
     &   cos_1, cos_2, cos_3, cos_4, sin_1, sin_2, sin_3, sin_4,
     &   vol_1, vol_2, vol_3, vol_4, dl12, dl14, dl23, dl34
      real*8, external :: funsq

      do j=1,Nt
         do i=1,iplas1
            DLr(i,j)=SQRT( (r(i+1,j)-r(i,j))**2+(z(i+1,j)-z(i,j))**2 )
            St(i,j)=DLr(i,j)*(r(i+1,j)+r(i,j))*0.5d0
         enddo
      enddo

      do j=1,Nt1
         do i=2,iplas
            DLt(i,j)=SQRT( (r(i,j+1)-r(i,j))**2+(z(i,j+1)-z(i,j))**2 )
            Sr(i,j)=DLt(i,j)*(r(i,j+1)+r(i,j))*0.5d0
         enddo
      enddo

      do j=1,Nt1
         do i=1,iplas1

            r1=r(i,j)
            r2=r(i+1,j)
            r3=r(i+1,j+1)
            r4=r(i,j+1)
            r0=(r1+r2+r3+r4)*0.25d0

            r12=(r1+r2)*0.5d0
            r23=(r3+r2)*0.5d0
            r34=(r3+r4)*0.5d0
            r14=(r1+r4)*0.5d0

            z1=z(i,j)
            z2=z(i+1,j)
            z3=z(i+1,j+1)
            z4=z(i,j+1)
            z0=(z1+z2+z3+z4)*0.25d0

            z12=(z1+z2)*0.5d0
            z23=(z3+z2)*0.5d0
            z34=(z3+z4)*0.5d0
            z14=(z1+z4)*0.5d0

            dl12=DLr(i,j)
            dl14=DLt(i,j)
            dl23=DLt(i+1,j)
            dl34=DLr(i,j+1)

            sq_1=funsq(r1,r12,r0,r14,z1,z12,z0,z14)
            sq_2=funsq(r2,r23,r0,r12,z2,z23,z0,z12)
            sq_3=funsq(r3,r34,r0,r23,z3,z34,z0,z23)
            sq_4=funsq(r4,r14,r0,r34,z4,z14,z0,z34)

            dr12=r2-r1
            dz12=z2-z1
            dr14=r4-r1
            dz14=z4-z1
            dr23=r3-r2
            dz23=z3-z2
            dr34=r3-r4
            dz34=z3-z4

            DLm2=dl12*dl23
            cos_2=(dr12*dr23+dz12*dz23)/DLm2
            sin_2=1.d0-cos_2**2
            vol_2=DLm2*SQRT(sin_2)*(r1+r2+r4)/12.d0

            DLm3=dl34*dl23
            cos_3=(dr34*dr23+dz34*dz23)/DLm3
            sin_3=1.d0-cos_3**2
            vol_3=DLm3*SQRT(sin_3)*(r2+r3+r4)/12.d0

            if(i.ne.1) then
               DLm4=dl34*dl14
               cos_4=(dr34*dr14+dz34*dz14)/DLm4
               sin_4=1.d0-cos_4**2
               vol_4=DLm4*SQRT(sin_4)*(r1+r3+r4)/12.d0
               DLm1=dl12*dl14
               cos_1=(dr12*dr14+dz12*dz14)/DLm1
               sin_1=1.d0-cos_1**2
               vol_1=DLm1*SQRT(sin_1)*(r1+r2+r4)/12.d0
            else
               cos_1=0.d0
               sin_1=1.d0
               vol_1=0.d0

               cos_4=0.d0
               sin_4=1.d0
               vol_4=0.d0
            endif

            cos2(i,j)=cos_2
            sin2(i,j)=sin_2
            vol2(i,j)=vol_2

            cos3(i,j)=cos_3
            sin3(i,j)=sin_3
            vol3(i,j)=vol_3

            cos4(i,j)=cos_4
            sin4(i,j)=sin_4
            vol4(i,j)=vol_4

            cos1(i,j)=cos_1
            sin1(i,j)=sin_1
            vol1(i,j)=vol_1

            vol(i,j)=vol_1+vol_2+vol_3+vol_4

            sq1(i,j)=sq_1
            sq2(i,j)=sq_2
            sq3(i,j)=sq_3
            sq4(i,j)=sq_4

            S(i,j)=sq_1+sq_2+sq_3+sq_4

         enddo
      enddo

      return
      end subroutine metrix

!----------------------------------------------------------------
      subroutine matcof

      use sp_parameters, only: nrp, ntp
      use compol, only: nt1, iplas1, st, sr, vol1, vol2, vol3, vol4,
     & cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4

      implicit none

      integer :: i, j
      real*8 :: s12, s14, s34, s23
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      common/comaaa/ a12, a23, a34, a14, a13, a24

      do j=1,nt1
         do i=1,iplas1
            s12=St(i,j)
            s14=Sr(i,j)
            s34=St(i,j+1)
            s23=Sr(i+1,j)

            if(i.ne.1) then

       a12(i,j)=( (-1.d0/s12 + cos1(i,j)/s14)*vol1(i,j)/sin1(i,j) +
     &            (-1.d0/s12 - cos2(i,j)/s23)*vol2(i,j)/sin2(i,j) )/s12

       a23(i,j)=( (-1.d0/s23 - cos2(i,j)/s12)*vol2(i,j)/sin2(i,j) +
     &            (-1.d0/s23 + cos3(i,j)/s34)*vol3(i,j)/sin3(i,j) )/s23

       a34(i,j)=( (-1.d0/s34 + cos3(i,j)/s23)*vol3(i,j)/sin3(i,j) +
     &            (-1.d0/s34 - cos4(i,j)/s14)*vol4(i,j)/sin4(i,j) )/s34

       a14(i,j)=( (-1.d0/s14 + cos1(i,j)/s12)*vol1(i,j)/sin1(i,j) +
     &            (-1.d0/s14 - cos4(i,j)/s34)*vol4(i,j)/sin4(i,j) )/s14

       a13(i,j)=  (cos2(i,j)/(s12*s23))*vol2(i,j)/sin2(i,j) +
     &            (cos4(i,j)/(s34*s14))*vol4(i,j)/sin4(i,j)

       a24(i,j)=  (-cos1(i,j)/(s12*s14))*vol1(i,j)/sin1(i,j) +
     &            (-cos3(i,j)/(s34*s23))*vol3(i,j)/sin3(i,j)

            else

       a12(i,j)=( (-1.d0/s12 - cos2(i,j)/s23)*vol2(i,j)/sin2(i,j) )/s12

       a23(i,j)=( (-1.d0/s23 - cos2(i,j)/s12)*vol2(i,j)/sin2(i,j) +
     &            (-1.d0/s23 + cos3(i,j)/s34)*vol3(i,j)/sin3(i,j) )/s23

       a34(i,j)=( (-1.d0/s34 + cos3(i,j)/s23)*vol3(i,j)/sin3(i,j) )/s34

       a13(i,j)=  (cos2(i,j)/(s12*s23))*vol2(i,j)/sin2(i,j)

       a24(i,j)=  (-cos3(i,j)/(s34*s23))*vol3(i,j)/sin3(i,j)

            endif
         enddo
      enddo

      return
      end subroutine matcof

!----------------------------------------------------------------
      subroutine procof(icq, cur_mu)

      use sp_parameters, only: amu0
      use keys, only: kstep
      use compol, only: ngav, tok, erru, iplas, iter, dpdpsi, dfdpsi,
     & f, fvac, flx_fi, psia, psim, psi_eav, psibon0, q, rm 

      implicit none

      integer, intent(in) :: icq
      real*8, intent(out) :: cur_mu

      real*8 :: pspl_av

      if(ngav.eq.0) cur_mu=tok*amu0
      call psib_pla(pspl_av)

      if(ngav.eq.1) then
         call procof_I
      elseif(ngav.eq.2) then
         call procof_fl(icq)
      elseif(ngav.eq.-10 .AND. erru.lt.5.d-3) then
         call promat_j(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &     psim, iter, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-1 .AND. erru.lt.5.d-3) then
         call promat(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm, 
     &     psim, iter, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-2 .AND. erru.lt.5.d-3) then
         call promat(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &     psim, iter, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-3 .AND. erru.lt.5.d-3) then
         cur_mu=tok*amu0
         call promat_I(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &     psim, iter, kstep, fvac, psibon0, cur_mu)
      endif

      return
      end subroutine procof

!----------------------------------------------------------------
      subroutine procof_fl(icq)

      use sp_parameters, only: nrp, ntp
      use compol, only: iplas, iplas1, iter, itin, nt1, q, f, rm,
     &  psia, psiax, psip, dpsda,
     & r, s, sq1, sq2, sq3, sq4, dpdpsi, dwdpsi, dfdpsi 

      implicit none

      integer, intent(in) :: icq

      integer :: i, j, npro
      real*8 :: r0, r1, r2, r3, r4, sqk, zavrk, zapro, zdelsk, qk,
     &   af0, afmn, afpl, sa0, samn, sapl, za0, zamn, zapl, qmn, qpl,
     &   a1, a2, a3, a4, a5, a6, a7, a8, a9, cappl, capmn, ac0, acj,
     &   cappa, skcen, dpsids, qn, q14, avr14, psiai, vrh, wgt, dmon,
     &   zdelv, zdelv3, zdelsc, zavrc
      real*8, dimension(nrp) :: amn, a0, apl, capp, avrc, delsc, 
     &   avrk, delsk, delv, bmn, b0, bpl, wrk1, wrk2, rhs, delv3,
     &   dfdpsn, delf, qs, qsn, psfn, dpsdas, psf, sqtor
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/compsf/ psf, sqtor

      do i=1,iplas
         dfdpsn(i)=dfdpsi(i)
         if(iter.eq.1) dfdpsn(i)=0.d0
         qs(i)=q(i)
         psfn(i)=psia(i)
         dpsdas(i)=dpsda(i)
      enddo

      do i=1,iplas
         qsn(i)=qs(i)
      enddo

      do i=1,iplas1
         zdelsc=0.d0
         zavrc=0.d0
         do j=2,Nt1
            r1=r(i,j)
            r2=r(i+1,j)
            r3=r(i+1,j+1)
            r4=r(i,j+1)
            r0=(r1+r2+r3+r4)*0.25d0
            zdelsc=zdelsc+s(i,j)
            zavrc=zavrc+s(i,j)/r0
         enddo
         avrc(i)=zavrc/zdelsc
         delsc(i)=zdelsc
      enddo

      sqtor(1)=0.d0

      do i=2,iplas
         sqtor(i)=sqtor(i-1)+delsc(i-1)
      enddo

      do i=2,Iplas
         zdelv=0.d0
         zdelv3=0.d0
         zdelsk=0.d0
         zavrk=0.d0
         zapro=0.d0
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
            zdelv3=zdelv3 +sqk*r0**3
            zavrk=zavrk+sqk/r0
            zapro=zapro
         enddo
         avrk(i)=zavrk/zdelsk
         delsk(i)=zdelsk
         delv(i)=zdelv
         delv3(i)=zdelv3
         rhs(i-1)=0.d0
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

         Qmn=-qs(i-1)/(avrc(i-1)*delsc(i-1))
         Qpl=-qs(i) / (avrc(i)*delsc(i))
         Qk=(Qmn*dpsdas(i-1)+Qpl*dpsdas(i))/(dpsdas(i-1)+dpsdas(i))
         afmn=Qk*Qmn*avrk(i)*delsk(i)
         afpl=Qk*Qpl*avrk(i)*delsk(i)
         af0=-(afmn+afpl)
         amn(i)=samn - afmn
         a0(i) =sa0  - af0
         apl(i)=sapl - afpl
         capp(i)=sapl*delsc(i)
         if(i.eq.iplas-1) then
            cappl=sapl*delsc(i)
            capmn=samn*delsc(i-1)
         endif
      enddo

      npro=iplas-2

      do i=1,npro
         if(i.ne.1) then
            bmn(i)=amn(i+1)
         else
            bmn(i)=0.d0
         endif
         b0(i)=-a0(i+1)
         if(i.ne.npro) then
            bpl(i)=apl(i+1)
         else
            bpl(i)=0.d0
         endif

         rhs(i)=rhs(i)-dpdpsi(i+1)*delv(i+1)-dwdpsi(i+1)*delv3(i+1)
      enddo

      rhs(1)=rhs(1)+amn(2)*psiax
      rhs(npro)=rhs(npro)+apl(iplas-1)*psip

      call prog1d(npro,psf(2),bpl,bmn,b0,rhs,wrk1,wrk2)

      psf(1)=psiax
      psf(iplas)=psip

      psfn(1)=1.d0
      psfn(iplas)=0.d0

      do i=2,iplas1
         psfn(i)=(psf(i)-psip)/(psiax-psip)
      enddo

      do i=1,iplas1
         f(i)=-qs(i)*(psf(i+1)-psf(i))/(delsc(i)*avrc(i))
      enddo

      do i=2,iplas1
         Qmn=-qs(i-1)/(avrc(i-1)*delsc(i-1))
         Qpl=-qs(i)/(avrc(i)*delsc(i))
         Qk=(Qmn*dpsdas(i-1)+Qpl*dpsdas(i))/(dpsdas(i-1)+dpsdas(i))
         dfdpsi(i)=Qk*(f(i)-f(i-1))
      enddo

! equation in central node

      ac0=0.d0
      skcen=0.d0

      do j=2,nt1
         acj=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         ac0=ac0-acj
         skcen=skcen+sq1(1,j)+sq4(1,j-1)
      enddo

      dfdpsi(1)=
     & ( (ac0*(psiax-psf(2)))/skcen-rm*dpdpsi(1)-dwdpsi(1)*rm**3 )*rm

! equation on bound node

      i=iplas
      samn=0.d0
      sa0 =0.d0

      do j=2,nt1
         a1=a13(i-1,j-1)
         a2=a34(i-1,j-1)+a12(i-1,j)
         a3=a24(i-1,j)
         zamn=a1+a2+a3
         za0 =-zamn
         samn=samn+zamn
         sa0 =sa0 +za0
      enddo

      call extrpc( sqtor(i),capp(i),
     &             sqtor(i),sqtor(i-1),sqtor(i-2),sqtor(i-3),
     &                      capp(i-1) ,capp(i-2) ,capp(i-3)  )

      cappa=capp(i)

      if(icq.eq.0) then
         DpsiDs=(
     &     samn*(psip-psf(i-1))+avrk(i)*delsk(i)*dfdpsi(i) +
     &    delv(i)*dpdpsi(i)+delv3(i)*dwdpsi(i) ) /cappa

         f(i)=SQRT( f(i-1)**2+dfdpsi(i)*(psip-psf(i-1)) )
         q(i)=-f(i)*avrk(i)/DpsiDs
         write(6,*) 'q(iplas)=',q(i)
      endif

      Qn=-q(i)/avrk(i)
      Qmn=-q(i-1)/avrc(i-1)

      q14=-0.5d0*(q(i)+q(i-1))

      avr14=avrk(i)

      f(i) = ( samn*(psip-psf(i-1)) +
     &     delv(i)*dpdpsi(i)+delv3(i)*dwdpsi(i)-Q14*f(i-1) )
     &     /( cappa/Qn-Q14 )

      dfdpsi(i)=(q14/avr14)*(f(i)-f(i-1))/delsk(i)

      psiai=psia(i)-0.25d0*(psia(i)- psia(i-1))

      call extrp2( psiai,dfdpsi(i),
     &             psia(i-3),psia(i-2),psia(i-1),
     &             dfdpsi(i-3) ,dfdpsi(i-2) ,dfdpsi(i-1)  )

      f(i)=SQRT(f(i-1)**2+dfdpsi(i)*(psip- psf(i-1)))

      if(itin.le.2) then
         vrh=1.00d0
      else
         vrh=1.00d0
      endif

      do i=1,iplas
         dfdpsi(i)=vrh*dfdpsi(i)+(1.d0-vrh)*dfdpsn(i)
      enddo

      do i=1,iplas1
         delf(i)=dfdpsi(i)-dfdpsn(i)
      enddo
      delf(iplas)=0.d0

      do i=2,iplas1
         wgt=1.0d0
         dmon=(delf(i+1)-delf(i))/(delf(i)-delf(i-1))
         if(dmon.lt.0.d0) wgt=0.5d0
         delf(i)=wgt*delf(i)+(1.d0-wgt)*
     &     (delf(i-1)*(psf(i+1)-psf(i)) + 
     &      delf(i+1)*(psf(i)-psf(i-1)))/(psf(i+1)-psf(i-1))
      enddo

      do i=2,iplas1
         dfdpsi(i)=dfdpsn(i)+delf(i)
      enddo

      return
      end subroutine procof_fl

!----------------------------------------------------------------
      subroutine procof_I

      use sp_parameters, only: nrp, ntp, amu0
      use compol, only: iplas, iplas1, itin, nt1, q, f, rm,
     &  psia, psip, psim, dpsda,
     & r, s, sq1, sq2, sq3, sq4, dpdpsi, dwdpsi, dfdpsi, tok

      implicit none

      integer :: i, j
      real*8 :: r0, r1, r2, r3, r4, sqk, zavrk, zdelsk, qk,
     &   af0, afmn, afpl, sa0, samn, sapl, za0, zamn, zapl, qmn, qpl,
     &   a1, a2, a3, a4, a5, a6, a7, a8, a9, ac0, acj,
     &   skcen, wgt, dmon,
     &   zdelv, zdelv3, zdelsc, zavrc, dftor, delpsi, deldf, zdfdps,
     &   f2n1, qcappl, qcapmn, errdf, zff
      real*8, dimension(nrp) :: amn, a0, apl, avrc, delsc, 
     &   avrk, delsk, delv, rhs, delv3,
     &   dfdpsn, delf, qcapp, psf, sqtor, dfdpsw
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/compsf/ psf, sqtor

      do i=1,iplas
         dfdpsw(i)=dfdpsi(i)
         dfdpsn(i)=dfdpsi(i)
         psf(i)=psia(i)*(psim-psip)+psip
      enddo

!   dS(i+1/2) and <1/R>(i+1/2)

      do i=1,iplas1
          zdelsc=0.d0
          zavrc=0.d0
          do j=2,Nt1
             r1=r(i,j)
             r2=r(i+1,j)
             r3=r(i+1,j+1)
             r4=r(i,j+1)

             r0=(r1+r2+r3+r4)*0.25d0

             zdelsc=zdelsc+s(i,j)
             zavrc=zavrc+s(i,j)/r0
         enddo
         avrc(i)=zavrc/zdelsc
         delsc(i)=zdelsc
      enddo

      sqtor(1)=0.d0

      do i=2,iplas
         sqtor(i)=sqtor(i-1)+delsc(i-1)
      enddo

! dS(i), dV(i) and <1/R>(i)

      do i=2,Iplas
         zdelv=0.d0
         zdelv3=0.d0
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
            zdelv3=zdelv3 +sqk*r0**3
            zavrk=zavrk+sqk/r0
         enddo
         avrk(i)=zavrk/zdelsk
         delsk(i)=zdelsk
         delv(i)=zdelv
         delv3(i)=zdelv3
         rhs(i-1)=0.d0  !right hand side for averaged equation
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

         Qmn=-q(i-1)/(avrc(i-1)*delsc(i-1))
         Qpl=-q(i) / (avrc(i)*delsc(i))
         Qk=(Qmn*dpsda(i-1)+Qpl*dpsda(i))/(dpsda(i-1)+dpsda(i))
         afmn=Qk*Qmn*avrk(i)*delsk(i)
         afpl=Qk*Qpl*avrk(i)*delsk(i)
         af0=-(afmn+afpl)
         amn(i)=samn !- afmn
         a0(i) =sa0  !- af0
         apl(i)=sapl !- afpl
         Qcapp(i)=Qpl/sapl

         if(i.eq.iplas-1) then
            Qcappl=Qpl/sapl
            Qcapmn=Qmn/samn
         endif

      enddo

      if(itin.eq.1) then
         errdf=0.d0
         do i=2,iplas1
            zff=amn(i)*psf(i-1)+a0(i)*psf(i)+apl(i)*psf(i+1)
            zdfdps=( zff-dpdpsi(i)*delv(i)-dwdpsi(i)*delv3(i) )
     &          /(avrk(i)*delsk(i))

            dfdpsi(i)=zdfdps
	    deldf=ABS(zdfdps-dfdpsn(i))
            errdf=MAX(deldf,errdf)
            dfdpsn(i)=dfdpsi(i)
            dfdpsi(i)=0.d0
         enddo
      endif

      do

         errdf=0.d0
         f2n1=(amu0*tok*Qcapp(iplas-1))**2
	 f(iplas-1)=SQRT(f2n1)

         do i=iplas-1,2,-1
            f(i-1)=SQRT( f(i)**2-dfdpsi(i)*(psf(i+1)-psf(i-1)) )
         enddo

         psf(iplas)=0.d0

         do i=iplas1,1,-1
            dftor=f(i)*delsc(i)*avrc(i)
            delpsi=dftor/q(i)
            psf(i)=psf(i+1)+delpsi
         enddo

         do i=2,iplas1
            zff=amn(i)*psf(i-1)+a0(i)*psf(i)+apl(i)*psf(i+1)
            dfdpsi(i)=( zff-dpdpsi(i)*delv(i)-dwdpsi(i)*delv3(i) )
     &             /(avrk(i)*delsk(i))
            deldf=ABS(dfdpsi(i)-dfdpsn(i))
            errdf=MAX(deldf,errdf)
            dfdpsn(i)=dfdpsi(i)
         enddo

         if(errdf.le.1.d-7) EXIT
      enddo

      i=iplas

      call extrp2( psia(i), dfdpsi(i),
     &             psia(i-3), psia(i-2), psia(i-1),
     &             dfdpsi(i-3), dfdpsi(i-2), dfdpsi(i-1)  )

! equation in central node

      ac0=0.d0
      skcen=0.d0

      do j=2,nt1
         acj=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         ac0=ac0-acj
         skcen=skcen+sq1(1,j)+sq4(1,j-1)
      enddo

      dfdpsi(1)=( ac0*(psf(1)-psf(2))/skcen
     &           -rm*dpdpsi(1)-dwdpsi(1)*rm**3 )*rm

      do i=1,iplas1
         delf(i)=dfdpsi(i)-dfdpsw(i)
      enddo
      delf(iplas)=0.d0

      do i=2,iplas1
         wgt=1.0d0
         dmon=(delf(i+1)-delf(i))/(delf(i)-delf(i-1))
         if(dmon.lt.0.d0) wgt=0.5d0
         delf(i)=wgt*delf(i)+(1.d0-wgt)*
     & (delf(i-1)*(psf(i+1)-psf(i))+delf(i+1)*(psf(i)-psf(i-1)))/
     &          (psf(i+1)-psf(i-1))
      enddo

      do i=2,iplas1
         dfdpsi(i)=dfdpsw(i)+delf(i)
      enddo

      return
      end subroutine procof_I
