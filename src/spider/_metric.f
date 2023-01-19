      subroutine f_metric

      use sp_parameters

      implicit none

      integer :: i, j

      real*8 ::  dl12, dl14, dl23, dl34,
     &    r0,  r1,  r2,  r3,  r4,  r12,  r23,  r34,  r14,
     &    z0,  z1,  z2,  z3,  z4,  z12,  z23,  z34,  z14,
     &   dr0, dr1, dr2, dr3, dr4, dr12, dr23, dr34, dr14,
     &   dz0, dz1, dz2, dz3, dz4, dz12, dz23, dz34, dz14
      real*8 :: funsq

      include 'compol.inc'

      do i=1,Nr1
         do j=1,Nt
            DLr(i,j)=dsqrt( (r(i+1,j)-r(i,j))**2+(z(i+1,j)-z(i,j))**2 )
            St(i,j)=DLr(i,j)*(r(i+1,j)+r(i,j))*0.5d0
         enddo
      enddo

      do i=2,Nr
         do j=1,Nt1
            DLt(i,j)=dsqrt( (r(i,j+1)-r(i,j))**2+(z(i,j+1)-z(i,j))**2 )
            Sr(i,j)=DLt(i,j)*(r(i,j+1)+r(i,j))*0.5d0
         enddo
      enddo

      do i=1,Nr1
         do j=1,Nt1
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

            sq1(i,j)=funsq(r1,r12,r0,r14,z1,z12,z0,z14)
            sq2(i,j)=funsq(r2,r23,r0,r12,z2,z23,z0,z12)
            sq3(i,j)=funsq(r3,r34,r0,r23,z3,z34,z0,z23)
            sq4(i,j)=funsq(r4,r14,r0,r34,z4,z14,z0,z34)
            S(i,j)=funsq(r1,r2,r3,r4,z1,z2,z3,z4)

            vol1(i,j)=funsq(r1,r2,r4,r4,z1,z2,z4,z4)*(r1+r2+r4)/6.d0
            vol2(i,j)=funsq(r1,r2,r3,r3,z1,z2,z3,z3)*(r1+r2+r3)/6.d0
            vol3(i,j)=funsq(r2,r3,r4,r4,z2,z3,z4,z4)*(r2+r3+r4)/6.d0
            vol4(i,j)=funsq(r1,r3,r4,r4,z1,z3,z4,z4)*(r1+r3+r4)/6.d0
            vol(i,j)=(vol1(i,j)+vol2(i,j)+vol3(i,j)+vol4(i,j))

            dr12=r2-r1
            dz12=z2-z1
            dl12= (dr12**2+dz12**2)

            dr14=r4-r1
            dz14=z4-z1
            dl14= (dr14**2+dz14**2)

            dr23=r3-r2
            dz23=z3-z2
            dl23= (dr23**2+dz23**2)

            dr34=r3-r4
            dz34=z3-z4
            dl34= (dr34**2+dz34**2)

            cos2(i,j)=(dr12*dr23+dz12*dz23)/dsqrt(dl12*dl23)
            sin2(i,j)=1.d0-cos2(i,j)**2

            cos3(i,j)=(dr34*dr23+dz34*dz23)/dsqrt(dl34*dl23)
            sin3(i,j)=1.d0-cos3(i,j)**2

            if(i.eq.1) CYCLE

            cos4(i,j)=(dr34*dr14+dz34*dz14)/dsqrt(dl34*dl14)
            sin4(i,j)=1.d0-cos4(i,j)**2

            cos1(i,j)=(dr12*dr14+dz12*dz14)/dsqrt(dl12*dl14)
            sin1(i,j)=1.d0-cos1(i,j)**2

         enddo
      enddo

      return
      end subroutine f_metric

!----------------------------------------------------------------
      subroutine f_matcof

      use sp_parameters

      implicit none

      integer :: i, j
      real*8 :: s12, s14, s34, s23
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      include 'compol.inc'

      common/comaaa/ a12, a23, a34, a14, a13, a24

      do i=1,nr1
         do j=1,nt1
            s12=St(i,j)
            s14=Sr(i,j)
            s34=St(i,j+1)
            s23=Sr(i+1,j)

            if(i.ne.1) then
               a12(i, j) = ( 
     &            (-1.d0/s12+cos1(i,j)/s14)*vol1(i,j)/sin1(i,j) +
     &            (-1.d0/s12-cos2(i,j)/s23)*vol2(i,j)/sin2(i,j) )/s12

               a23(i, j) = ( 
     &            (-1.d0/s23-cos2(i,j)/s12)*vol2(i,j)/sin2(i,j) +
     &            (-1.d0/s23+cos3(i,j)/s34)*vol3(i,j)/sin3(i,j) )/s23

               a34(i, j) = ( 
     &            (-1.d0/s34+cos3(i,j)/s23)*vol3(i,j)/sin3(i,j) +
     &            (-1.d0/s34-cos4(i,j)/s14)*vol4(i,j)/sin4(i,j) )/s34

               a14(i, j) = ( 
     &            (-1.d0/s14+cos1(i,j)/s12)*vol1(i,j)/sin1(i,j) +
     &            (-1.d0/s14-cos4(i,j)/s34)*vol4(i,j)/sin4(i,j) )/s14

               a13(i,j) = (cos2(i,j)/(s12*s23))*vol2(i,j)/sin2(i,j) +
     &                    (cos4(i,j)/(s34*s14))*vol4(i,j)/sin4(i,j)

               a24(i, j) = (-cos1(i,j)/(s12*s14))*vol1(i,j)/sin1(i,j) +
     &                     (-cos3(i,j)/(s34*s23))*vol3(i,j)/sin3(i,j)

            else

               a12(i, j) = ( (-1.d0/s12 - cos2(i,j)/s23) * 
     &                          vol2(i,j)/sin2(i,j) )/s12

               a23(i, j) = ( 
     &            (-1.d0/s23-cos2(i,j)/s12)*vol2(i,j)/sin2(i,j) +
     &            (-1.d0/s23+cos3(i,j)/s34)*vol3(i,j)/sin3(i,j) )/s23

               a34(i, j) = ( (-1.d0/s34 + cos3(i,j)/s23) *
     &                          vol3(i,j)/sin3(i,j) )/s34

               a13(i, j) =  (cos2(i,j)/(s12*s23))*vol2(i,j)/sin2(i,j)
               a24(i, j) = (-cos3(i,j)/(s34*s23))*vol3(i,j)/sin3(i,j)

            endif
         enddo
      enddo

      return
      end subroutine f_matcof

!----------------------------------------------------------------
      subroutine f_matrix

!    IL  -- number of equation (matrix line)
!    IM  -- element number in one-dimensional array A (a(im))
! IA(IL) -- number of first nonzero element in line IL
! JA(IM) -- number of matrix column for element a(im)

      use sp_parameters

      implicit none

      integer :: i, j, il, im, km
      real*8 :: a1, a2, a3, a4, a5, a6, a7, a8, a9
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24
      integer, external :: numlin

      include 'compol.inc'
      include 'compol_add.inc'

      common /comaaa/ a12, a23, a34, a14, a13, a24

! equation for central point
      il=1
      im=1
      ia(1)=1
      ja(1)=1
      a(1)=0.d0

      do j=2,nt1
         im=im+1
         ja(im)=numlin(2,j,nr,nt)
         a(im)=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         a(1)=a(1)-a(im)
      enddo

! Main loop
      do i=2,nr1
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

            if(i.ne.nr1) then
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

         do j=3,nt2
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

               if(i.ne.nr1) then
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
         enddo

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

            if(i.ne.nr1) then
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

      enddo

      neq=il
      il=il+1
      ia(il)=im+1

      if( (itin/nitdel*nitdel+nitbeg) .EQ. itin ) then
         do km=1,im 
            aop0(km)=a(km)
         enddo
      elseif(itin.gt.nitbeg) then
         do km=1,im 
            daop(km)=a(km)-aop0(km)
         enddo
      endif

      return
      end subroutine f_matrix

!----------------------------------------------------------------
      subroutine f_procof(icq, cur_mu)

      use sp_parameters
      use keys, only: kstep

      implicit none

      integer, intent(in) :: icq
      real*8, intent(out) :: cur_mu

      real*8 :: tok_mu, psi_bon, psex_av, pspl_av, 
     &   bps_bon, cps_bon, dps_bon

      include 'compol.inc'

      tok_mu=tok*amu0

      call bongri
      call f_psib_ext(psex_av)
      call f_psib_pla(pspl_av)
      call cof_bon(cps_bon,bps_bon,dps_bon)
      psi_bon=psi_eav+pspl_av
      psibon=psi_bon
      psi_eav=psex_av

      if(ngav.eq.0) cur_mu=tok_mu

      if(ngav.eq.1) then
         call f_procof_I(icq)
      elseif(ngav.eq.2) then
         call f_procof_fl(icq)
      elseif(ngav.eq.-10 .AND. erru.lt.5.d-3) then
         call promat_j(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &     psipla, itin, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-1 .AND. erru.lt.5.d-3) then
         call promat(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm, 
     &     psipla, itin, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-2 .AND. erru.lt.5.d-2) then
         call promat(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm, 
     &     psipla, itin, kstep, fvac, psi_eav, pspl_av, psibon0, cur_mu)
      elseif(ngav.eq.-3 .AND. erru.lt.5.d-3) then
         call promat_I(iplas, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &     psipla, itin, kstep, fvac, psi_eav, pspl_av, psibon0, tok_mu)
      endif

      return
      end subroutine f_procof

!----------------------------------------------------------------
      subroutine f_procof_fl(icq)

      use sp_parameters

      implicit none

      integer, intent(in) :: icq

      integer :: i, j, npro
      real*8 :: a1, a2, a3, a4, a5, a6, a7, a8, a9, r0, r1, r2, r3, r4,
     &   cappa, cappl, capmn, ac0, acj, aexcen, avr14, skcen, dpsids, 
     &   psiai, psimx, qk, qn, q14, sqk, 
     &   zapro, zavrc, zavrk, zdelsc, zdelsk, zdelv, 
     &   af0, afpl, afmn, sa0, sapl, samn, za0, zapl, zamn, qmn, qpl
      real*8, dimension(nrp) :: amn, a0, apl, capp, avrc, delsc, avrk, 
     &   delsk, delv, bmn, b0, bpl, wrk1, wrk2, psf, rhs, sqtor
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      include 'compol.inc'
      include 'compol_add.inc'

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/compsf/ psf, sqtor

      psimx=psip+psipla

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
            zavrk=zavrk+sqk/r0
            zapro=zapro+aex(i,j)
         enddo

           avrk(i)=zavrk/zdelsk
           delsk(i)=zdelsk
           delv(i)=zdelv
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

         Qmn=-q(i-1)/(avrc(i-1)*delsc(i-1))
         Qpl=-q(i) / (avrc(i)*delsc(i))
         Qk=0.5d0*(Qmn+Qpl)
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
         rhs(i)=-dpdpsi(i+1)*delv(i+1)+rhs(i)
      enddo

      rhs(1)=rhs(1)+amn(2)*psimx
      rhs(npro)=rhs(npro)+apl(iplas-1)*psip

      call prog1d(npro,psf(2),bpl,bmn,b0,rhs,wrk1,wrk2)

      psf(1)=psimx
      psf(iplas)=psip

      do i=1,iplas1
         f(i)=-q(i)*(psf(i+1)-psf(i))/(delsc(i)*avrc(i))
      enddo

      do i=2,iplas1
         Qmn=-q(i-1)/(avrc(i-1)*delsc(i-1))
         Qpl=-q(i)/(avrc(i)*delsc(i))
         Qk=0.5d0*(qmn+qpl)
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

      aexcen=aex(1,2)
      dfdpsi(1) = ( (ac0*(psimx-psf(2)) )/skcen-rm*dpdpsi(1))*rm

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

      call extrpc( sqtor(i), capp(i), sqtor(i), sqtor(i-1), sqtor(i-2),
     &             sqtor(i-3), capp(i-1), capp(i-2), capp(i-3)  )

      cappa=capp(i)

      if(icq.eq.0) then
         DpsiDs = ( samn*(psip-psf(i-1))+avrk(i)*delsk(i)*dfdpsi(i) +
     &            delv(i)*dpdpsi(i) ) /cappa

         f(i)=dsqrt( f(i-1)**2+dfdpsi(i)*(psip-psf(i-1)) )
         q(i)=-f(i)*avrk(i)/DpsiDs
         write(*,*) 'q(iplas)=',q(i)
      endif

      Qn=-q(i)/avrk(i)
      Qmn=-q(i-1)/avrc(i-1)
      q14=-0.5d0*(q(i)+q(i-1))
      avr14=avrk(i)

      f(i) = (samn*(psip-psf(i-1)) + delv(i)*dpdpsi(i) - Q14*f(i-1)) /
     &       ( cappa/Qn - Q14 )

      dfdpsi(i) = (q14/avr14)*(f(i) - f(i-1))/delsk(i)
      psiai = psia(i)-0.25d0*(psia(i)- psia(i-1))
 
      call extrp2( psiai, dfdpsi(i), psia(i-3), psia(i-2), psia(i-1),
     &                    dfdpsi(i-3), dfdpsi(i-2), dfdpsi(i-1)  )

      f(i)=dsqrt(f(i-1)**2+dfdpsi(i)*(psip- psf(i-1)))

      return
      end subroutine f_procof_fl

!----------------------------------------------------------------
      subroutine f_procof_i(icq)

      use sp_parameters

      implicit none

      integer, intent(in) :: icq

      integer :: i, j
      real*8 :: a1, a2, a3, a4, a5, a6, a7, a8, a9, r0, r1, r2, r3, r4,
     &   zavrc, zavrk, zdelsc, zdelsk, zdelv, zff, zdfdps,
     &   deldf, delpsi, dftor, errdf, f2n1, sqk, qk, skcen,
     &   af0, afpl, afmn, sa0, sapl, samn, za0, zapl, zamn, qmn, qpl,
     &   qcappl, qcapmn, ac0, acj
     
      real*8, dimension(nrp) :: amn, a0, apl, Qcapp, capp, avrc, avrk, 
     &   delsc, delsk, delv, bmn, b0, bpl, wrk1, wrk2, psf, rhs, 
     &   sqtor, dfdpsn
      real*8, dimension(nrp, ntp) :: a12, a23, a34, a14, a13, a24

      include 'compol.inc'
      include 'compol_add.inc'

      common/comaaa/ a12, a23, a34, a14, a13, a24
      common/compsf/ psf, sqtor

      do i=1,iplas
         dfdpsn(i)=dfdpsi(i)
         psf(i)=psia(i)*(psim-psip)+psip 
      enddo

! dS(i+1/2) and <1/R>(i+1/2)

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
         avrk(i)=zavrk/zdelsk
         delsk(i)=zdelsk
         delv(i)=zdelv
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

      do

         if(itin.eq.1) then
            errdf=0.d0
            do i=2,iplas1
               zff=amn(i)*psf(i-1)+a0(i)*psf(i)+apl(i)*psf(i+1)
               zdfdps=( zff-dpdpsi(i)*delv(i) )/(avrk(i)*delsk(i))
               dfdpsi(i)=zdfdps
	       deldf=dabs(zdfdps-dfdpsn(i))
               errdf=dmax1(deldf,errdf)
               dfdpsn(i)=dfdpsi(i)
            enddo
	 endif

         errdf=0.d0
         f2n1=(tok*Qcapp(iplas-1))**2
         f(iplas-1)=dsqrt(f2n1)
         do i=iplas-1,2,-1
            f(i-1)=dsqrt( f(i)**2-dfdpsi(i)*(psf(i+1)-psf(i-1)) )
         enddo

         psf(iplas)=0.d0
	      
         do i=iplas1,1,-1
            dftor=f(i)*delsc(i)*avrc(i)
	    delpsi=dftor/q(i)
            psf(i)=psf(i+1)+delpsi 
         enddo
         do i=2,iplas1
            zff=amn(i)*psf(i-1)+a0(i)*psf(i)+apl(i)*psf(i+1)
            dfdpsi(i)=( zff-dpdpsi(i)*delv(i) )/(avrk(i)*delsk(i))
	    deldf=dabs(dfdpsi(i)-dfdpsn(i))
            errdf=dmax1(deldf,errdf)
            dfdpsn(i)=dfdpsi(i)
         enddo

         if(errdf.le.1.d-7) EXIT
      enddo
   
      i=iplas
	   
      call extrp2( psia(i), dfdpsi(i), psia(i-3), psia(i-2),
     &             psia(i-1), dfdpsi(i-3), dfdpsi(i-2), dfdpsi(i-1) )

! equation in central node

      ac0=0.d0
      skcen=0.d0

      do j=2,nt1
         acj=a12(1,j)+a24(1,j)+a34(1,j-1)+a13(1,j-1)
         ac0=ac0-acj
         skcen=skcen+sq1(1,j)+sq4(1,j-1)
      enddo

      dfdpsi(1) = ( ac0*(psf(1) - psf(2))/skcen - rm*dpdpsi(1) )*rm

      return
      end subroutine f_procof_i
