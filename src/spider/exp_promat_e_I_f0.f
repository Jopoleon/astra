      subroutine promat(n, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm, 
     &                  psim, iter, k_step, fvac, psi_eav, pspl_av,
     &                  psibon0, cur_mu)
!      diffusion & equilibrium both equations solvinq 
!      Newton metod with linearisation
! error-correcting

      use sp_parameters, only: nrp
      use keys, only: kastr, key_fixfree

      implicit none

      integer, intent(in) :: n, iter, k_step
      real*8, intent(in) :: rm, psim, fvac, psi_eav, pspl_av
      real*8, intent(in), dimension(n) :: flx_fi, psia
      real*8, intent(out) :: cur_mu, psibon0
      real*8, intent(out), dimension(n) :: dfdpsi, dpdpsi, f, q

      integer :: i, k, ks, key_boncon, k_step_prev, jumpstep
      real*8 :: del_psb, psi_bn1, rm0, ac0n, skcen0, xa, ffprim,
     &   alfa22n, epsel, dt, fbnd, skcen, dpsdt, d2psi, d2fi,
     &   ac0, fk0, rik0, qk0, dri0, df0, psi0n, fi0n, fin, fin1,
     &   cap_psi, cap_fi, capfi, capsi, pdv, pbb, cbutn,
     &   pspl_tst, psi_ext, psiext, psi_bnd, psin, psin1,
     &   bps_bon, cps_bon, dps_bon, av, bv, cv,
     &   a11, a12, a21, a22, c11, c12, c21, c22, h1, h2,
     &   ds1, dtsig, deltamax, ddivmax, abd, aby, znv, f_wght,
     &   dfldt, fn, fn1, rin, rin1, dsn, dsn1,
     &   al22n, al22n1, al22x, al33n, al33n1
      real*8, dimension(512) :: wrk
      real*8, dimension(nrp) :: psi, fi, ri, alfa22, alfa33, alp33k, ds,
     &   dsk, dv, dvk, dpsi, dfi, dri, df, fk, rik, qk, delf, dfdpsn, 
     &   psi_d, fi_d, f_d, ri_d, ps_pnt,
     &   psi0, fi0, f0, ri0, q0, dpsidt, dfidt,
     &   sigma, cbut_b
      real*8, dimension(nrp, 2) :: y, h, zn
      real*8, dimension(nrp, 4) :: a, b, c
      real*8, external :: fun_jb, fun_jb_a, fun_sig, fun_sig_a

      common /selcon/ psi_d, fi_d, f_d, ri_d, ps_pnt, del_psb, psi_bn1
      common /savt0/ psi0, fi0, f0, ri0, q0, dpsidt, dfidt, 
     &               rm0, ac0n, skcen0
      common /com_but/ sigma, cbut_b

      save k_step_prev

      key_boncon=1
      if(key_boncon.eq.0) then
         dpdpsi(n)=0.d0
      endif
      epsel = 1.0d-10
      dt = 1.d0
      cap_psi=0.0d0   !1.d0
      cap_fi =0.d0
      fbnd = fvac

      call cof_bon(cps_bon,bps_bon,dps_bon)
      pspl_tst = -cps_bon*(psia(n) - psia(n-1))*psim +
     &            bps_bon*dfdpsi(n) + dps_bon*dpdpsi(n)

      psi_ext=psi_eav
      psi_bnd=psi_ext + pspl_av

! metric coefficients
      call metcof(alfa22, alfa33, ds, dv, ac0, skcen, ri,
     &            alp33k, dsk, dvk)

      do i=1,n
         dfdpsn(i)=dfdpsi(i)
         psi(i)=psim*psia(i)+psi_bnd
         fi(i)=flx_fi(i)
         if(k_step.ne.k_step_prev) then
            jumpstep=k_step - k_step_prev
            del_psb=psi_bn1 - psibon0
            if((jumpstep.eq.1 .AnD. k_step.gt.1) .OR. 
     &          key_fixfree.eq.0) then
               psibon0=psi_bn1	   	            
            endif
            psi0(i)=psim*psia(i)+psibon0
            fi0(i)=fi(i)
            F0(i) =f(i)
            rI0(i) =rI(i)
            Q0(i) =q(i)
            dpsidt(i)=0.d0
            dfidt(i )=0.d0
         endif
         xa=dsqrt(flx_fi(i)/flx_fi(n))	  	  
         if(kastr.eq.1) then
            sigma(i) =fun_sig(xa,i)  
            cbut_b(i) =fun_jb(xa,i)  
         else
            sigma(i) =fun_sig_a(xa,i)  
            cbut_b(i) =fun_jb_a(xa,i)  
         endif
      enddo

      if(key_boncon.eq.0) then
         cbut_b(n)=0.d0
         sigma(n)=-1.d-6
      endif

      if(iter.eq.1) then
         rm0=rm
         skcen0=skcen
         ac0n=ac0
      endif	  

      do i=1,n
         wrk(i)=dsqrt(dabs(alfa22(i)))   
      enddo

      do i=1,n-1
         f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)   
         ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)   
         q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))   
      enddo

      F(N)=fbnd 
      ffprim = rm*( ac0*(psi(1)-psi(2))/skcen - rm*dpdpsi(1) )
      delf(1) = ffprim - dfdpsi(1) 
      dfdpsi(1)=ffprim
      do i=2,n-1
         ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))   
         delf(i) = ffprim - dfdpsi(i) 
         dfdpsi(i) =  ffprim  
      enddo
	  
      alfa22n=alfa22(n)
      dfdpsi(n) =  (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))

      do ks=1,100             ! iteration cycle begins
         do i=2,n-1
            dpsi(i)=0.5d0*(psi(i+1)-psi(i-1))
            dfi(i) =0.5d0*(fi(i+1)-fi(i-1))
            dri(i) =ri(i)-ri(i-1)
            df(i)  =f(i)-f(i-1)
            Fk(i) =0.5d0*(f(i)+f(i-1))
            rIk(i) =0.5d0*(rI(i)+rI(i-1))
            Qk(i) =-(dsk(i)/alp33k(i))*(f(i)+f(i-1))/(psi(i+1)-psi(i-1))
         enddo

! matrix completing:   a(i)*y(i-1)-c(i)*y(i)+b(i)*y(i+1)=-h(i)
         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
               capfi= cap_fi 
            endif
            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif

            d2psi=0.5d0*(psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi=0.5d0*(fi(i+1)-2.d0*fi(i)+fi(i-1))        

            a(i,1) = 0.5d0*dfidt(i) - ( Fk(i)+0.5d0*df(i) ) *
     &   	        ( alfa22(i-1)/ds(i-1) )/sigma(i) -
     &   	     0.5d0*dfidt(i)*capsi 

            a(i,2) = -0.5d0*dpsidt(i) + ( rIk(i)+0.5d0*dri(i) ) *
     &	                   ( alfa33(i-1)/ds(i-1) )/sigma(i) +
     &   	      0.5d0*dpsidt(i)*capfi 

            a(i,3) = -0.5d0*dri(i)/dpsi(i) + alfa22(i-1)/ds(i-1)

            a(i,4)= Qk(i)*alfa33(i-1)/ds(i-1) +
     &              0.5d0*(df(i)/dpsi(i))*(dsk(i)/alp33k(i)) *
     &              (alfa33(i-1)/ds(i-1))

            b(i,1) = -0.5d0*dfidt(i) + ( 0.5d0*df(i)-Fk(i) ) *
     &                       ( alfa22(i)/ds(i) )/sigma(i) -
     &   	      0.5d0*dfidt(i)*capsi 

            b(i,2) = 0.5d0*dpsidt(i) + ( rIk(i)-0.5d0*dri(i) )*
     &                        ( alfa33(i)/ds(i) )/sigma(i) +
     &   	     0.5d0*dpsidt(i)*capfi 

            b(i,3) = 0.5d0*dri(i)/dpsi(i) + alfa22(i)/ds(i)

            b(i,4) = Qk(i)*alfa33(i)/ds(i) - 0.5d0*(df(i)/dpsi(i)) *
     &               (dsk(i)/alp33k(i))*(alfa33(i)/ds(i))

            c(i,1) = -dfi(i)/dt - ( Fk(i) * 
     &                (alfa22(i)/ds(i) + alfa22(i-1)/ds(i-1)) -
     &                0.5d0*df(i) *
     &                (alfa22(i)/ds(i) - alfa22(i-1)/ds(i-1))) /
     &                sigma(i) - dfidt(i)*capsi

            c(i,2) = dpsi(i)/dt + ( 0.5d0*dri(i)*
     &	           (alfa33(i-1)/ds(i-1) - alfa33(i)/ds(i)) +
     &             rIk(i)*(alfa33(i)/ds(i) + alfa33(i-1)/ds(i-1)) )
     &             /sigma(i) + dpsidt(i)*capfi

            c(i,3) = alfa22(i)/ds(i) + alfa22(i-1)/ds(i-1)

            c(i,4) = Qk(i)*( alfa33(i)/ds(i)+alfa33(i-1)/ds(i-1) ) +
     &             0.5d0*(df(i)/dpsi(i))*(dsk(i)/alp33k(i)) *
     &             ( alfa33(i-1)/ds(i-1)-alfa33(i)/ds(i) )

            Fk0  = 0.5d0*(f0(i)+f0(i-1))
            rIk0 = 0.5d0*(rI0(i)+rI0(i-1))
            Qk0  = 0.5d0*(q0(i)+q0(i-1))
            dri0 = ri0(i)-ri0(i-1)
            df0  = f0(i)-f0(i-1)  

            h(i,1)= dpsidt(i)*dfi(i) - dfidt(i)*dpsi(i) +
     &	            dpsidt(i)*capfi*d2fi - dfidt(i)*capsi*d2psi +
     &              ( cbut_b(i)*dvk(i) )/sigma(i) -
     &	            (f(i)*f(i-1))*( ri(i)/f(i)-ri(i-1)/f(i-1) )/sigma(i)

            h(i,2)= dri(i)-dfdpsi(i)*dsk(i)/alp33k(i)-dpdpsi(i)*dvk(i)
         enddo

! boundary equations
         b(1,1) = rm*ac0/(sigma(1)*skcen)
         b(1,2) = 0.d0
         b(1,3) = 0.d0
         b(1,4) = 0.d0

         c(1,1) = rm*ac0/(sigma(1)*skcen)-1.d0/dt
         c(1,2) = 0.d0
         c(1,3) = 0.d0
         c(1,4) = 1.d0

         h(1,1)= dpsidt(1)-rm*ac0*(psi(1)-psi(2))/(sigma(1)*skcen)
     &                  +cbut_b(1)*rm**2/f(1)/sigma(1)
         h(1,2)= 0.d0

         if(key_boncon.eq.0) then
            a(n,1) = cps_bon
            a(n,2) = 0.d0
            a(n,3) = 0.d0
            a(n,4) = 1.d0

            c(n,1) = cps_bon + 1.d0
            c(n,2) = 0.d0
            c(n,3) = 0.d0
            c(n,4) = 1.d0

            h(n,1) = psi_ext+cps_bon*psi(n-1)-(cps_bon+1.d0)*psi(n) +
     &              dps_bon*dpdpsi(n)
            h(n,2) = fbnd*ds(n-1)/alfa33(n-1)-(fi(n)-fi(n-1))
         else
            psi0n=psi0(n)
            psin=psi(n)
            psin1=psi(n-1)
            fi0n=fi0(n)
            fin=fi(n)
            fin1=fi(n-1)
            fn=fbnd
            fn1=f(n-1)
            rin=ri(n)
            rin1=ri(n-1)
            al22n=alfa22(n)
            al22n1=alfa22(n-1)
            al22x=1.d1
            al33n=alp33k(n)
            al33n1=alfa33(n-1)
            dsn=dsk(n)
            dsn1=ds(n-1)
            dtsig=-dt/sigma(n)/dsn
            av=cps_bon
            bv=-bps_bon
            cv=-dps_bon
            psiext=psi_ext
            pdv=dpdpsi(n)*(psin-psin1)*dvk(n)/dsk(n)
            pbb=dpdpsi(n)*(psin-psin1)*cv
            cbutn=cbut_b(n)*dvk(n)

      a11= -(2*(rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n)/
     &  al33n/(-fn1*dtsig*al22n+fin-fi0n)**2*al22n**2*fn*dtsig-2*rin1)/
     &  al22n*al22n1/dsn1

      a12= -(2*(rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n)*
     &  *2/al33n**2/(-fn1*dtsig*al22n+fin-fi0n)**3*al22n**2*dtsig+2*fn1/
     &  al33n)*al33n1/dsn1

      c11= -(2*(rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n)/
     &  al33n/(-fn1*dtsig*al22n+fin-fi0n)**2*al22n**2*fn*dtsig-2*rin1)/
     &  al22n*al22n1/dsn1+2*(rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+
     &  cbutn*dtsig*al33n)/al33n**2/(-fn1*dtsig*al22n+fin-fi0n)**2*
     &  al22n*fn

      c12= -(2*(rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n
     &  )**2/al33n**2/(-fn1*dtsig*al22n+fin-fi0n)**3*al22n**2*dtsig+
     &  2*fn1/al33n)*al33n1/dsn1+2*
     &  (rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n)**2/
     &  al33n**2/(-fn1*dtsig*al22n+fin-fi0n)**3*al22n

      a21= -2*av*(psin-psin1)-psin+psiext
      a22= 2*bv*fn1*al33n1/dsn1

      c21= -(1+av)*(psin-psin1)-psin+psiext-av*(psin-psin1)
      c22= 2*bv*fn1*al33n1/dsn1

      h1= ((rin1*fn*dtsig*al33n-fn*psin+fn*psi0n+cbutn*dtsig*al33n)**2/
     &  al33n**2/(-fn1*dtsig*al22n+fin-fi0n)**2*al22n**2-rin1**2)/al22n-
     &  (fn**2-fn1**2)/al33n-pdv

      h2= (psin-psiext+av*(psin-psin1))*(psin-psin1)+
     &  bv*(fn**2-fn1**2)+pbb

      a(n,1)= a11
      a(n,2)= a12
      a(n,3)= a21  ! 0.d0
      a(n,4)= a22  ! 0.d0 

      c(n,1)= c11
      c(n,2)= c12
      c(n,3)= c21  ! -1.d0
      c(n,4)= c22  !  0.d0

      h(n,1)=h1
      h(n,2)=h2   !0.d0

         endif

         call mtrx_prog (n,y,a,b,c,h)

         do i=1,n
            psi(i)= psi(i) + y(i,1)   ! y(i,1) = delta(psi(i))
            fi(i) = fi(i)  + y(i,2)   ! y(i,2) = delta(fi (i))
            dpsidt(i)=(psi(i)-psi0(i))/dt
            dfidt(i )=(fi(i)-fi0(i))/dt
         enddo

         do i=1,n-1
            f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)   
            ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)   
            q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))   
         enddo

         ffprim=rm*(ac0*(psi(1)-psi(2))/skcen-rm*dpdpsi(1))
         delf(1) = ffprim - dfdpsi(1) 
         dfdpsi(1)=ffprim
         do i=2,n-1
            ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))   
            delf(i) = ffprim - dfdpsi(i) 
            dfdpsi(i) =  ffprim  
         enddo
         alfa22n=alfa22(n)
         F(N)=fbnd 
         rI(n)=dsqrt( ri(n-1)**2+alfa22n/alp33k(n)*(F(N)**2-f(n-1)**2)
     &      +alfa22n*dpdpsi(n)*(psi(n)-psi(n-1))*dvk(n)/dsk(n) )
         dfdpsi(n) =  (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))

         deltamax = 0.d0
         do k=1,2
            do i=1,n
               aby = dabs(y(i,k))
               if (aby.ge.deltamax) deltamax=aby
            enddo
         enddo

         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
               capfi= cap_fi 
            endif

            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif
            d2psi= (psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi= (fi(i+1)-2.d0*fi(i)+fi(i-1))        

            znv = (dpsidt(i)*(fi(i+1) - fi(i-1)) -
     &             dfidt(i)*(psi(i+1) - psi(i-1))) -
     &             (f(i)*f(i-1))*(ri(i)/f(i) - ri(i-1)/f(i-1) ) * 
     &             2.d0/sigma(i) + dpsidt(i)*capfi*d2fi - 
     &             dfidt(i)*capsi*d2psi

            zn(i,1) = 0.5d0*znv

            zn(i,2) = (ri(i)-ri(i-1)) - dfdpsi(i)*dsk(i)/alp33k(i) -
     &                 dpdpsi(i)*dvk(i)
         enddo

         ddivmax = 0.d0
         do k=1,2
            do i=1,n
               abd = dabs(zn(i,k))
               if (abd.ge.ddivmax) ddivmax=abd
            enddo
         enddo

         cur_mu=ri(n)

         if (deltamax.le.epsel) exit

         if(ks.eq.199) then
            write(*,*) 'no newton iterations convergence in promat' 
            write(*,*) 'execution was terminated'
            write(*,*) 'ks=',ks
            stop
         endif	
             	    
      enddo
!---------------------------------------------------------------------

      k_step_prev=k_step
      psi_bn1=psi(n)

      f_wght=0.75d0
      do i=1,n
         dfdpsi(i)=f_wght*dfdpsi(i)+(1.d0-f_wght)*dfdpsn(i)        
      enddo

! Self-consistency check

      do i=1,n
         psi_d(i)=psi(i)-psi(n)         
         fi_d(i)=fi(i)         
         f_d(i)=f(i)         
         ri_d(i)=ri(i)         
         ps_pnt(i)=dpsidt(i)         
      enddo

      return
      end subroutine promat

!----------------------------------------------------------------
      subroutine promat_I(n, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm, 
     &                    psim, iter, k_step, fvac, psi_eav, pspl_av,
     &                    psibon0, tok)
!      diffusion & equilibrium both equations solvinq 
!      Newton metod with linearisation
! error-correcting 

      use sp_parameters, only: nrp
      use keys, only: kastr, key_fixfree

      implicit none

      integer, intent(in) :: n, iter, k_step
      real*8, intent(in) :: rm, psim, fvac, pspl_av, tok
      real*8, intent(in), dimension(n) :: flx_fi, psia, dpdpsi
      real*8, intent(out) :: psibon0, psi_eav
      real*8, intent(out), dimension(n) :: dfdpsi, f, q

      integer :: i, k, ks, key_boncon, k_step_prev, jumpstep
      real*8 :: del_psb, psi_bn1, rm0, ac0n, skcen0, xa, ffprim,
     &   alfa22n, epsel, dt, fbnd, skcen, dpsdt, d2psi, d2fi,
     &   ac0, fk0, rik0, qk0, dri0, df0, psi0n, fi0n, fin, fin1,
     &   cap_psi, cap_fi, capfi, capsi, pdv, pbb, cbutn,
     &   pspl_tst, psi_ext, psiext, psi_bnd, psin, psin1,
     &   bps_bon, cps_bon, dps_bon, av, bv, cv,
     &   a11, a12, a21, a22, c11, c12, c21, c22, h1, h2,
     &   ds1, dtsig, deltamax, ddivmax, abd, aby, znv, f_wght,
     &   dfldt, fn, fn1, rin, rin1, dsn, dsn1,
     &   al22n, al22n1, al22x, al33n, al33n1
      real*8, dimension(512) :: wrk
      real*8, dimension(nrp) :: psi, fi, ri, alfa22, alfa33, alp33k, ds,
     &   dsk, dv, dvk, dpsi, dfi, dri, df, fk, rik, qk, delf, dfdpsn, 
     &   psi_d, fi_d, f_d, ri_d, ps_pnt,
     &   psi0, fi0, f0, ri0, q0, dpsidt, dfidt,
     &   sigma, cbut_b
      real*8, dimension(nrp, 2) :: y, h, zn
      real*8, dimension(nrp, 4) :: a, b, c
      real*8, external :: fun_jb, fun_jb_a, fun_sig, fun_sig_a

      common /selcon/ psi_d, fi_d, f_d, ri_d, ps_pnt, del_psb, psi_bn1
      common /savt0/ psi0, fi0, f0, ri0, q0, dpsidt, dfidt, 
     &               rm0, ac0n, skcen0
      common /com_but/ sigma, cbut_b

      save k_step_prev

      epsel = 1.0d-10
      dt = 1.d0
      cap_psi=0.d0   !1.d0
      cap_fi =0.d0
      fbnd = fvac

      if(k_step.eq.1)then
         psi_bnd=psibon0
      else
         psi_bnd=psi_bn1
      endif

      psi_eav=psi_bnd-pspl_av
         	  
! metric coefficients
      call metcof(alfa22, alfa33, ds, dv, ac0, skcen, ri,
     &            alp33k, dsk, dvk)

      do i=1,n
         dfdpsn(i)=dfdpsi(i)
         psi(i)=psim*psia(i)+psi_bnd
         fi(i)=flx_fi(i)
         if(k_step.ne.k_step_prev) then
            jumpstep=k_step - k_step_prev
            del_psb=psi_bn1-psibon0
            if((jumpstep.eq.1 .AnD. k_step.gt.1) .OR. 
     &         key_fixfree.eq.0) then
               psibon0=psi_bn1
            endif
            psi0(i)=psim*psia(i)+psibon0
            fi0(i)=fi(i)
            F0(i) =f(i)
            rI0(i) =rI(i)
            Q0(i) =q(i)
            dpsidt(i)=0.d0
            dfidt(i )=0.d0
         endif

         xa=dsqrt(flx_fi(i)/flx_fi(n))	  	  
         if(kastr.eq.1) then
            sigma(i) =fun_sig(xa,i)  
            cbut_b(i) =fun_jb(xa,i)  
         else
            sigma(i) =fun_sig_a(xa,i)  
            cbut_b(i) =fun_jb_a(xa,i)  
         endif

      enddo

      if(iter.eq.1) then
         rm0=rm
         skcen0=skcen
         ac0n=ac0
      endif

      do i=1,n
         wrk(i)=dsqrt(dabs(alfa22(i)))
      enddo
      do i=1,n-1
         fi(i)=fi(n)*(i-1.d0)**2/(n-1.d0)**2
      enddo
      do i=1,n-1
         f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)
         ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)
         q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))
      enddo

      F(N)=fbnd
      rI(n)=tok

      ffprim=rm*(ac0*(psi(1)-psi(2))/skcen-rm*dpdpsi(1))
      delf(1) = ffprim - dfdpsi(1)
      dfdpsi(1)=ffprim
      do i=2,n-1
         ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))
         delf(i) = ffprim - dfdpsi(i)
         dfdpsi(i) =  ffprim
         dfdpsn(i) =alp33k(i)*(ri(i)-ri(i-1))/dsk(i)
      enddo
          
      alfa22n=alfa22(n)
      dfdpsi(n) =  (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))

      do ks = 1,100      ! iteration cycle
         do i=2,n-1
            dpsi(i)=0.5d0*(psi(i+1)-psi(i-1))
            dfi(i) =0.5d0*(fi(i+1)-fi(i-1))
            dri(i) =ri(i)-ri(i-1)
            df(i)  =f(i)-f(i-1)  
            Fk(i) =0.5d0*(f(i)+f(i-1))
            rIk(i) =0.5d0*(rI(i)+rI(i-1))
            Qk(i) =-(dsk(i)/alp33k(i))*(f(i)+f(i-1))/(psi(i+1)-psi(i-1))
         enddo

! matrix completing:   a(i)*y(i-1)-c(i)*y(i)+b(i)*y(i+1)=-h(i)
         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
               capfi= cap_fi 
            endif

            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif

            d2psi=0.5d0*(psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi=0.5d0*(fi(i+1)-2.d0*fi(i)+fi(i-1))        

            a(i,1) = 0.5d0*dfidt(i)-( Fk(i)+0.5d0*df(i) ) *
     &               ( alfa22(i-1)/ds(i-1) )/sigma(i)     -
     &               0.5d0*dfidt(i)*capsi 

            a(i,2) = -0.5d0*dpsidt(i)+( rIk(i)+0.5d0*dri(i) ) *
     &	                   ( alfa33(i-1)/ds(i-1) )/sigma(i)   +
     &   	      0.5d0*dpsidt(i)*capfi 

            a(i,3) = -0.5d0*dri(i)/dpsi(i) + alfa22(i-1)/ds(i-1)

            a(i,4) = Qk(i)*alfa33(i-1)/ds(i-1)                +
     &               0.5d0*(df(i)/dpsi(i))*(dsk(i)/alp33k(i)) *
     &               (alfa33(i-1)/ds(i-1))

            b(i,1) = -0.5d0*dfidt(i)+( 0.5d0*df(i)-Fk(i) ) *
     &                ( alfa22(i)/ds(i) )/sigma(i)         -
     &   	      0.5d0*dfidt(i)*capsi 

            b(i,2) = 0.5d0*dpsidt(i)+( rIk(i)-0.5d0*dri(i) ) *
     &               ( alfa33(i)/ds(i) )/sigma(i)            +
     &   	     0.5d0*dpsidt(i)*capfi 

            b(i,3) = 0.5d0*dri(i)/dpsi(i)+alfa22(i)/ds(i)

            b(i,4) = Qk(i)*alfa33(i)/ds(i) - 0.5d0*(df(i)/dpsi(i)) *
     &               (dsk(i)/alp33k(i))*(alfa33(i)/ds(i))

            c(i,1) = -dfi(i)/dt - ( Fk(i) *
     &               (alfa22(i)/ds(i) + alfa22(i-1)/ds(i-1)) -
     &               0.5d0*df(i) *
     &               (alfa22(i)/ds(i)-alfa22(i-1)/ds(i-1)))/sigma(i) -
     &               dfidt(i)*capsi

            c(i,2) = dpsi(i)/dt +
     &	      ( 0.5d0*dri(i)*(alfa33(i-1)/ds(i-1)-alfa33(i)/ds(i)) +
     &        rIk(i)*(alfa33(i)/ds(i)+alfa33(i-1)/ds(i-1)) )/sigma(i) +
     &        dpsidt(i)*capfi

            c(i,3) = alfa22(i)/ds(i) + alfa22(i-1)/ds(i-1)

            c(i,4)= Qk(i)*( alfa33(i)/ds(i)+alfa33(i-1)/ds(i-1) ) +
     &              0.5d0*(df(i)/dpsi(i))*(dsk(i)/alp33k(i)) *
     &              ( alfa33(i-1)/ds(i-1)-alfa33(i)/ds(i) )

            Fk0 =0.5d0*(f0(i)+f0(i-1))
            rIk0 =0.5d0*(rI0(i)+rI0(i-1))
            Qk0 =0.5d0*(q0(i)+q0(i-1))
            dri0 =ri0(i)-ri0(i-1)
            df0  =f0(i)-f0(i-1)  

            h(i,1)= dpsidt(i)*dfi(i)     - dfidt(i)*dpsi(i) +
     &	            dpsidt(i)*capfi*d2fi - dfidt(i)*capsi*d2psi +
     &              ( cbut_b(i)*dvk(i) )/sigma(i) -
     &	          (f(i)*f(i-1))*( ri(i)/f(i)-ri(i-1)/f(i-1) )/sigma(i)

            h(i,2)= dri(i)-dfdpsi(i)*dsk(i)/alp33k(i)-dpdpsi(i)*dvk(i)
         enddo

! boundary equations
         b(1,1) = rm*ac0/(sigma(1)*skcen)
         b(1,2) = 0.d0
         b(1,3) = 0.d0
         b(1,4) = 0.d0

         c(1,1) = rm*ac0/(sigma(1)*skcen)-1.d0/dt
         c(1,2) = 0.d0
         c(1,3) = 0.d0
         c(1,4) = 1.d0

         h(1,1) = dpsidt(1)-rm*ac0*(psi(1)-psi(2))/(sigma(1)*skcen) +
     &            cbut_b(1)*rm**2/f(1)/sigma(1)
         h(1,2) = 0.d0

         psi0n=psi0(n)
         psin=psi(n)
         psin1=psi(n-1)
         fi0n=fi0(n)
         fin=fi(n)
         fin1=fi(n-1)
         fn=fbnd
         fn1=f(n-1)
         rin=ri(n)
         rin1=ri(n-1)
         al22n=alfa22(n)
         al22n1=alfa22(n-1)
         al22x=1.d1
         al33n=alp33k(n)
         al33n1=alfa33(n-1)
         dsn=dsk(n)
         dsn1=ds(n-1)
         dtsig=-dt/sigma(n)/dsn

         pdv=dpdpsi(n)*(psin-psin1)*dvk(n)/dsk(n)
         cbutn=cbut_b(n)*dvk(n)

         a11= -fn*al22n1/dsn1
         a12= rin*al33n1/dsn1

         c11= -fn*al22n1/dsn1
         c12= rin*al33n1/dsn1

         a21= 2*rin1/al22n*al22n1/dsn1
         a22= -2*fn1/al33n*al33n1/dsn1

         c21= 2*rin1/al22n*al22n1/dsn1
         c22= -2*fn1/al33n*al33n1/dsn1

         h1= (fin-fi0n)/dtsig*rin/al22n-(psin-psi0n)/dtsig*fn/al33n +
     &       cbutn - rin*fn1 + rin1*fn

         h2= (rin**2-rin1**2)/al22n - (fn**2-fn1**2)/al33n - pdv

         a(n,1)= a11
         a(n,2)= a12
         a(n,3)= a21
         a(n,4)= a22

         c(n,1)= c11
         c(n,2)= c12
         c(n,3)= c21
         c(n,4)= c22

         h(n,1)=h1
         h(n,2)=h2

         call mtrx_prog(n, y, a, b, c, h)

         do i=1,n
            psi(i)= psi(i) + y(i,1)   ! y(i,1) = delta(psi(i))
            fi(i) = fi(i)  + y(i,2)   ! y(i,2) = delta(fi (i))
            dpsidt(i)=(psi(i)-psi0(i))/dt
            dfidt(i )=(fi(i)-fi0(i))/dt
         enddo

         do i=1,n-1
            f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)   
            ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)   
            q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))   
         enddo

         ffprim=rm*(ac0*(psi(1)-psi(2))/skcen-rm*dpdpsi(1))
         delf(1) = ffprim - dfdpsi(1) 
         dfdpsi(1)=ffprim
         do i=2,n-1
            ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))   
            delf(i) = ffprim - dfdpsi(i) 
            dfdpsi(i) =  ffprim  
            dfdpsn(i) =alp33k(i)*(ri(i)-ri(i-1))/dsk(i) 
         enddo

         alfa22n=alfa22(n)
         dfdpsi(n) =  (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))
         deltamax = 0.d0
         do k=1,2
            do i=1,n
               aby = dabs(y(i,k))
               if (aby.ge.deltamax) deltamax=aby
            enddo
         enddo

         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
               capfi= cap_fi 
            endif
            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif

            d2psi= (psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi= (fi(i+1)-2.d0*fi(i)+fi(i-1))        

            znv = (dpsidt(i)*(fi(i+1) - fi(i-1))   - 
     &             dfidt(i)*(psi(i+1) - psi(i-1))) -
     &            (f(i)*f(i-1))*( ri(i)/f(i)-ri(i-1)/f(i-1) )*2.d0 /
     &            sigma(i)    +
     &	          dpsidt(i)*capfi*d2fi - dfidt(i)*capsi*d2psi

            zn(i,1)=0.5d0*znv

            zn(i,2)= (ri(i)-ri(i-1)) - dfdpsi(i)*dsk(i)/alp33k(i) -
     &               dpdpsi(i)*dvk(i)
         enddo

         ddivmax = 0.d0
         do k=1,2
            do i=1,n
               abd = dabs(zn(i,k))
               if (abd.ge.ddivmax) ddivmax=abd
            enddo
         enddo

         if (deltamax.le.epsel) EXIT

         if(ks.eq.199) then
            write(*,*) 'no newton iterations convergence in promat' 
            write(*,*) 'execution was terminated'
            write(*,*) 'ks=',ks
            stop
         endif	
             	    
      enddo                     ! iteration cycle ends

!---------------------------------------------------------------------

      k_step_prev=k_step
      psi_bn1=psi(n)

! Self-consistency check

      do i=1,n
         psi_d(i)=psi(i)-psi(n)         
         fi_d(i)=fi(i)         
         f_d(i)=f(i)         
         ri_d(i)=ri(i)         
         ps_pnt(i)=dpsidt(i)         
      enddo

      return
      end subroutine promat_I

!----------------------------------------------------------------
      subroutine mtrx_prog (n, y, a, b, c, f)
! matrix progonka for two equations

      use sp_parameters, only: nrp

      implicit none

      integer, intent(in) :: n
      real*8, intent(in), dimension(nrp, 4) :: a, b, c 
      real*8, intent(in), dimension(nrp, 2) :: f
      real*8, intent(out), dimension(nrp, 2) :: y

      integer :: j
      real*8 :: det, c1i1, c1i2, c1i3, c1i4, fab1, fab2,
     &   cmb1, cmb2, cmb3, cmb4, cmbi1, cmbi2, cmbi3, cmbi4
      real*8, dimension(nrp, 4) :: alfa
      real*8, dimension(nrp, 2) :: beta

!------------------------------------------------------------------
!    a(i)*y(i-1)-c(i)*y(i)+b(i)*y(i+1)=-f(i), i=2,n-1
!   -c(1)*y(1)+b(1)*y(2)=-f(1), a(n)*y(n-1)-c(n)*y(n)=-f(n)
!------------------------------------------------------------------
! a(i), b(i), c(i) -> matrix [2*2]; y(i), f(i) -> vector [2]
! {a(i)11  a(i)12  a(i)21  a(i)22} -> {a(i)1  a(i)2  a(i)3  a(i)4}
!------------------------------------------------------------------
      det=c(1,1)*c(1,4)-c(1,2)*c(1,3)   ! inverse c(1)
      c1i1= c(1,4)/det
      c1i2=-c(1,2)/det
      c1i3=-c(1,3)/det
      c1i4= c(1,1)/det

      alfa(1,1) = c1i1*b(1,1)+c1i2*b(1,3)
      alfa(1,2) = c1i1*b(1,2)+c1i2*b(1,4)
      alfa(1,3) = c1i3*b(1,1)+c1i4*b(1,3)
      alfa(1,4) = c1i3*b(1,2)+c1i4*b(1,4)

      beta(1,1) = c1i1*f(1,1)+c1i2*f(1,2)
      beta(1,2) = c1i3*f(1,1)+c1i4*f(1,2)

      do j=2,n
         cmb1 = c(j,1) - a(j,1)*alfa(j-1,1)-a(j,2)*alfa(j-1,3)
         cmb2 = c(j,2) - a(j,1)*alfa(j-1,2)-a(j,2)*alfa(j-1,4)
         cmb3 = c(j,3) - a(j,3)*alfa(j-1,1)-a(j,4)*alfa(j-1,3)
         cmb4 = c(j,4) - a(j,3)*alfa(j-1,2)-a(j,4)*alfa(j-1,4)

         det=cmb1*cmb4-cmb2*cmb3   !  inverse {c(i)-a(i)*alfa(i-1)}
         cmbi1= cmb4/det
         cmbi2=-cmb2/det
         cmbi3=-cmb3/det
         cmbi4= cmb1/det

         alfa(j,1) = cmbi1*b(j,1)+cmbi2*b(j,3)
         alfa(j,2) = cmbi1*b(j,2)+cmbi2*b(j,4)
         alfa(j,3) = cmbi3*b(j,1)+cmbi4*b(j,3)
         alfa(j,4) = cmbi3*b(j,2)+cmbi4*b(j,4)

         fab1 = f(j,1) + a(j,1)*beta(j-1,1)+a(j,2)*beta(j-1,2)
         fab2 = f(j,2) + a(j,3)*beta(j-1,1)+a(j,4)*beta(j-1,2)

         beta(j,1) = cmbi1*fab1+cmbi2*fab2
         beta(j,2) = cmbi3*fab1+cmbi4*fab2
      enddo

      y(n,1) = beta(n,1)
      y(n,2) = beta(n,2)

      do j=n-1,1,-1
         y(j,1) = alfa(j,1)*y(j+1,1)+alfa(j,2)*y(j+1,2) + beta(j,1)
         y(j,2) = alfa(j,3)*y(j+1,1)+alfa(j,4)*y(j+1,2) + beta(j,2)
      enddo

      return
      end subroutine mtrx_prog

!----------------------------------------------------------------
      real*8 function fun_sig(a, i)

      use sp_parameters, only: amu0
      use sigcd, only: C_sig
      use tim, only: dtim, ctim

      implicit none

      integer, intent(in) :: i
      real*8, intent(in) :: a

      fun_sig=-C_sig(i)*amu0/dtim

      return
      end function fun_sig

!----------------------------------------------------------------
      real*8 function fun_sig_a(a, i)

      use sp_parameters, only: pi
      use tim, only: dtim

      implicit none

      integer, intent(in) :: i
      real*8, intent(in) :: a
      real*8 :: a0, a1, t, t0, t1, fun_sig

      a0=0.0d0
      a1=1.0d0
      t0 =1.0d3
      t1 =1.0d2

      if(a.le.a0) then
         t= t0
      elseif(a.gt.a0 .AND. a.lt.a1) then
         t=t1 + (t0-t1)*(1.d0-((a-a0)/(a1-a0))**2)
      elseif(a.ge.a1) then
         t=t1
      endif

      fun_sig=-(8.d0*pi/3.d0)*dabs(t)**1.5d0
      fun_sig_a=fun_sig*(1.d-4/dtim)

!temporary
      if(a.ge.0.5d0) fun_sig_a=fun_sig_a*0.25d0

      return
      end function fun_sig_a

!----------------------------------------------------------------
      real*8 function fun_jb(a, i)

      use sp_parameters, only: nrp, amu0
      use sigcd, only: C_bts, c_driv

      implicit none

      integer, intent(in) :: i
      real*8, intent(in) :: a

      fun_jb=(C_bts(i)+C_driv(i))*amu0

      return
      end function fun_jb

!----------------------------------------------------------------
      real*8 function fun_jb_a(a, i)

      use sp_parameters, only: nrp, pi
      use jb, only: Bj_av

      implicit none

      integer, intent(in) :: i
      real*8, intent(in) :: a

      real*8 :: ampl, a0, w, fun

      ampl=0.0d0
      a0=0.50d0
      w=0.05d0
      fun=BJ_av(1)*( 1.d0-dtanh(((a0-a)/w)**2) )
      fun_jb_a=ampl*fun

      return
      end function fun_jb_a

!----------------------------------------------------------------
      subroutine promat_j(n, dfdpsi, dpdpsi, f, flx_fi, psia, q, rm,
     &                    psim, iter, k_step, fvac, psi_eav, pspl_av,
     &                    psibon0, cur_mu)
! diffusion & equilibrium both equations solvinq 
! Newton metod with linearisation

      use sp_parameters, only: nrp
      use jb, only: Bj_av

      implicit none

      integer, intent(in) :: n, iter, k_step
      real*8, intent(in) :: psim, fvac, psi_eav, pspl_av, psibon0
      real*8, intent(in), dimension(n) :: dpdpsi, flx_fi, psia
      real*8, intent(out) :: cur_mu
      real*8, intent(out), dimension(n) :: dfdpsi, f, q

      integer :: i, k, ks, k_step_prev
      real*8 :: abd, aby, ac0n, skcen0, epsel, dt, ddivmax, deltamax,
     &   cap_psi, cap_fi, capfi, capsi, cbutn,
     &   bps_bon, cps_bon, dps_bon, ac0, skcen, xa, av, bv, pdv,
     &   pspl_tst, psi_bnd, psin, psi0n, psin1, psiext, psi_ext,
     &   fi0n, fin, fin1, fn, fn1, fbnd,
     &   rin, rin1, rm0, rm, ffprim, dpsdt, dfldt, dsn, dsn1, dtsig,
     &   d2fi, d2psi, fk0, rik0, qk0, dri0, df0, znv,
     &   al22n, al22n1, alfa22n, al33n, al33n1,
     &   a11, a12, c11, c12, h1
      real*8, dimension(512) :: wrk
      real*8, dimension(nrp, 2) :: y, h, zn
      real*8, dimension(nrp, 4) :: a, b, c
      real*8, dimension(nrp) :: psi, fi, ri, alfa22, alfa33, alp33k, 
     &   ds, dsk, dv, dvk, dpsi, dfi, dri, df, Fk, rIk, Qk, delf, dfdpsn
      real*8, dimension(nrp) :: sigma, cbut_b,
     &   psi0, fi0, f0, ri0, q0, dpsidt, dfidt

      common /savt0/ psi0, fi0, f0, ri0, q0, 
     &               dpsidt, dfidt, rm0,ac0n,skcen0
      common /com_but/ sigma, cbut_b

      save k_step_prev

      epsel = 1.0d-11
      dt = 1.d20
      cap_psi=0.d0   !1.d0
      cap_fi =0.d0
      fbnd = fvac

      call cof_bon(cps_bon, bps_bon, dps_bon)

      pspl_tst = -cps_bon*(psia(n) - psia(n-1))*psim +
     &            bps_bon*dfdpsi(n)+dps_bon*dpdpsi(n)
      psi_bnd=0.d0

! metric coefficients
      call metcof(alfa22, alfa33, ds, dv, ac0, skcen, ri, alp33k, 
     &            dsk, dvk)

      do i=1,n
         dfdpsn(i)=dfdpsi(i)
         psi(i)=psim*psia(i)+psi_bnd
         fi(i)=flx_fi(i)
         if(k_step.ne.k_step_prev) then
            psi0(i)=psim*psia(i)+psibon0
            fi0(i)=fi(i)
            F0(i) =f(i)
            rI0(i) =rI(i)
            Q0(i) =q(i)
            dpsidt(i)=(psi(i)-psi0(i))/dt
            dfidt(i )=0.d0
         endif
         xa=dfloat(i-1)/dfloat(n-1)
         sigma(i) =1.d0
         cbut_b(i) =BJ_av(i)
      enddo

      if(iter.eq.1) then
         rm0=rm
         skcen0=skcen
         ac0n=ac0
      endif

      do i=1,n
         wrk(i)=dsqrt(dabs(alfa22(i)))
      enddo
      do i=1,n-1
         fi(i)=fi(n)*(i-1.d0)**2/(n-1.d0)**2
      enddo
      do i=1,n-1
         f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)
         ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)
         q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))
      enddo

      F(N)=fbnd
      ffprim=rm*(ac0*(psi(1)-psi(2))/skcen-rm*dpdpsi(1))
      delf(1) = ffprim - dfdpsi(1)
      dfdpsi(1)=ffprim
      do i=2,n-1
         ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))
         delf(i) = ffprim - dfdpsi(i)
         dfdpsi(i) =  ffprim
         dfdpsn(i) =alp33k(i)*(ri(i)-ri(i-1))/dsk(i)
      enddo
   
      alfa22n=alfa22(n)
      dfdpsi(n) = (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))


      do ks = 1,100 
               
         do i=2,n-1
            dpsi(i)= 0.5d0*(psi(i+1)-psi(i-1))
            dfi(i) = 0.5d0*(fi(i+1)-fi(i-1))
            dri(i) = ri(i)-ri(i-1)
            df(i)  = f(i)-f(i-1)  
            Fk(i)  = 0.5d0*(f(i)+f(i-1))
            rIk(i) = 0.5d0*(rI(i)+rI(i-1))
            Qk(i)  = -(dsk(i)/alp33k(i)) * 
     &                (f(i)+f(i-1))/(psi(i+1)-psi(i-1))
         enddo

! matrix completing:   a(i)*y(i-1)-c(i)*y(i)+b(i)*y(i+1)=-h(i)
         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
              capfi= cap_fi 
            endif

            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif

            d2psi=0.5d0*(psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi=0.5d0*(fi(i+1)-2.d0*fi(i)+fi(i-1))        

            a(i,1) = 0.5d0*dfidt(i) - ( Fk(i) + 0.5d0*df(i) ) *
     &   	     ( alfa22(i-1)/ds(i-1) )/sigma(i)         -
     &   	     0.5d0*dfidt(i)*capsi 

            a(i,2) = -0.5d0*dpsidt(i) + ( rIk(i) + 0.5d0*dri(i) ) *
     &	              ( alfa33(i-1)/ds(i-1) )/sigma(i)            +
     &   	      0.5d0*dpsidt(i)*capfi 

            a(i,3) = -0.5d0*dri(i)/dpsi(i) + alfa22(i-1)/ds(i-1)

            a(i,4) = Qk(i)*alfa33(i-1)/ds(i-1) + 0.5d0 *
     &        (df(i)/dpsi(i))*(dsk(i)/alp33k(i))*(alfa33(i-1)/ds(i-1))

            b(i,1) = -0.5d0*dfidt(i)+( 0.5d0*df(i)-Fk(i) )  *
     &                ( alfa22(i)/ds(i) )/sigma(i)          -
     &   	      0.5d0*dfidt(i)*capsi 

            b(i,2) = 0.5d0*dpsidt(i)+( rIk(i)-0.5d0*dri(i) ) *
     &                        ( alfa33(i)/ds(i) )/sigma(i)   +
     &   	     0.5d0*dpsidt(i)*capfi

            b(i,3) = 0.5d0*dri(i)/dpsi(i) + alfa22(i)/ds(i)

            b(i,4) = Qk(i)*alfa33(i)/ds(i) - 0.5d0 *
     &           (df(i)/dpsi(i))*(dsk(i)/alp33k(i))*(alfa33(i)/ds(i))

            c(i,1) = -dfi(i)/dt - ( Fk(i) *
     &               (alfa22(i)/ds(i) + alfa22(i-1)/ds(i-1))   -
     &               0.5d0*df(i)                               *
     &               ( alfa22(i)/ds(i)-alfa22(i-1)/ds(i-1) ) ) /
     &               sigma(i) - dfidt(i)*capsi

            c(i,2) = dpsi(i)/dt + 
     &	       ( 0.5d0*dri(i)*(alfa33(i-1)/ds(i-1)-alfa33(i)/ds(i)) +
     &         rIk(i)*(alfa33(i)/ds(i)+alfa33(i-1)/ds(i-1)) )/sigma(i) +
     &         dpsidt(i)*capfi

            c(i,3) = alfa22(i)/ds(i)+alfa22(i-1)/ds(i-1)

            c(i,4) = Qk(i)*( alfa33(i)/ds(i) + alfa33(i-1)/ds(i-1) ) +
     &               0.5d0*(df(i)/dpsi(i))*(dsk(i)/alp33k(i))        *
     &               ( alfa33(i-1)/ds(i-1) - alfa33(i)/ds(i) )

            Fk0  = 0.5d0*(f0(i)+f0(i-1))
            rIk0 = 0.5d0*(rI0(i)+rI0(i-1))
            Qk0  = 0.5d0*(q0(i)+q0(i-1))
            dri0 = ri0(i)-ri0(i-1)
            df0  = f0(i)-f0(i-1)  

            h(i,1) = dpsidt(i)*dfi(i)-dfidt(i)*dpsi(i) +
     &	             dpsidt(i)*capfi*d2fi-dfidt(i)*capsi*d2psi -
     &	         (f(i)*f(i-1))*( ri(i)/f(i)-ri(i-1)/f(i-1) )/sigma(i) +
     &	         ( cbut_b(i)*dvk(i) )/sigma(i)

            h(i,2)= dri(i)-dfdpsi(i)*dsk(i)/alp33k(i)-dpdpsi(i)*dvk(i)

         enddo

! boundary equations
         b(1,1) = rm*ac0/(sigma(1)*skcen)
         b(1,2) = 0.d0
         b(1,3) = 0.d0
         b(1,4) = 0.d0

         c(1,1) = rm*ac0/(sigma(1)*skcen)-1.d0/dt
         c(1,2) = 0.d0
         c(1,3) = 0.d0
         c(1,4) = 1.d0

         h(1,1) = dpsidt(1)-rm*ac0*(psi(1)-psi(2))/(sigma(1)*skcen)
     &                  +cbut_b(1)*rm**2/f(1)/sigma(1)
         h(1,2) = 0.d0

         psi0n=psi0(n)
         psin=psi(n)
         psin1=psi(n-1)
         fi0n=fi0(n)
         fin=fi(n)
         fin1=fi(n-1)
         fn=fbnd
         fn1=f(n-1)
         rin=ri(n)
         rin1=ri(n-1)
         al22n=alfa22(n)
         al22n1=alfa22(n-1)
         al33n=alp33k(n)
         al33n1=alfa33(n-1)
         dsn=dsk(n)
         dsn1=ds(n-1)
         dtsig=-dt/sigma(n)/dsn
         av=cps_bon
         bv=-bps_bon
         psiext=psi_ext
         pdv=dpdpsi(n)*(psin-psin1)*dvk(n)/dsk(n)
         cbutn=cbut_b(n)*dvk(n)

         a11 = -(2*(cbutn+rin1*fn)/fn1**2*fn-2*rin1)/al22n*al22n1/dsn1
         a12 = -(-2*(cbutn+rin1*fn)**2/fn1**3/al22n + 2*fn1/al33n) *
     &         al33n1/dsn1

         c11 = -(2*(cbutn+rin1*fn)/fn1**2*fn-2*rin1)/al22n*al22n1/dsn1
         c12 = -(-2*(cbutn+rin1*fn)**2/fn1**3/al22n+2*fn1/al33n) *
     &         al33n1/dsn1

         h1 = ((cbutn+rin1*fn)**2/fn1**2-rin1**2)/al22n - 
     &        (fn**2-fn1**2)/al33n - pdv

         a(n,1)= a11
         a(n,2)= a12
         a(n,3)= 0.d0
         a(n,4)= 0.d0 

         c(n,1)= c11
         c(n,2)= c12
         c(n,3)=-1.d0
         c(n,4)= 0.d0

         h(n,1)=h1
         h(n,2)=0.d0
! end include

         call mtrx_prog (n,y,a,b,c,h)

         do i=1,n
            psi(i)= psi(i) + y(i,1)   ! y(i,1) = delta(psi(i))
            fi(i) = fi(i)  + y(i,2)   ! y(i,2) = delta(fi (i))
            dpsidt(i)=(psi(i)-psi0(i))/dt
            dfidt(i )=(fi(i)-fi0(i))/dt
         enddo
         do i=1,n-1
            f(i)=alfa33(i)*(fi(i+1)-fi(i))/ds(i)   
            ri(i)=alfa22(i)*(psi(i+1)-psi(i))/ds(i)   
            q(i)=-(fi(i+1)-fi(i))/(psi(i+1)-psi(i))   
         enddo
         ffprim=rm*(ac0*(psi(1)-psi(2))/skcen-rm*dpdpsi(1))
         delf(1) = ffprim - dfdpsi(1) 
         dfdpsi(1)=ffprim
         do i=2,n-1
            ffprim = (f(i)**2-f(i-1)**2)/(psi(i+1)-psi(i-1))   
            delf(i) = ffprim - dfdpsi(i) 
            dfdpsi(i) =  ffprim  
            dfdpsn(i) =alp33k(i)*(ri(i)-ri(i-1))/dsk(i) 
         enddo
         alfa22n=alfa22(n)
         F(N)=fbnd 
         dfdpsi(n) =  (f(n)**2-f(n-1)**2)/(psi(n)-psi(n-1))
         deltamax = 0.d0
         do k=1,2
            do i=1,n
               aby = dabs(y(i,k))
               if (aby.ge.deltamax) deltamax=aby
            enddo
         enddo

         do i=2,n-1
            dpsdt=dpsidt(i)
            if(dpsdt.GE.0.d0) then
               capfi=-cap_fi 
            else
               capfi= cap_fi 
            endif
            dfldt=dfidt(i)
            if(dfldt.LE.0.d0) then
               capsi=-cap_psi 
            else
               capsi= cap_psi 
            endif
            d2psi= (psi(i+1)-2.d0*psi(i)+psi(i-1))        
            d2fi= (fi(i+1)-2.d0*fi(i)+fi(i-1))        

            znv = ( dpsidt(i)*(fi(i+1) - fi(i-1)) - 
     &              dfidt(i)*(psi(i+1) - psi(i-1)) ) -
     &             (f(i)*f(i-1))*( ri(i)/f(i) - ri(i-1)/f(i-1) )*2.d0 /
     &             sigma(i) + dpsidt(i)*capfi*d2fi - 
     &                         dfidt(i)*capsi*d2psi

            zn(i,1) = 0.5d0*znv

            zn(i,2) = ri(i) - ri(i-1) - dfdpsi(i)*dsk(i)/alp33k(i) -
     &                dpdpsi(i)*dvk(i)
         enddo

         ddivmax = 0.d0
         do k=1,2
            do i=1,n
               abd = dabs(zn(i,k))
               if (abd.ge.ddivmax) ddivmax=abd
            enddo
         enddo

         if (deltamax.le.epsel) EXIT

         rI(n) = dsqrt(ri(n-1)**2 + alfa22(n)/alp33k(n) * 
     &           (F(N)**2 - f(n-1)**2) + alfa22n*dpdpsi(n) *
     &           (psi(n)-psi(n-1))*dvk(n)/dsk(n) )

         cur_mu=ri(n)

         if(ks.eq.99) then
            write(*,*) 'no newton iterations convergence in promat' 
            write(*,*) 'execution was terminated'
            write(*,*) 'ks=',ks
            stop
         endif
             	    
      enddo

      k_step_prev=k_step

      return
      end subroutine promat_j

!----------------------------------------------------------------
      subroutine put_key_fix(k)

      use keys, only: key_fixfree

      implicit none

      integer, intent(in) :: k

      key_fixfree=k

      return
      end subroutine put_key_fix
