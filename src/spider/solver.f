      subroutine solve(isol, wdm)

      use sp_parameters
      use comblc, only: ni, ni1, nj, nj1, itin, r, z, curf

      implicit none

      integer, intent(in) :: isol
      real*8, intent(inout), dimension(nip, njp) :: wdm

      integer :: com_solver, i, j
      real*8 :: time_beg, time_end, wght
      real*8, dimension(nip) :: axm, ax0, axp
      real*8, dimension(nip, njp) :: wdmx, fwrk, ffwrk

      common /comffs/ axm, ax0, axp
      common /com_solver_ef/ com_solver

      do i=1,ni
         do j=1,nj
            wdmx(i,j)=wdm(i,j)
         enddo
      enddo

      do i=1,ni
         do j=1,nj
            fwrk(i,j)=-curf(i,j)*r(i)
         enddo
      enddo

      if (com_solver.eq.0) then
         call foursol(ni, nj, wdm, fwrk, ffwrk, r, z, axm, ax0, axp)
      endif
      if (com_solver.eq.1) then
       	 call foursol_ef(ni, nj, wdm, fwrk, ffwrk, r, z, axm, ax0, axp)
      endif
      if (com_solver.eq.2) then
       	 call foursol_ws(ni, nj, wdm, fwrk, ffwrk, r, z, axm, ax0, axp)
      endif

      wght=1.d0
      if(itin.le.2 ) wght=1.0d0
      do i=2,ni1
         do j=2,nj1
            wdm(i,j)=wdm(i,j)*wght+(1.d0-wght)*wdmx(i,j)
         enddo
      enddo

      return
      end subroutine solve

!----------------------------------------------------------------
      subroutine inverse_spid(a, c, n)
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! c(n,n) - inverse matrix of A

      implicit none 

      integer, intent(in) :: n
      double precision, intent(out), dimension(n, n) :: c
      double precision, intent(inout), dimension(n, n) :: a

      integer :: i, j, k
      double precision coeff
      double precision, dimension(n) :: b, d, x
      double precision, dimension(n, n) :: L, U

      L=0.0
      U=0.0
      b=0.0

! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            if (isnan(coeff)) write(*,*) k,coeff
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            enddo
         enddo
      enddo

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
         L(i,i) = 1.0
      enddo
! U matrix is the upper triangular part of A
      do j=1,n
         do i=1,j
            U(i,j) = a(i,j)
         enddo
      enddo

! Step 3: compute columns of the inverse matrix C
      do k=1,n
         b(k)=1.0
         d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
         do i=2,n
            d(i)=b(i)
            do j=1,i-1
               d(i) = d(i) - L(i,j)*d(j)
            enddo
         enddo
! Step 3b: Solve Ux=d using the back substitution
         x(n)=d(n)/U(n,n)
         do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
               x(i)=x(i)-U(i,j)*x(j)
            enddo
            x(i) = x(i)/u(i,i)
         enddo
! Step 3c: fill the solutions x(n) into column k of C
         do i=1,n
            c(i,k) = x(i)
         enddo
         b(k)=0.0
      enddo

      return
      end subroutine inverse_spid

!----------------------------------------------------------------
      subroutine foursol(ni, nj, u, f, ff, x, y, axm, ax0, axp)

      use sp_parameters, only: nip, njp, pi

      implicit none

      integer, intent(in) :: ni, nj
      real*8, intent(in), dimension(nip) :: x, axm, axp, ax0
      real*8, intent(in), dimension(njp) :: y
      real*8, intent(in), dimension(nip, njp) :: f
      real*8, intent(out), dimension(nip, njp) :: u, ff

      integer :: i, j, k, j_callax, m, md, mm, ni1, nj1
      real*8 :: hy, slam
      real*8, dimension(njp) :: w, a
      real*8, dimension(nip) :: ap, bp, cp, fp, alf, bet, u1

      data j_callax/0/
      save j_callax
      
!power nj=2*mm definition
      md=nj-1
      do m=1,100
         mm=m
         if(md/2*2.ne.md) then
            write(*,*) 'nj-1 must be power of 2.program is terminated '
            write(*,*) 'nj=', nj
            stop
         endif							 
         md=md/2	 
         if(md.eq.1) EXIT
      enddo

      nj1=nj-1
      ni1=ni-1

      hy=y(2)-y(1)

      do i=2,ni1
         do j=2,nj1
            u(i,j)=f(i,j)
         enddo
      enddo

      do i=2,ni-1
         u(i,2)=u(i,2)-u(i,1)/hy**2
         u(i,nj-1)=u(i,nj-1)-u(i,nj)/hy**2
      enddo

      do j=2,nj-1
         u(2,j)=u(2,j)-u(1,j)*axm(2)
         u(ni-1,j)=u(ni-1,j)-u(ni,j)*axp(ni1)
      enddo

      do i=2,ni1
         do j=2,nj1
            w(j)=u(i,j)
         enddo

         call rft1(mm,w,a,3)

         do k=2,nj1
            ff(i,k)=a(k)
         enddo
      enddo

      do k=2,nj-1

         slam=4.d0*dsin(0.5d0*(k-1.d0)*pi/nj1)**2/hy**2

         do i=2,ni-1
            ap(i-1)=axp(i)
            bp(i-1)=axm(i)
            cp(i-1)=ax0(i)+slam
            fp(i-1)=-ff(i,k)          
         enddo 

         ap(ni-2)=0.d0	   
         bp(1)=0.d0
      	
         call pROG1D(ni-2,u1(2), Ap,Bp,Cp,Fp, ALF,BET)

         do i=2,ni-1
            ff(i,k)=u1(i)	 
         enddo 
      		        
      enddo

      do i=2,ni-1
         do k=2,nj-1
            a(k)=ff(i,k)
         enddo
         call rft1(mm,a,w,6)
         do j=2,nj-1
           u(i,j)=w(j)
         enddo
      enddo

      return
      end subroutine foursol

!----------------------------------------------------------------
      subroutine foursol_ws(ni, nj, u, f, ff, x, y, axm, ax0, axp)

      use sp_parameters, only: nip, njp, pi

      implicit none

      integer, intent(in) :: ni, nj
      real*8, intent(in), dimension(nip) :: x, axm, axp, ax0
      real*8, intent(in), dimension(njp) :: y
      real*8, intent(in), dimension(nip, njp) :: f
      real*8, intent(out), dimension(nip, njp) :: u, ff

      integer :: i, j, k, j_callax, m, md, mm, ni1, nj1
      real*8 :: hy, slam
      real*8, dimension(njp) :: w, a
      real*8, dimension(nip) :: ap, bp, cp, fp, alf, bet, u1

      data j_callax/0/
      save j_callax

! power nj=2*mm definition
      md=nj-1
      do m=1,100
         mm=m
         if(md/2*2.ne.md) then
            write(*,*) 'nj-1 must be power of 2.program is terminated '
            write(*,*) 'nj=', nj
            stop
         endif							 
         md=md/2	 
         if(md.eq.1) EXIT
      enddo
	  
      nj1=nj-1
      ni1=ni-1

      hy=y(2)-y(1)

      do i=2,ni1
         do j=2,nj1
            u(i,j)=f(i,j)
         enddo
      enddo

      do i=2,ni-1
         u(i,2)=u(i,2)-u(i,1)/hy**2
         u(i,nj-1)=u(i,nj-1)-u(i,nj)/hy**2
      enddo

      do j=2,nj-1
         u(2,j)=u(2,j)-u(1,j)*axm(2)
         u(ni-1,j)=u(ni-1,j)-u(ni,j)*axp(ni1)
      enddo

      do i=2,ni1
         do j=2,nj1
            w(j)=u(i,j)
         enddo
         call rft1(mm,w,a,3)
         do k=2,nj1
            ff(i,k)=a(k)
         enddo
      enddo

      do k=2,nj-1
         slam=4.d0*dsin(0.5d0*(k-1.d0)*pi/nj1)**2/hy**2
         do i=2,ni-1
            ap(i-1)=axp(i)
            bp(i-1)=axm(i)
            cp(i-1)=ax0(i)+slam
            fp(i-1)=-ff(i,k)          
         enddo 

         ap(ni-2)=0.d0	   
         bp(1)=0.d0

         call pROG1D(ni-2,u1(2), Ap,Bp,Cp,Fp, ALF,BET)

         do i=2,ni-1
            ff(i,k)=u1(i)
         enddo    		        
      enddo

      do i=2,ni-1
         do k=2,nj-1
            a(k)=ff(i,k)	 
         enddo 

         call rft1(mm,a,w,6)

         do j=2,nj-1
            u(i,j)=w(j)
         enddo 
      enddo 

      return
      end subroutine foursol_ws

!----------------------------------------------------------------
      subroutine foursol_ef(ni,nj,u,f,ff,x,y,axm,ax0,axp)

      use sp_parameters, only: nip, njp, pi

      implicit none

      integer, intent(in) :: ni, nj
      real*8, intent(in), dimension(nip) :: x, axm, axp, ax0
      real*8, intent(in), dimension(njp) :: y
      real*8, intent(in), dimension(nip, njp) :: f
      real*8, intent(out), dimension(nip, njp) :: u, ff

      integer :: i, j, k, j_callax, m, md, mm, ni1, nj1, jk1
      integer, dimension(63, 63, 63) :: ifmatrixx
      real*8 :: hy, slam
      real*8, dimension(njp) :: w, a
      real*8, dimension(nip) :: ap, bp, cp, fp, alf, bet, u1
      real*8, dimension(63, 63, 63) :: fmatrixx

      data j_callax/0/
      save j_callax

! power nj=2*mm definition

      md=nj-1
      do m=1,100
         mm=m
         if(md/2*2.ne.md) then
            write(*,*) 'nj-1 must be power of 2.program is terminated '
            write(*,*) 'nj=', nj
            stop
         endif							 
         md=md/2	 
         if(md.eq.1) EXIT
      enddo
 
      nj1=nj-1
      ni1=ni-1

      hy=y(2)-y(1)

      do i=2,ni1
         do j=2,nj1
            u(i,j)=f(i,j)
         enddo
      enddo

      do i=2,ni-1
         u(i,2)=u(i,2)-u(i,1)/hy**2
         u(i,nj-1)=u(i,nj-1)-u(i,nj)/hy**2
      enddo

      do j=2,nj-1
         u(2,j)=u(2,j)-u(1,j)*axm(2)
         u(ni-1,j)=u(ni-1,j)-u(ni,j)*axp(ni1)
      enddo

      do i=2,ni1
         do j=2,nj1
            w(j)=u(i,j)
         enddo
         call rft1(mm,w,a,3)
         do k=2,nj1
            ff(i,k)=a(k)
         enddo
      enddo

      jk1=63
      if (j_callax.eq.0) then
         fmatrixx=0.
      	 ifmatrixx=0.
         do k=2,nj-1
            slam=4.d0*dsin(0.5d0*(k-1.d0)*pi/nj1)**2/hy**2
            do i=1,ni-2
               if (i.lt.ni-2) fmatrixx(i,i+1,k-1)=axp(i+1)
               if (i.gt.1) fmatrixx(i,i-1,k-1)=axm(i+1)
               fmatrixx(i,i,k-1)=-ax0(i+1)-slam
            enddo 

      	    call inverse_spid(fmatrixx(1:jk1,1:jk1,k-1),
     &         ifmatrixx(1:jk1,1:jk1,k-1),jk1)	
         enddo
         j_callax=1
      endif

      do k=2,nj-1
         u1(1:jk1)=matmul(ifmatrixx(1:jk1,1:jk1,k-1),
     &      ff(2:jk1+1,k))	 
            ff(2:jk1+1,k)=u1(1:jk1)	 
      enddo
 		        
      do i=2,ni-1
         do k=2,nj-1
            a(k)=ff(i,k)	 
         enddo
         call rft1(mm,a,w,6)
         do j=2,nj-1
            u(i,j)=w(j)
         enddo 
      enddo 

      return
      end subroutine foursol_ef

!----------------------------------------------------------------
      subroutine RFT1 (MM, X, Y, MODE)
! Fast Fourier Transform in 1d

      implicit none

      integer, intent(in) :: MM, MODE
      real*8, intent(in), dimension(*) :: X
      real*8, intent(out), dimension(*) :: Y

      integer :: i, N, N1, N2, M, MD, NC, NS, ND, NR
      real*8 :: F, RTTWO, W(2048)
      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, ND, NR
      common /FWORK/ W

! Fixed data preparation
      MD = MM
      if(MODE .NE. 1 .AND. MODE .NE. 4) MD = MM + 1
      call SETFT (MD)
      N1 = 2 ** MM + 1
! Copy input vector into working storage
      do I = 1, N1
         W(ND+I) = X(I) * F
      enddo

      select case(MODE)

      case(1) ! General analysis
         W(NR+N1) = 0.0d0
 
      case(2) ! Dcosine analysis or synthesis
         call RCA

      case(3) ! Dsine analysis or synthesis
         call RSA
         W(NR+1) = X(1)
         W(NR+N1) = X(N1)

      case(4) ! General synthesis
         call RPS
         W(NR+N1) = W(NR+1)

      end select

! Copy into output vector
      do I = 1, N1
         Y(I) = W(NR+I)
      enddo

      return
      end subroutine RFT1

!-----------------------------------------------------------------------
      subroutine RPA
! Real analysis

      implicit none

      integer :: L, IT, NQ, NT, NIR, NOR, NI, NII, NOI, NO1, NO2, 
     &   NI1, NI2
      real*8 :: CC, SS, RE, AI
      integer :: N, N2, M, MD, NC, NS, NDATA, NRES
      real*8 :: F, RTTWO, W(2048)

      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, NDATA, NRES
      common /FWORK/ W

      NQ = N2
      NIR = NDATA
      NOR = NDATA + N
      do L = 1, M
         NII = NIR + N2
         NOI = NOR + N2
! TPRIME = 0.
         do IT = 1, NQ
            W(NOR+IT) = W(NIR+IT) + W(NIR+IT+NQ)
            W(NOI+IT) = W(NIR+IT) - W(NIR+IT+NQ)
         enddo
! 0 .LT. TPRIME .LE. NP / 2
         NO1 = NQ
         NO2 = N2 - NQ
         if(NO1 - NO2) 20, 40, 60
   20    continue
         NI1 = NO1 + NO1
         NI2 = NI1 + NQ
         CC = W(NC+NO1)
         SS = W(NS+NO1)
         do IT = 1, NQ
            RE = CC * W(NIR+NI2+IT) - SS * W(NII+NI2+IT)
            AI = SS * W(NIR+NI2+IT) + CC * W(NII+NI2+IT)
            W(NOR+NO1+IT) = W(NIR+NI1+IT) + RE
            W(NOR+NO2+IT) = W(NIR+NI1+IT) - RE
            W(NOI+NO1+IT) =   W(NII+NI1+IT) + AI
            W(NOI+NO2+IT) = - W(NII+NI1+IT) + AI
         enddo
         NO1 = NO1 + NQ
         NO2 = NO2 - NQ
         if(NO1 - NO2) 20, 40, 60
! TPRIME = NP/2.
   40    continue
         do IT = 1, NQ
            W(NOR+NO1+IT) = W(NII   +IT)
            W(NOI+NO1+IT) = W(NII+NQ+IT)
         enddo
   60    continue
         NT = NIR
         NIR = NOR
         NOR = NT
         NQ = NQ / 2
      enddo

      NRES = NIR

      return
      end subroutine RPA

!--------------------------------------------------------------------
      subroutine RPS
! Real synthesis

      implicit none

      integer :: NQ, NT, NIR, L, NII, NOI, NOR, IT, NI1, NI2, NO1, NO2
      real*8 :: CC, SS, RE, AI
      integer :: N, N2, M, MD, NC, NS, NDATA, NRES
      real*8 :: F, RTTWO, W(2048)

      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, NDATA, NRES
      common /FWORK/ W

      NQ = 1
      NIR = NDATA
      NOR = NDATA + N
      L = M
      do
         NII = NIR + N2
         NOI = NOR + N2
! TPRIME = 0.
         do IT = 1, NQ
            W(NOR+IT) = W(NIR+IT) + W(NII+IT)
            W(NOR+IT+NQ) = W(NIR+IT) - W(NII+IT)
         enddo
! 0 .LT. TPRIME .LT. NP/2.
         NI1 = NQ
         NI2 = N2 - NQ
         if(NI1 - NI2) 40, 80, 120
   40    continue
         NO1 = NI1 + NI1
         NO2 = NO1 + NQ
         CC = W(NC+NI1)
         SS = W(NS+NI1)
         do IT = 1, NQ
            W(NOR+NO1+IT) = W(NIR+NI1+IT) + W(NIR+NI2+IT)
            RE            = W(NIR+NI1+IT) - W(NIR+NI2+IT)
            W(NOI+NO1+IT) = W(NII+NI1+IT) - W(NII+NI2+IT)
            AI            = W(NII+NI1+IT) + W(NII+NI2+IT)
            W(NOR+NO2+IT) = CC * RE + SS * AI
            W(NOI+NO2+IT) = CC * AI - SS * RE
         enddo
         NI1 = NI1 + NQ
         NI2 = NI2 - NQ
         if(NI1 - NI2) 40, 80, 120
   80    continue
         do IT = 1, NQ
            W(NOI   +IT) = W(NIR+NI1+IT) * 2.0d0
            W(NOI+NQ+IT) = W(NII+NI1+IT) * 2.0d0
         enddo
  120    continue
         NT = NIR
         NIR = NOR
         NOR = NT
         NQ = NQ + NQ
         L = L - 1
         if(L .LE. 0) EXIT
      enddo
      NRES = NIR

      return
      end subroutine RPS

!----------------------------------------------------------------
      subroutine RSA
! Real odd analysis or synthesis

      implicit none

      integer :: L, IT, NQ, NI, NR, NT, NIB, NOB, NIM, NIE, NQH, NQL, 
     &   NOM, NO1, NO2
      real*8 :: CC, SS, RE, AI
      integer :: N, N2, M, MD, NC, NS, NDATA, NRES
      real*8 :: F, RTTWO, W(2048)

      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, NDATA, NRES
      common /FWORK/ W

      NQ = N2
      NIB = NDATA + 1
      NOB = NDATA + N2 + 2
      do L = 1, M
         NIM = NIB + NQ
         NIE = NIM + NQ
         NQH = NQ / 2
         NQL = (NQ - 1) / 2
         NOM = NOB + NQL + 1
! TPRIME = 0.
         if(NQL .NE. 0) then
            do IT = 1, NQL
               W(NOB+IT) = W(NIB+IT) - W(NIM-IT)
               W(NOM+IT) = W(NIB+IT) + W(NIM-IT)
            enddo
         endif
         if(NQH .NE. NQL) W(NOM) = 2.0d0 * W(NIB+NQH)
! 0 .LT. TPRIME .LT. NP/2.
         NO1 = NQ
         NO2 = N2 - NQ
         if(NO1 - NO2) 40, 80, 120
   40    continue
         NI = NO1 + NO1
         W(NOB+NO1) =   W(NIB+NI) + W(NIM+NI)
         W(NOB+NO2) = - W(NIB+NI) + W(NIM+NI)
         if(NQL .NE. 0)then
            CC = W(NC+NO1)
            SS = W(NS+NO1)
            do IT=1,NQL
               RE = CC * W(NIE+NI-IT) - SS * W(NIM+NI-IT)
               AI = SS * W(NIE+NI-IT) + CC * W(NIM+NI-IT)
               W(NOM+NO1+IT) = W(NIM+NI+IT) - RE
               W(NOM+NO2+IT) = W(NIM+NI+IT) + RE
               W(NOB+NO1+IT) =   W(NIB+NI+IT) + AI
               W(NOB+NO2+IT) = -  W(NIB+NI+IT) + AI
            enddo
         endif
         if(NQH .NE. NQL) then
            NR = NO1 / 2
            W(NOM+NO1)=2.d0*(W(NS+NR)*W(NIM+NI+NQH) + 
     &                       W(NC+NR)*W(NIB+NI+NQH))
            W(NOM+NO2)=2.d0*(W(NC+NR)*W(NIM+NI+NQH) - 
     &                       W(NS+NR)*W(NIB+NI+NQH))
         endif
         NO1 = NO1 + NQ
         NO2 = NO2 - NQ
         if(NO1 - NO2) 40, 80, 120
   80    continue
         W(NOB+NO1) = W(NIM)
         if(NQL .NE. 0) then
            do IT = 1, NQL
               W(NOM+NO1+IT) = W(NIM+IT)
               W(NOB+NO1+IT) = W(NIE-IT)
            enddo
         endif
         if(NQH .NE. NQL) W(NOM+NO1) = RTTWO * W(NIM+NQH)
  120    continue
         NT = NIB
         NIB = NOB
         NOB = NT
         NQ = NQH
      enddo
      NRES = NIB - 1

      return
      end subroutine RSA

!----------------------------------------------------------------
      subroutine RCA
! Real odd analysis or synthesis

      implicit none

      integer :: L, IT, NQ, NI, NR, NT, NIB, NOB, NIM, NIE, NQH, NQL, 
     &   NOM, NO1, NO2
      real*8 :: CC, SS, RE, AI
      integer :: N, N2, M, MD, NC, NS, NDATA, NRES
      real*8 :: F, RTTWO, W(2048)

      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, NDATA, NRES
      common /FWORK/ W

      NQ = N2
      NIB = NDATA + 1
      NOB = NDATA + N2 + 2
      do L = 1, M
         NIM = NIB + NQ
         NIE = NIM + NQ
         NQH = NQ / 2
         NQL = (NQ - 1) / 2
         NOM = NOB + NQL + 1
! TPRIME = 0.
         W(NOB   ) = W(NIB) + W(NIM)
         W(NOB+N2) = W(NIB) - W(NIM)
         if(NQL .NE. 0) then
            do IT = 1, NQL
               W(NOB+IT) = W(NIB+IT) + W(NIM-IT)
               W(NOM+IT) = W(NIB+IT) - W(NIM-IT)
            enddo
         endif
         if(NQH .NE. NQL) W(NOM) = W(NIB+NQH) * 2.d0
! 0 .LT. TPRIME .LT. NP / 2
         NO1 = NQ
         NO2 = N2 - NQ
         if(NO1 - NO2) 40, 80, 120
   40    continue
         NI = NO1 + NO1
         W(NOB+NO1) = W(NIB+NI) + W(NIM+NI)
         W(NOB+NO2) = W(NIB+NI) - W(NIM+NI)
         if(NQL .NE. 0) then
            CC = W(NC+NO1)
            SS = W(NS+NO1)
            do IT = 1, NQL
               RE = CC * W(NIM+NI-IT) - SS * W(NIE+NI-IT)
               AI = SS * W(NIM+NI-IT) + CC * W(NIE+NI-IT)
               W(NOB+NO1+IT) = W(NIB+NI+IT) + RE
               W(NOB+NO2+IT) = W(NIB+NI+IT) - RE
               W(NOM+NO1+IT) =   W(NIM+NI+IT) - AI
               W(NOM+NO2+IT) = - W(NIM+NI+IT) - AI
            enddo
         endif
         if(NQH .NE. NQL) then
            NR = NO1 / 2
            W(NOM+NO1)=2.d0*(W(NC+NR)*W(NIB+NI+NQH) - 
     &                       W(NS+NR)*W(NIM+NI+NQH))
            W(NOM+NO2)=2.d0*(W(NS+NR)*W(NIB+NI+NQH) + 
     &                       W(NC+NR)*W(NIM+NI+NQH))
         endif
         NO1 = NO1 + NQ
         NO2 = NO2 - NQ
         if(NO1 - NO2) 40, 80, 120
! TPRIME = NP/2
   80    continue
         W(NOB+NO1) = W(NIB+N2)
         if(NQL .NE. 0) then
            do IT = 1, NQL
               W(NOB+NO1+IT) =   W(NIM+IT)
               W(NOM+NO1+IT) = - W(NIE-IT)
            enddo
         endif
         if(NQH .NE. NQL) W(NOM+NO1) = RTTWO * W(NIM+NQH)
  120    continue
         NT = NIB
         NIB = NOB
         NOB = NT
         NQ = NQH
      enddo

      NRES = NIB - 1

      return
      end subroutine RCA

!----------------------------------------------------------------
      subroutine SETFT (MM)
! Prepare fixed data for Fast Fourier Transform

      use sp_parameters, only: twopi

      implicit none

      integer, intent(in) :: MM
      integer :: I, N1
      real*8 :: DA, A
      integer :: N, N2, M, MD, NC, NS, NDATA, NRES
      real*8 :: F, RTTWO, W(2048)

      common /FDATA/ F, RTTWO, N, N2, M, NC, NS, NDATA, NRES
      common /FWORK/ W

      if(MM .LT. 0) return

      M = MM
      N = 2 ** M
      N2 = N / 2
      F = 1.d0 / dFLOAT (N)
      DA = twopi * F
      RTTWO = 2.0d0
      RTTWO = DSQRT (RTTWO)
      F = DSQRT (F)
      N1 = (N2 - 1) / 2
      NC = 0
      NS = N1
      NDATA = N1 + N1
      if(N1 .LE. 0) return
      do I = 1, N1
         A = dFLOAT(I) * DA
         W(I) = DCOS(A)
         W(N1+I) = DSIN(A)
      enddo

      return
      end subroutine SETFT





