!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine remesh(erro,errpsi,imov)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'
	common
     *  /c_kpr/kpr

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

             iswtch=0
	     
	     if(iter.ne.1) then

               if(imax.eq.1 .and. erro.lt.5.d-4) iswtch=1

               if(kpr.eq.1) then
              write(*,*) 'mesh:imax jmax iswth',imax,jmax,iswtch
               endif	     
	     else

               iswtch=0
	       
	     endif

             if(ngav.le.0) iswtch=0

             if(imax.ge.2) then
                call prgrid(imax,jmax,erro,imov,errpsi)
             else
                if(iswtch.eq.0) then
	             call regrid0(erro,imov,errpsi)
                else
	             call regrid(erro,imov,errpsi)
                endif
             endif

		 return
		 end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine prgrid(imax,jmax,erro,imov,errpsi)

           include 'double.inc'
           include 'dim.inc'
           parameter(nshp=ntp+1)
           include 'compol.inc'

         dimension ron(ntp),rop(ntp),rn(ntp),zn(ntp),tetn(ntp)
         dimension rob(ntp),teti(ntp),roi(ntp)
         dimension xs(nshp),ys(nshp),fun(nshp),dp(5)
	common
     *  /c_kpr/kpr

!          atan(xx)=datan(xx)
          acos(xx)=dacos(xx)
          cos(xx)=dcos(xx)
          sin(xx)=dsin(xx)
          sqrt(xx)=dsqrt(xx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! new position of magnetic axis !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(kpr.eq.1) then
           write(*,*) 'prgrid:::'
            endif
           rm0=r(1,1)
           zm0=z(1,1)
           psim0=psi(1,1)

           nsh=1

           xs(nsh)=r(imax,jmax)
           ys(nsh)=z(imax,jmax)
           fun(nsh)=psi(imax,jmax)

          if(imax.eq.1) then
           do j=2,nt1
           nsh=nsh+1
           xs(nsh)=r(2,j)
           ys(nsh)=z(2,j)
           fun(nsh)=psi(2,j)
           enddo
          elseif(imax.eq.2) then

           nsh=nsh+1
           xs(nsh)=r(1,2)
           ys(nsh)=z(1,2)
           fun(nsh)=psi(1,2)
           do j=2,nt1
           nsh=nsh+1
           xs(nsh)=r(3,j)
           ys(nsh)=z(3,j)
           fun(nsh)=psi(3,j)
           enddo

          else

        do k=-1,1
          i= imax+k
        do l=-1,1
          j= jmax+l
          if(i.ne.imax .OR. j.ne.jmax) then
           nsh=nsh+1
           xs(nsh)=r(i,j)
           ys(nsh)=z(i,j)
           fun(nsh)=psi(i,j)
	    endif
        enddo
        enddo
           
          endif

          call deriv5(xs,ys,fun,nsh,5,dp)

        DET = dp(3)*dp(5) - dp(4)**2
        Rm = Xs(1) + ( dp(2)*dp(4) - dp(1)*dp(5) )/DET
        Zm = Ys(1) + ( dp(1)*dp(4) - dp(2)*dp(3) )/DET

          erro=dsqrt( (rm-rm0)**2+(zm-zm0)**2 )

           rmx=r(imax,jmax)
           zmx=z(imax,jmax)

       psim=fun(1)+ dp(1)*(rm-rmx) + dp(2)*(zm-zmx)
     +          + 0.5d0*dp(3)*(rm-rmx)*(rm-rmx)
     +          +       dp(4)*(rm-rmx)*(zm-zmx)
     +          + 0.5d0*dp(5)*(zm-zmx)*(zm-zmx)

            if(kpr.eq.1) then
         write(6,*) 'rm zm psim',rm,zm,psim
         write(6,*) 'prgrid:errma',erro
            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! nomalized poloidal flux definition !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          errpsi=0.d0

         do i=1,iplas
         do j=1,nt

	    psnn=psin(i,j)
          psin(i,j)=(psi(i,j)-psip)/(psim-psip)
		delpsn=dabs(psin(i,j)-psnn)

           if(delpsn.gt.errpsi) then
             ierm=i
             jerm=j
	       errpsi=delpsn
           endif

		!errpsi=dmax1(errpsi,delpsn)

         enddo
         enddo
        if(imov.eq.0) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! definition of new angle grid teta(j) !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,nt
           tetn(j)=teta(j)
         enddo

        do j=1,nt
         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm
         robj=sqrt(drx**2+dzx**2)
         rob(j)=robj
         tetp=acos(drx/robj)
          if(dzx.lt.0.d0) then
          teta(j)=-tetp
          else
          teta(j)=tetp
          endif
        enddo

            do j=2,nt  !<---

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

	      enddo      !<---

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi

ccc    grid moving along rais

          psinm0=psin(1,2)

	    do i=3,iplas1
	     if(psia(i).lt.psinm0) then
	      !ibeg=i
	      ibeg=i+1
	      go to 10
	     endif
	    enddo

 10       continue

        do i=ibeg,iplas-1
          zps =psia(i)
         do j=1,nt     !<------ j-loop
          do is=1,iplas1
            ps0 =psin(is,j)
            pspl=psin(is+1,j)
            ro0 =ro(is,j)
            ropl=ro(is+1,j)
 	     if(zps.le.ps0 .and. zps.ge.pspl) then
		  isc=is
		  go to 822
	     endif
          enddo

 822    continue

            grad=-(pspl-ps0)/(ropl-ro0)

	     zro=ro0-(psia(i)-ps0)/grad
           ron(j)=zro
         enddo   !<------ j-loop

        do j=1,nt

            rnj=rm0+ron(j)*cos(tetn(j))
            znj=zm0+ron(j)*sin(tetn(j))

         drx=rnj-rm
         dzx=znj-zm

         roi(j)=sqrt(drx**2+dzx**2)
         tetp=acos(drx/roi(j))
          if(dzx.lt.0.d0) then
          teti(j)=-tetp
          else
          teti(j)=tetp
          endif

        enddo

        do j=2,nt  !<---

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

	  enddo      !<---

           teti(1)=teti(nt1)-2.d0*pi
           teti(nt)=teti(2)+2.d0*pi

           roi(1)=roi(nt1)
           roi(nt)=roi(2)

        do j=2,nt1 !<------ j-loop

              tetv=teta(j)

            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi

          do jj=1,nt1

           tt0 =teti(jj)
           ttpl=teti(jj+1)

           ro0 =roi(jj)
           ropl=roi(jj+1)

	       if(tetv.le.ttpl .and. tetv.ge.tt0) then
              zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
              go to 1047
	       endif
          enddo

 1047      continue

		 rop(j)=zro

           r(i,j)=rop(j)*cos(teta(j))+rm
           z(i,j)=rop(j)*sin(teta(j))+zm

         enddo  !<------  j-loop

        enddo  ! i-loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           i=2
           zps =psia(i)
           ps0=psip+zps*(psim-psip)

             do j=1,nt

         ttj=teta(j)
        ro2j=(ps0-psim)/( 0.5d0*dp(3)*cos(ttj)**2
     +                +     dp(4)*cos(ttj)*sin(ttj)
     +                +   0.5d0*dp(5)*sin(ttj)**2 )

        rodeb=ro2j

          rop(j)=dsqrt(ro2j)

             r(i,j)=rm+rop(j)*cos(teta(j))
             z(i,j)=zm+rop(j)*sin(teta(j))

             enddo  ! j-loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(ibeg.gt.3) then

         do j=2,nt1

          ros1=(r(2,j)-rm)**2+(z(2,j)-zm)**2
          ros2=(r(ibeg,j)-rm)**2+(z(ibeg,j)-zm)**2
          psa1=psia(2)
          psa2=psia(ibeg)

          do i=3,ibeg-1

           psai=psia(i)
           rosi=( (psa2-psai)*ros1+(psai-psa1)*ros2 )/(psa2-psa1)
           rosi=dsqrt(rosi)

           r(i,j)=rm+rosi*cos(teta(j))
           z(i,j)=zm+rosi*sin(teta(j))

          enddo

         enddo

        endif

         do i=2,iplas1
          do j=2,nt1

            ro(i,j)=dsqrt((r(i,j)-rm)**2+(z(i,j)-zm)**2)

         enddo
          enddo

            do j=2,nt1

            ro(iplas,j)=rob(j)
            r(iplas,j)=rob(j)*cos(teta(j))+rm
            z(iplas,j)=rob(j)*sin(teta(j))+zm

           ro(1,j)=0.d0
           r(1,j)=rm
           z(1,j)=zm

            enddo

          do 25 i=1,iplas

           ro(i,1)=ro(i,nt1)
           ro(i,nt)=ro(i,2)

           r(i,1)=r(i,nt1)
           r(i,nt)=r(i,2)

           z(i,1)=z(i,nt1)
           z(i,nt)=z(i,2)

 25       continue

           do 30 i=1,nr
           do 30 j=1,nt

            psin(i,j)=psia(i)
            psi(i,j)=psip+psin(i,j)*(psim-psip)

 30       continue

        call wrb

       return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine regrid(erro,imov,errpsi)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

         dimension ron(ntp),rop(ntp),rn(ntp),zn(ntp),tetn(ntp)
         dimension rob(ntp),teti(ntp),roi(ntp)
         dimension roplt(nrp,ntp)
         dimension dp(5)
	common
     *  /c_kpr/kpr

!         xzer(x1,x2,v1,v2) = (x1*v2-x2*v1)/(v2-v1)

         cos(xx)=dcos(xx)
         sin(xx)=dsin(xx)
!         atan(xx)=datan(xx)
!         asin(xx)=dasin(xx)
         sqrt(xx)=dsqrt(xx)
         acos(xx)=dacos(xx)

          if(kpr.eq.1) then
           write(*,*) 'regrid:::'
          endif
           if(ngav.eq.0) then
            alfa=1.0d0
           else
            alfa=1.0d0
           endif

          rmold=rm
          zmold=zm

          if(kpr.eq.1) then
           write(*,*) 'rm,zm,psim',rm,zm,psim
          endif

             call axdef(rma,zma,psima,dp)

          if(kpr.eq.1) then
           write(*,*) 'rma,zma,psima',rma,zma,psima
          endif
	if (isnan(rma)) then
			write(*,*) 'mag axis major radius is NaN'
			call err_catch_a
	endif
	     rm=rma
	     zm=zma
	     psim=psima

             erro=dsqrt((rmold-rm)**2+(zmold-zm)**2)

          if(kpr.eq.1) then
           write(6,*) 'erro mag.axis',erro
          endif

           errpsi=0.d0

         do i=1,iplas
         do j=1,nt

	    psnn=psin(i,j)
          psin(i,j)=(psi(i,j)-psip)/(psim-psip)
		delpsn=dabs(psin(i,j)-psnn)

           if(delpsn.gt.errpsi) then

             ierm=i
             jerm=j
	       errpsi=delpsn

           endif

		!errpsi=dmax1(errpsi,delpsn)

         enddo
         enddo

          if(kpr.eq.1) then
         write(*,*) 'i,j errpsi',ierm,jerm,errpsi,(psim-psip)
          endif
         if(imov.eq.0) return

cc      definition of new angle grid teta(j)

         do j=1,nt
           tetn(j)=teta(j)
         enddo

        do j=1,nt

         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm

         robj=sqrt(drx**2+dzx**2)
         rob(j)=robj
         tetp=acos(drx/robj)
          if(dzx.lt.0.d0) then
          teta(j)=-tetp
          else
          teta(j)=tetp
          endif

        enddo

            do j=2,nt  !<---

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

	      enddo      !<---

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           i=2

           zps =psia(i)
           ps0=psip+zps*(psim-psip)

             do j=1,nt

         ttj=teta(j)

        ro2j=(ps0-psim)/( 0.5d0*dp(3)*cos(ttj)**2
     +                +     dp(4)*cos(ttj)*sin(ttj)
     +                + 0.5d0*dp(5)*sin(ttj)**2 )

          rop(j)=dsqrt(ro2j)

             r(i,j)=rm+rop(j)*cos(teta(j))
             z(i,j)=zm+rop(j)*sin(teta(j))

             enddo  ! j-loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ccc    grid moving along rais

          do i=2,iplas-1

           zps =psia(i)

          do j=1,nt     !<------ j-loop

!              do is=1,iplas1
!
!           ps0 =psin(is,j)
!           !psmn=psin(is-1,j)
!           pspl=psin(is+1,j)
!
!           ro0 =ro(is,j)
!           !romn=ro(is-1,j)
!           ropl=ro(is+1,j)
!
!	              if(zps.le.ps0 .and. zps.ge.pspl) then
!              zro=((pspl-zps)*ro0-(ps0-zps)*ropl)/(pspl-ps0)
!                go to 822
!	              endif
!            !zro=qvadin(zps, psmn,ps0,pspl, romn,ro0,ropl)
!                 enddo
! 822              continue

            ps0 =psin(i,j)
            psmn=psin(i-1,j)
            pspl=psin(i+1,j)

           ro0 =ro(i,j)
           romn=ro(i-1,j)
           ropl=ro(i+1,j)

             gradpl=-(pspl-ps0)/(ropl-ro0)
             gradmn=-(ps0-psmn)/(ro0-romn)

             grad=dmax1(gradpl,gradmn)

           if( ngav.eq.0 )  then
               alfa=1.5d0
           else
               alfa=1.0d0+0.45d0*dabs(q(i)-q(1))/q(1)
           end if

           if(alfa.lt.1.5d0) alfa=1.5d0 

	     if(ngav.gt.0) grad=grad*alfa

	      zro=ro0-(psia(i)-ps0)/grad

           ron(j)=zro

          enddo  !<------ j-loop

        do j=1,nt

            rnj=rmold+ron(j)*cos(tetn(j))
            znj=zmold+ron(j)*sin(tetn(j))

         drx=rnj-rm
         dzx=znj-zm

         roi(j)=sqrt(drx**2+dzx**2)
         tetp=acos(drx/roi(j))
          if(dzx.lt.0.d0) then
          teti(j)=-tetp
          else
          teti(j)=tetp
          endif

        enddo

            do j=2,nt  !<---

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

	      enddo      !<---

           teti(1)=teti(nt1)-2.d0*pi
           teti(nt)=teti(2)+2.d0*pi

           roi(1)=roi(nt1)
           roi(nt)=roi(2)

          do j=2,nt1 !<------ j-loop

              tetv=teta(j)

            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi

                do jj=1,nt1

           tt0 =teti(jj)
           !ttmn=teti(j-1)
           ttpl=teti(jj+1)

           ro0 =roi(jj)
           !romn=roi(j-1)
           ropl=roi(jj+1)

	            if(tetv.le.ttpl .and. tetv.ge.tt0) then
              zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
               go to 1047
	            endif
                enddo

 1047      continue

           ! zro=qvadin(teta(j), ttmn,tt0,ttpl, romn,ro0,ropl)

		 rop(j)=zro

                    roerr=dabs(ron(j)-ro(i,j))

	              if(roerr.gt.erro) then
	                 erro=roerr
                      ierm=i
                      jerm=j
	              endif

             roplt(i,j)=ron(j)-ro(i,j)

           r(i,j)=rop(j)*cos(teta(j))+rm
           z(i,j)=rop(j)*sin(teta(j))+zm

          enddo  !<------  j-loop

          enddo  ! i-loop

          if(kpr.eq.1) then
         write(6,*) 'erro max.',erro
         write(6,*) 'i,j erro',ierm,jerm,erro
          endif
             do i=2,iplas1
             do j=2,nt1

            ro(i,j)=dsqrt((r(i,j)-rm)**2+(z(i,j)-zm)**2)
            ronor(i,j)=ro(i,j)/rob(j)

             enddo  !<------  j-loop
             enddo  ! i-loop

            do j=2,nt1

           ro(iplas,j)=rob(j)

            r(iplas,j)=rob(j)*cos(teta(j))+rm
            z(iplas,j)=rob(j)*sin(teta(j))+zm


           ro(1,j)=0.d0

           r(1,j)=rm
           z(1,j)=zm

            enddo

          do 25 i=1,iplas

           ro(i,1)=ro(i,nt1)
           ro(i,nt)=ro(i,2)

           ronor(i,1)=ronor(i,nt1)
           ronor(i,nt)=ronor(i,2)

           r(i,1)=r(i,nt1)
           r(i,nt)=r(i,2)

           z(i,1)=z(i,nt1)
           z(i,nt)=z(i,2)

           roplt(i,1)=roplt(i,nt1)
           roplt(i,nt)=roplt(i,2)

 25       continue

          errm=erro

           do 30 i=1,nr
           do 30 j=1,nt
              psin(i,j)=psia(i)
              psi(i,j)=psip+psin(i,j)*(psim-psip)
  30       continue

        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine regrid0(erro,imov,errpsi)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

         dimension ron(ntp),rop(ntp),rn(ntp),zn(ntp),tetn(ntp)
         dimension rob(ntp),teti(ntp),roi(ntp)
         dimension roplt(nrp,ntp)
         dimension dp(5)
	common
     *  /c_kpr/kpr

!         xzer(x1,x2,v1,v2) = (x1*v2-x2*v1)/(v2-v1)

         cos(xx)=dcos(xx)
         sin(xx)=dsin(xx)
!         atan(xx)=datan(xx)
!         asin(xx)=dasin(xx)
         acos(xx)=dacos(xx)
         sqrt(xx)=dsqrt(xx)

          if(kpr.eq.1) then
           write(*,*) 'regrid0:::'
          endif
           if(ngav.eq.0) then
            alfa=1.0d0
           else
            alfa=1.0d0
           endif

          rmold=rm
          zmold=zm

          if(kpr.eq.1) then
           write(*,*) 'rm,zm,psim',rm,zm,psim
          endif
             call axdef(rma,zma,psima,dp)

          if(kpr.eq.1) then
           write(*,*) 'rma,zma,psima',rma,zma,psima

          endif
	if (isnan(rma))	then
			write(*,*) 'mag axis major radius is NaN'
			call err_catch_a
	endif
	     rm=rma
	     zm=zma
	     psim=psima

             erro=dsqrt((rmold-rm)**2+(zmold-zm)**2)

          if(kpr.eq.1) then
          write(*,*) 'erro mag.axis',erro
          endif

           errpsi=0.d0

         do i=1,iplas
         do j=1,nt

	    psnn=psin(i,j)
          psin(i,j)=(psi(i,j)-psip)/(psim-psip)
		delpsn=dabs(psin(i,j)-psnn)

           if(delpsn.gt.errpsi) then

             ierm=i
             jerm=j
	       errpsi=delpsn

           endif

		!errpsi=dmax1(errpsi,delpsn)

         enddo
         enddo

          if(kpr.eq.1) then
         write(*,*) 'i,j errpsi',ierm,jerm,errpsi,(psim-psip)
          endif

         if(imov.eq.0) return

cc      definition of new angle grid teta(j)

         do j=1,nt

           tetn(j)=teta(j)

         enddo

        do j=1,nt

         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm

         robj=sqrt(drx**2+dzx**2)
         rob(j)=robj
         tetp=acos(drx/robj)
          if(dzx.lt.0.d0) then
          teta(j)=-tetp
          else
          teta(j)=tetp
          endif

        enddo

            do j=2,nt  !<---

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

	      enddo      !<---

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           i=2

           zps =psia(i)
           ps0=psip+zps*(psim-psip)

             do j=1,nt

         ttj=teta(j)

        ro2j=(ps0-psim)/( 0.5d0*dp(3)*cos(ttj)**2
     +                +     dp(4)*cos(ttj)*sin(ttj)
     +                + 0.5d0*dp(5)*sin(ttj)**2 )

          rop(j)=dsqrt(ro2j)

             r(i,j)=rm+rop(j)*cos(teta(j))
             z(i,j)=zm+rop(j)*sin(teta(j))

             enddo  ! j-loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ccc    grid moving along rais

          do i=3,iplas-1
          !do i=2,iplas-1

           zps =psia(i)

          do j=1,nt     !<------ j-loop

              do is=1,iplas1

           ps0 =psin(is,j)
!           !psmn=psin(is-1,j)
           pspl=psin(is+1,j)

           ro0 =ro(is,j)
!           !romn=ro(is-1,j)
           ropl=ro(is+1,j)

	              if(zps.le.ps0 .and. zps.ge.pspl) then
              zro=((pspl-zps)*ro0-(ps0-zps)*ropl)/(pspl-ps0)
                go to 822
	              endif
!            !zro=qvadin(zps, psmn,ps0,pspl, romn,ro0,ropl)
                 enddo
 822              continue


!            ps0 =psin(i,j)
!            psmn=psin(i-1,j)
!            pspl=psin(i+1,j)

!           ro0 =ro(i,j)
!           romn=ro(i-1,j)
!           ropl=ro(i+1,j)

!             gradpl=-(pspl-ps0)/(ropl-ro0)
!             gradmn=-(ps0-psmn)/(ro0-romn)

!             grad=dmax1(gradpl,gradmn)

!           if( ngav.eq.0 )  then
!               alfa=1.5d0
!           else
!               alfa=1.0d0+0.45d0*dabs(q(i)-q(1))/q(1)
!           end if

!                 if(alfa.lt.1.5d0) alfa=1.5d0 

!	     if(ngav.gt.0) grad=grad*alfa

!	      zro=ro0-(psia(i)-ps0)/grad

           ron(j)=zro

          enddo  !<------ j-loop

        do j=1,nt

            rnj=rmold+ron(j)*cos(tetn(j))
            znj=zmold+ron(j)*sin(tetn(j))

         drx=rnj-rm
         dzx=znj-zm

         roi(j)=sqrt(drx**2+dzx**2)
         tetp=acos(drx/roi(j))
          if(dzx.lt.0.d0) then
          teti(j)=-tetp
          else
          teti(j)=tetp
          endif

        enddo

            do j=2,nt  !<---

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

         if(teti(j).lt.teti(j-1)) then
            teti(j)=teti(j)+2.d0*pi
         endif

	      enddo      !<---


           teti(1)=teti(nt1)-2.d0*pi
           teti(nt)=teti(2)+2.d0*pi

           roi(1)=roi(nt1)
           roi(nt)=roi(2)

          do j=2,nt1 !<------ j-loop

              tetv=teta(j)

            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.lt.teti(1)) tetv=tetv+2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi
            if(tetv.gt.teti(nt)) tetv=tetv-2.d0*pi

                do jj=1,nt1

           tt0 =teti(jj)
           !ttmn=teti(j-1)
           ttpl=teti(jj+1)

           ro0 =roi(jj)
           !romn=roi(j-1)
           ropl=roi(jj+1)

	            if(tetv.le.ttpl .and. tetv.ge.tt0) then
              zro=((ttpl-tetv)*ro0+(tetv-tt0)*ropl)/(ttpl-tt0)
               go to 1047
	            endif
                enddo

 1047      continue

           ! zro=qvadin(teta(j), ttmn,tt0,ttpl, romn,ro0,ropl)

		 rop(j)=zro

                    roerr=dabs(ron(j)-ro(i,j))

	              if(roerr.gt.erro) then
	                 erro=roerr
                      ierm=i
                      jerm=j
	              endif

             roplt(i,j)=ron(j)-ro(i,j)

           r(i,j)=rop(j)*cos(teta(j))+rm
           z(i,j)=rop(j)*sin(teta(j))+zm

          enddo  !<------  j-loop

          enddo  ! i-loop

          if(kpr.eq.1) then
         write(*,*) 'erro max.',erro
         write(*,*) 'i,j erro',ierm,jerm,erro
          endif
             do i=2,iplas1
             do j=2,nt1

            ro(i,j)=dsqrt((r(i,j)-rm)**2+(z(i,j)-zm)**2)
            ronor(i,j)=ro(i,j)/rob(j)

             enddo  !<------  j-loop
             enddo  ! i-loop

            do j=2,nt1

           ro(iplas,j)=rob(j)

            r(iplas,j)=rob(j)*cos(teta(j))+rm
            z(iplas,j)=rob(j)*sin(teta(j))+zm


           ro(1,j)=0.d0

           r(1,j)=rm
           z(1,j)=zm

            enddo

          do 25 i=1,iplas

           ro(i,1)=ro(i,nt1)
           ro(i,nt)=ro(i,2)

           ronor(i,1)=ronor(i,nt1)
           ronor(i,nt)=ronor(i,2)

           r(i,1)=r(i,nt1)
           r(i,nt)=r(i,2)

           z(i,1)=z(i,nt1)
           z(i,nt)=z(i,2)

           roplt(i,1)=roplt(i,nt1)
           roplt(i,nt)=roplt(i,2)

 25       continue

              errm=erro

           do 30 i=1,nr
           do 30 j=1,nt

            psin(i,j)=psia(i)
            psi(i,j)=psip+psin(i,j)*(psim-psip)

  30       continue

           !call wrb


        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_1(igdf,nstep)
         
         use bnd_modul       

           include 'double.inc'
           include 'dim.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           include 'compol.inc'
!           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
	common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh

           real*8 robn(nbtabp),tetbn(nbtabp)
           real*8 RRK(nb4),CCK(nb4),WRK(nb6),CWK(4)

!            cos(xx)=dcos(xx)
!            acos(xx)=dacos(xx)
!            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)

!      if(.not.allocated(rbtab))
!     %allocate( rbtab(nbtab), zbtab(nbtab) )
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nbsh.eq.1) then
        
!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname)
!         !open(1,file='tab_bnd.dat')
!	       read(1,*) nbtab
!	      do ib=1,nbtab
!	         read(1,*) rbtab(ib),zbtab(ib)
!	      enddo
!         close(1)

           !call arc_x_bnd(nt,nbtab,rbtab,zbtab)  !!!!!!!!
           call arc_x_bnd(nt)  !!!!!!!!

!!!!!! new position for magn.axis

            rbomax=0.d0
            rbomin=rm
            zbomax=zm
            zbomin=zm

         do ib=1,nbtab
            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
            if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
            if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)
	   enddo

           rc0new=0.5d0*(rbomax+rbomin)
           zc0new=0.5d0*(zbomax+zbomin)

           drc0=rc0new-rc0
           dzc0=zc0new-zc0

            rm=rm+drc0
            zm=zm+dzc0

        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !call arc_x_bnd(nt,nbtab,rbtab,zbtab)  !!!?!!!!

            !if(kstep.eq.5) call bnd_mov 
             !call bnd_mov 

            rbomax=0.d0
            rbomin=rm
            zbomax=zm
            zbomin=zm

         do ib=1,nbtab
          !  write(6,*) 'ib:',ib
           drx=rbtab(ib)-rm
           dzx=zbtab(ib)-zm

           robn(ib)=sqrt(drx**2+dzx**2)
            tetp=dacos(drx/robn(ib))
           if(dzx.lt.0.d0) then
            tetbn(ib)=-tetp
           else
            tetbn(ib)=tetp
          endif
            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
            if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
            if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)
	   enddo

           rc0=0.5d0*(rbomax+rbomin)
           zc0=0.5d0*(zbomax+zbomin)
              ! write(*,*) 'grid1:rc0',rc0

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo
         do ib=2,nbtab
          if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

          nbn=nbtab
        if(tetbn(1)+2.d0*pi-tetbn(nbn) .GT. 1.d-4) then
          nbn=nbn+1
          tetbn(nbn)=tetbn(1)+2.d0*pi
          robn(nbn)=robn(1)
        endif

         do j=1,nbn

            teta(j+1)=tetbn(j)
            ro(iplas,j+1)=robn(j)
	      ro(1,j)=0.d0

	   enddo

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi
	       ro(iplas,nt)=ro(iplas,2)
	       ro(iplas,1)=ro(iplas,nt1)

         do j=1,Nt
              tetp=teta(j)
	       ro(1,j)=0.d0
	       r(iplas,j)=rm+ro(iplas,j)*dcos(tetp)
	       z(iplas,j)=zm+ro(iplas,j)*dsin(tetp)

	   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           do i=2,iplas1
            do j=2,nt

	       ro(i,j)=ro(iplas,j)*ronor(i,j)

            enddo
           enddo

           do i=1,iplas1
            do j=2,nt1

	       r(i,j)=rm+ro(i,j)*dcos(teta(j))
	       z(i,j)=zm+ro(i,j)*dsin(teta(j))
	       psin(i,j)=psia(i)

            enddo
           enddo

            do 120 j=1,nt

             r(1,j)=rm
             z(1,j)=zm
             ro(1,j)=0.d0
             psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)
             psin(i,1)=psin(i,nt1)

             ro(i,nt)=ro(i,2)
             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             psin(i,nt)=psin(i,2)

 333    continue

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo

               !call wrb

	!pause 'grid_1:wrb'
              ! stop
              
      !deallocate( rbtab, zbtab )

          return
          end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_11(igdf,nstep)
         
         use bnd_modul       

           include 'double.inc'
           include 'dim.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           include 'compol.inc'
           !common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

           real*8 robn(nbtabp),tetbn(nbtabp)
           real*8 RRK(nb4),CCK(nb4),WRK(nb6),CWK(4)

!            cos(xx)=dcos(xx)
!            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)
            
!      if(.not.allocated(rbtab))
!     %allocate( rbtab(nbtab), zbtab(nbtab) )

!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname)
!         !open(1,file='tab_bnd.dat')
!	        read(1,*) nbtab
!	        do ib=1,nbtab
!	           read(1,*) rbtab(ib),zbtab(ib)
!	        enddo
!         close(1)

         do ib=1,nbtab
          !  write(6,*) 'ib:',ib
           drx=rbtab(ib)-rm
           dzx=zbtab(ib)-zm

           robn(ib)=sqrt(drx**2+dzx**2)
            tetp=dacos(drx/robn(ib))
           if(dzx.lt.0.d0) then
            tetbn(ib)=-tetp
           else
            tetbn(ib)=tetp
          endif
            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
	   enddo

           rc0=0.5d0*(rbomax+rbomin)
               !write(6,*) 'grid1:rc0',rc0

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo
         do ib=2,nbtab
          if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

          nbn=nbtab
        if(tetbn(1)+2.d0*pi-tetbn(nbn) .GT. 1.d-4) then
          nbn=nbn+1
          tetbn(nbn)=tetbn(1)+2.d0*pi
          robn(nbn)=robn(1)
        endif

        CALL E01BAF(Nbn,tetbn,robn,RRK,CCK,nbn+4,WRK,6*nbn+16,IFAIL)

         do j=1,Nt
           tetp=teta(j)
           if(tetp.lt.tetbn(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetbn(nbn)) then
             tetp=tetp-2.d0*pi
           endif
           if(tetp.lt.tetbn(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetbn(nbn)) then
             tetp=tetp-2.d0*pi
           endif

          CALL E02BCF(Nbn+4,RRK,CCK,tetp     ,0,CWk,IFAIL)

            ro(iplas,j)=cwk(1)
	      ro(1,j)=0.d0

	       r(iplas,j)=rm+ro(iplas,j)*dcos(tetp)
	       z(iplas,j)=zm+ro(iplas,j)*dsin(tetp)

	   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           do i=2,iplas1
            do j=2,nt

	       ro(i,j)=ro(iplas,j)*ronor(i,j)

            enddo
           enddo

           do i=1,iplas1
            do j=2,nt1

	       r(i,j)=rm+ro(i,j)*dcos(teta(j))
	       z(i,j)=zm+ro(i,j)*dsin(teta(j))
	       psin(i,j)=psia(i)

            enddo
           enddo

            do 120 j=1,nt

             r(1,j)=rm
             z(1,j)=zm
             ro(1,j)=0.d0
             psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)
             psin(i,1)=psin(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             psin(i,nt)=psin(i,2)

 333    continue

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo

      !        call wrb

	!pause 'grid_1:wrb'
              ! stop

      !deallocate( rbtab, zbtab )
      
          return
          end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine grid_p0(igdf,nstep)
           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'
            erru=1.d0	    
           do j=1,nt
           psi(iplas,j)=0.d0 
           enddo 
          return
          end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine grid_p1(igdf,nstep)
           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

           do j=1,nt
           psi(iplas,j)=0.d0 
           enddo 

          return
          end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_0(igdf,nstep)
         
         use bnd_modul       

           include 'double.inc'
           include 'dim.inc'

          !INCLUDE 'parrec.inc'

          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           include 'compol.inc'
          ! common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

	common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh

           !include 'comrec.inc'

           real*8 drds(ntp),dzds(ntp),rcrv(ntp),rcrn(ntp)
           real*8 dls(ntp),drodt(ntp)
           real*8 ron(ntp),tetn(ntp)
           real*8 robn(nbtabp),tetbn(nbtabp)
           real*8 RRK(nb4),CCK(nb4),WRK(nb6),CWK(4)

	common
     *  /c_kpr/kpr

            frbon(r0,aa,tr,tet)=r0+(r0/aa)*dcos(tet+tr*dsin(tet))
            fzbon(r0,z0,aa,el,tet)=z0+(r0/aa)*el*dsin(tet)

!            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)
            
!       if(.not.allocated(rbtab))
!     %allocate( rbtab(nbtab), zbtab(nbtab) )
           

      if(nstep.eq.0) then

         if(igdf.eq.0)then

	        do i=1,iplas

           psia(i)=(iplas-i)/(iplas-1.d0)

              enddo

	        do i=1,iplas1

	     dpsda(i)=-1.d0

              enddo

         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3)then

	        do i=1,iplas
           psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
              enddo

	        do i=1,iplas1
           dpsda(i)=(1-2*i)/(iplas-1.d0)
              enddo

         endif

         endif


               rm=rm0    !
               zm=zm0    !  position of magn. axis

               psim=1.d0
               psip=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(nbsh.eq.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             dtet=2.d0*pi/(nt-2)
             teta(1)=-dtet

            do j=2,nt
             teta(j)=teta(j-1)+dtet
            enddo

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

           do j=1,nt

            tet=teta(j)

	   tri=0.5d0*(tr_up+tr_lw+(tr_up-tr_lw)*sin(tet))
	   ell=0.5d0*(el_up+el_lw+(el_up-el_lw)*sin(tet))

            rrr=frbon(rc0,asp0,tri,tet)
            r(iplas,j)=rrr
            zzz=fzbon(rc0,zc0,asp0,ell,tet)
            z(iplas,j)=zzz

	      roxx=sqrt((r(iplas,j)-rm)**2+(z(iplas,j)-zm)**2)
	      ro(iplas,j)=roxx
	      ro(1,j)=0.d0

           enddo

         !open(1,file='tab_bnd.dat') 
	      nbtab=nt1
	     !write(1,*) nbtab 
	      do ib=1,nbtab
	         rbtab(ib)= r(iplas,ib)
	         zbtab(ib)= z(iplas,ib)
	        !write(1,*) rbtab(ib),zbtab(ib) 
	      enddo
         !close(1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname)
!         !open(1,file='tab_bnd.dat')
!	     read(1,*) nbtab
!	     do ib=1,nbtab
!	        read(1,*) rbtab(ib),zbtab(ib)
!	     enddo
!         close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
c_ast!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !call arc_x_bnd(nt,nbtab,rbtab,zbtab)
            call arc_x_bnd(nt)

!         if(kpr.ge.0) then
!
!        write(fname,'(a,a)') path(1:kname),'bon_arc.dat'
!        open(1,file=fname)
!         !open(1,file='bon_arc.dat') 
!	     write(1,*) nbtab-1 
!	      do ib=1,nbtab-1
!	        write(1,*) rbtab(ib),zbtab(ib) 
!	      enddo
!         close(1) 
!
!           endif

            rbomax=0.d0
            rbomin=rm
            zbomax=zm
            zbomin=zm

         do ib=1,nbtab

           drx=rbtab(ib)-rm
           dzx=zbtab(ib)-zm

           robn(ib)=sqrt(drx**2+dzx**2)

            tetp=dacos(drx/robn(ib))

           if(dzx.lt.0.d0) then
            !tetbn(ib)=2.d0*pi-tetp
            tetbn(ib)=-tetp
           else
            tetbn(ib)=tetp
          endif

            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)
            if(zbtab(ib).gt.zbomax) zbomax=zbtab(ib)
            if(zbtab(ib).lt.zbomin) zbomin=zbtab(ib)

        enddo

           rc0=0.5d0*(rbomax+rbomin)
           zc0=0.5d0*(zbomax+zbomin)

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

          nbn=nbtab

        if(tetbn(1)+2.d0*pi-tetbn(nbn) .GT. 1.d-4) then

          nbn=nbn+1

          tetbn(nbn)=tetbn(1)+2.d0*pi
          robn(nbn)=robn(1)

        endif

         do j=1,nbn

            teta(j+1)=tetbn(j)
            ro(iplas,j+1)=robn(j)
	      ro(1,j)=0.d0

	   enddo

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi

         do j=2,Nt1
              tetp=teta(j)
	       r(iplas,j)=rm+ro(iplas,j)*dcos(tetp)
	       z(iplas,j)=zm+ro(iplas,j)*dsin(tetp)

	   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!            do j=2,nt  !<---
!
!         if(teta(j).lt.teta(j-1)) then
!            teta(j)=teta(j)+pi
!         endif
!
!         if(teta(j).lt.teta(j-1)) then
!            teta(j)=teta(j)+pi
!         endif
!
!	      enddo      !<---

            do i=2,iplas1
            do j=2,nt
	         ro(i,j)=ro(iplas,j)*dsqrt(1.d0-psia(i))
            enddo
            enddo

           do i=1,iplas1
           do j=2,nt1
	       r(i,j)=rm+ro(i,j)*dcos(teta(j))
	       z(i,j)=zm+ro(i,j)*dsin(teta(j))
	       psin(i,j)=psia(i)
           enddo
           enddo

!         call wrb

            do 120 j=1,nt

              r(1,j)=rm
              z(1,j)=zm
              ro(1,j)=0.d0
              psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)
             psin(i,1)=psin(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             psin(i,nt)=psin(i,2)

 333    continue

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo

         call wrb

       !deallocate( rbtab, zbtab )
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_00(igdf,nstep)

         use bnd_modul
                
           include 'double.inc'
           include 'dim.inc'

          !INCLUDE 'parrec.inc'

          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           include 'compol.inc'
           !common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

	common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh

           !include 'comrec.inc'

           real*8 drds(ntp),dzds(ntp),rcrv(ntp),rcrn(ntp)
           real*8 dls(ntp),drodt(ntp)
           real*8 ron(ntp),tetn(ntp)
           real*8 robn(nbtabp),tetbn(nbtabp)
           real*8 RRK(nb4),CCK(nb4),WRK(nb6),CWK(4)

            frbon(r0,aa,tr,tet)=r0+(r0/aa)*dcos(tet+tr*dsin(tet))
            fzbon(r0,z0,aa,el,tet)=z0+(r0/aa)*el*dsin(tet)

!            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)

!      if(.not.allocated(rbtab))
!     %allocate( rbtab(nbtab), zbtab(nbtab) )

      if(nstep.eq.0) then

         if(igdf.eq.0)then

	        do i=1,iplas

           psia(i)=(iplas-i)/(iplas-1.d0)

              enddo

	        do i=1,iplas1

	     dpsda(i)=-1.d0

              enddo

         elseif(igdf.eq.1 .OR. igdf.eq.2)then

	        do i=1,iplas
           psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
              enddo

	        do i=1,iplas1
           dpsda(i)=(1-2*i)/(iplas-1.d0)
              enddo

         endif

         endif

             dtet=2.d0*pi/(nt-2)
             teta(1)=-dtet

            do j=2,nt
             teta(j)=teta(j-1)+dtet
            enddo

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

               rm=rm0    !
               zm=zm0    !  position of magn. axis

               psim=1.d0
               psip=0.d0

        if(nbsh.eq.0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           do j=1,nt

            tet=teta(j)

	   tri=0.5d0*(tr_up+tr_lw+(tr_up-tr_lw)*sin(tet))
	   ell=0.5d0*(el_up+el_lw+(el_up-el_lw)*sin(tet))

            rrr=frbon(rc0,asp0,tri,tet)
            r(iplas,j)=rrr
            zzz=fzbon(rc0,zc0,asp0,ell,tet)
            z(iplas,j)=zzz

	      roxx=sqrt((r(iplas,j)-rm)**2+(z(iplas,j)-zm)**2)
	      ro(iplas,j)=roxx
	      ro(1,j)=0.d0

           enddo

          !open(1,file='tab_bnd.dat')
	      nbtab=nt1
	     ! write(1,*) nbtab
	      do ib=1,nbtab
	      !   write(1,*) r(iplas,ib),z(iplas,ib)
	      enddo
         ! close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
!        open(1,file=fname)
!         !open(1,file='tab_bnd.dat')
!	     read(1,*) nbtab
!	     do ib=1,nbtab
!	        read(1,*) rbtab(ib),zbtab(ib)
!	     enddo
!         close(1)

            rbomax=0.d0
            rbomin=rm

         do ib=1,nbtab

          !  write(6,*) 'ib:',ib

           drx=rbtab(ib)-rm
           dzx=zbtab(ib)-zm

           robn(ib)=sqrt(drx**2+dzx**2)

            tetp=dacos(drx/robn(ib))

           if(dzx.lt.0.d0) then
            tetbn(ib)=-tetp
           else
            tetbn(ib)=tetp
          endif

            if(rbtab(ib).gt.rbomax) rbomax=rbtab(ib)
            if(rbtab(ib).lt.rbomin) rbomin=rbtab(ib)

        enddo

           rc0=0.5d0*(rbomax+rbomin)

C           write(6,*) 'Subr."grid_0": rc0 = ', rc0

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

         do ib=2,nbtab
	    if(tetbn(ib).lt.tetbn(ib-1)) tetbn(ib)=tetbn(ib)+2.d0*pi
	   enddo

          nbn=nbtab

        if(tetbn(1)+2.d0*pi-tetbn(nbn) .GT. 1.d-4) then

          nbn=nbn+1

          tetbn(nbn)=tetbn(1)+2.d0*pi
          robn(nbn)=robn(1)

        endif

           !write(6,*) tetbn
           !pause 'pause'

        CALL E01BAF(Nbn,tetbn,robn,RRK,CCK,nbn+4,WRK,6*nbn+16,IFAIL)

         !write(6,*) '*ifail=',ifail

         do j=1,Nt

           tetp=teta(j)

           if(tetp.lt.tetbn(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetbn(nbn)) then
             tetp=tetp-2.d0*pi
           endif

           if(tetp.lt.tetbn(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetbn(nbn)) then
             tetp=tetp-2.d0*pi
           endif

          CALL E02BCF(Nbn+4,RRK,CCK,tetp     ,0,CWk,IFAIL)

            ro(iplas,j)=cwk(1)
	      ro(1,j)=0.d0

	       r(iplas,j)=rm+ro(iplas,j)*dcos(tetp)
	       z(iplas,j)=zm+ro(iplas,j)*dsin(tetp)

       !write(6,*) 'ifail=',ifail,j

	   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do j=1,nt
         drx=r(iplas,j)-rm
         dzx=z(iplas,j)-zm
         drob=sqrt(drx**2+dzx**2)
         tetp=acos(drx/drob)
           if(dzx.lt.0.d0) then
            teta(j)=-tetp
           else
            teta(j)=tetp
          endif
        enddo

            do j=2,nt  !<---

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

         if(teta(j).lt.teta(j-1)) then
            teta(j)=teta(j)+2.d0*pi
         endif

	      enddo      !<---

           teta(1)=teta(nt1)-2.d0*pi
           teta(nt)=teta(2)+2.d0*pi


            do i=2,iplas1
            do j=2,nt
	         ro(i,j)=ro(iplas,j)*dsqrt(1.d0-psia(i))
            enddo
            enddo

           do i=1,iplas1
           do j=2,nt1
	       r(i,j)=rm+ro(i,j)*dcos(teta(j))
	       z(i,j)=zm+ro(i,j)*dsin(teta(j))
	       psin(i,j)=psia(i)
           enddo
           enddo

            do 120 j=1,nt

              r(1,j)=rm
              z(1,j)=zm
              ro(1,j)=0.d0
              psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)
             psin(i,1)=psin(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             psin(i,nt)=psin(i,2)

 333    continue

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo

	!do i=1,iplas
	!dpdpsi(i)=tabp(psia(i))
	!dfdpsi(i)=tabf(psia(i))
	!enddo

          call wrb

      !deallocate( rbtab, zbtab )
     
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine loop95_b(tetpol,ntet,ro0,u0)

         include 'double.inc'
         include 'parrc1.inc'
         include 'comrec.inc'
             common/limpotcom/  zlimpotential(nip,njp)

         dimension rxb(nbndp2),zxb(nbndp2)
         real*8 ut(nip,njp)
         real*8 roxb(nbndp2),tetxb(nbndp2)
         real*8 RRK(nbndp4),CCK(nbndp4),WRK(nbndp6)
         real*8 CWK(4),tetpol(*),ro0(*)
	real*8 utemps(nip,njp)
	common
     *  /c_kpr/kpr

         xzer(x1,x2,v1,v2) = (x1*v2-x2*v1)/(v2-v1)

C         write(6,*) '***loop95:enter'

        ! write(6,*) 'loop95:imax,jmax',imax,jmax
        ! write(6,*) 'loop95:rmax,zmax',rm,zm

          pi=3.14159265359d0

	 !if(u0.lt.2.d-4) u0=2.d-4

         imax1=imax-1
         jmax1=jmax-1

         do 10 i=imax1,imax
          if(xm.le.x(i+1) .AND. xm.gt.x(i)) icell=i
 10      continue

         do 20 j=jmax1,jmax
          if(ym.le.y(j+1) .AND. ym.gt.y(j)) jcell=j
 20      continue

        ! write(6,*) 'loop95:icell,jcell',icell,jcell

         i=icell
         j=jcell+1

         do 15 ii=1,ni
         do 15 jj=1,nj
          ut(ii,jj)=un(ii,jj)-u0
	utemps(ii,jj)=ut(ii,jj)

	if (zlimpotential(ii,jj).lt.0.)	
     & ut(ii,jj)=abs(ut(ii,jj))*zlimpotential(ii,jj)
 15      continue

         ig=1

 888     continue

                do 885 i=imax,ni1
          !   write(6,*) 'un(i,j)',un(i,j),un(i+1,j)
             if(ut(i,j)*ut(i+1,j).le.0.d0) then

              ic=i
              jc=j

        rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
        zxb(ig)=y(j)

              go to 886

             endif

 885            continue

          if(kpr.eq.1) then
          write(6,*) 'loop95: first point was not find'
          endif
            stop

 886            continue

              ic1=ic
              jc1=jc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          lin=1

 100      ig=ig+1

          i=ic
          j=jc

             if(ut(i+1,j)*ut(i,j).le.0.d0 .AND. lin.ne.1) then

              ic=i
              jc=j-1

              lin=3

              rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
              zxb(ig)=y(j)

             elseif(ut(i+1,j+1)*ut(i+1,j).le.0.d0.AND. lin.ne.2) then

              ic=i+1
              jc=j

              lin=4

              rxb(ig)=x(i+1)
              zxb(ig)=xzer(y(j+1),y(j),ut(i+1,j+1),ut(i+1,j))

             elseif(ut(i+1,j+1)*ut(i,j+1).le.0.d0 .AND. lin.ne.3) then

              ic=i
              jc=j+1

              lin=1

              rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j+1),ut(i,j+1))
              zxb(ig)=y(j+1)

             elseif(ut(i,j+1)*ut(i,j).le.0.d0.AND. lin.ne.4) then

              ic=i-1
              jc=j

              lin=2

              rxb(ig)=x(i)
              zxb(ig)=xzer(y(j+1),y(j),ut(i,j+1),ut(i,j))

             endif

c        write(6,*) 'loop:ic,jc',ic,jc,ig

            if(jc.eq.jc1 .AND. ic.eq.ic1) then

              ig=ig+1

              rxb(ig)=rxb(2)
              zxb(ig)=zxb(2)

              go to 1212

            endif

             go to 100

 1212        nxb=ig
	ut=utemps

            do 200 ig=1,nxb

         drx=rxb(ig)-xm
         dzx=zxb(ig)-ym

         tetp=datan(dzx/drx)

       if(ig.ne.1) then
         if(drx.lt.0.d0) tetp=tetp+pi
         deltet=tetp-tetxb(ig-1)
         if(deltet.lt.0.d0) tetp=tetp+2.d0*pi
       endif

         tetxb(ig)=tetp
         roxb(ig)=dsqrt(drx**2+dzx**2)

c         write(6,*) 'ig,ro,tet',ig,roxb(ig),tetxb(ig)

 200        continue

        CALL E01BAF(Nxb,tetxb,roxb,RRK,CCK,nxb+4,WRK,6*nxb+16,IFAIL)

c       write(6,*) 'ifail=',ifail

        DO 210 j=1,Ntet
           tetp=tetpol(j)
           if(tetp.lt.tetxb(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetxb(nxb)) then
             tetp=tetp-2.d0*pi
           endif
          CALL E02BCF(Nxb+4,RRK,CCK,tetp     ,0,CWk,IFAIL)
            ro0(j)=cwk(1)
c       write(6,*) 'ifail=',ifail,j
 210    CONTINUE

           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine loopL_b(tetpol,ntet,ro0,u0)

         include'double.inc'
         include'parrc1.inc'
         !include'comrec.inc'
         common /comind/ ni,nj,ni1,nj1,ni2,nj2,nbnd,nkin,nkout
         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                   ui(nip,njp),g(nip,njp),
     +                   ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon,ue_k(nip,njp,500),rcondzzz
	common/limpotcom/  zlimpotential(nip,njp)
		 

         dimension rxb(nbndp2),zxb(nbndp2)

         real*8 ut(nip,njp)

         real*8 roxb(nbndp2),tetxb(nbndp2)

         real*8 tetpol(*),ro0(*)
	real*8 utemps(nip,njp)

	common
     *  /c_kpr/kpr
         xzer(x1,x2,v1,v2) = (x1*v2-x2*v1)/(v2-v1)

C           write(6,*) '*** loopL_b: enter'

        ! write(6,*) 'loopL:imax,jmax',imax,jmax
        !write(6,*) 'loopL:rmax,zmax',rm,zm

          pi=3.14159265359d0

	 !if(u0.lt.2.d-4) u0=2.d-4

         imax1=imax-1
         jmax1=jmax-1

         do 10 i=imax1,imax
          if(xm.le.x(i+1) .AND. xm.gt.x(i)) icell=i
 10      continue

         do 20 j=jmax1,jmax
          if(ym.le.y(j+1) .AND. ym.gt.y(j)) jcell=j
 20      continue

        ! write(6,*) 'loopL:icell,jcell',icell,jcell

         i=icell
         j=jcell+1

         do 15 ii=1,ni
         do 15 jj=1,nj
          ut(ii,jj)=un(ii,jj)-u0
	utemps(ii,jj)=ut(ii,jj)

	if (zlimpotential(ii,jj).lt.0.)	
     & ut(ii,jj)=abs(ut(ii,jj))*zlimpotential(ii,jj)
 15      continue


!	open(32,file='psimap.dat')
!	write(32,*) x(1:ni),y(1:nj),ut(1:ni,1:nj)
!	close(32)

         ig=1

 888         continue

                do 885 i=imax,ni1

          !   write(6,*) 'un(i,j)',un(i,j),un(i+1,j)

             if(ut(i,j)*ut(i+1,j).le.0.d0) then

              ic=i
              jc=j

        rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
        zxb(ig)=y(j)

              go to 886
             endif

 885            continue
          if(kpr.eq.1) then
          write(*,*) 'loopL: first point was not find'
          endif
            stop

 886            continue

              ic1=ic
              jc1=jc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          lin=1

 100      ig=ig+1

          i=ic
          j=jc

             if(ut(i+1,j)*ut(i,j).le.0.d0 .AND. lin.ne.1) then

              ic=i
              jc=j-1

              lin=3

              rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j),ut(i,j))
              zxb(ig)=y(j)

             elseif(ut(i+1,j+1)*ut(i+1,j).le.0.d0.AND. lin.ne.2) then

              ic=i+1
              jc=j

              lin=4

              rxb(ig)=x(i+1)
              zxb(ig)=xzer(y(j+1),y(j),ut(i+1,j+1),ut(i+1,j))

             elseif(ut(i+1,j+1)*ut(i,j+1).le.0.d0 .AND. lin.ne.3) then

              ic=i
              jc=j+1

              lin=1

              rxb(ig)=xzer(x(i+1),x(i),ut(i+1,j+1),ut(i,j+1))
              zxb(ig)=y(j+1)

             elseif(ut(i,j+1)*ut(i,j).le.0.d0.AND. lin.ne.4) then

              ic=i-1
              jc=j

              lin=2

              rxb(ig)=x(i)
              zxb(ig)=xzer(y(j+1),y(j),ut(i,j+1),ut(i,j))

             endif


            if(jc.eq.jc1 .AND. ic.eq.ic1) then

              ig=ig+1

              rxb(ig)=rxb(2)
              zxb(ig)=zxb(2)

              go to 1212

            endif

!	write(687,*) ig,i,j,ic,jc,ut(i,j),u0,un(i,j)
	if (ig.eq.2056) then
		write(*,*) 'cwap boundary'
		call err_catch_a
	endif
	
             go to 100

 1212        nxb=ig
	ut=utemps
            do 200 ig=1,nxb

         drx=rxb(ig)-xm
         dzx=zxb(ig)-ym

         tetp=datan(dzx/drx)

       if(ig.ne.1) then
         if(drx.lt.0.d0) tetp=tetp+pi
         deltet=tetp-tetxb(ig-1)
         if(deltet.lt.0.d0) tetp=tetp+2.d0*pi
       endif

         tetxb(ig)=tetp

         roxb(ig)=dsqrt(drx**2+dzx**2)

c         write(6,*) 'ig,ro,tet',ig,roxb(ig),tetxb(ig)

 200        continue

c       write(6,*) 'ifail=',ifail

        DO 210 j=1,Ntet
           tetp=tetpol(j)
           if(tetp.lt.tetxb(1)) then
             tetp=tetp+2.d0*pi
           elseif(tetp.gt.tetxb(nxb)) then
             tetp=tetp-2.d0*pi
           endif

            do jw=1,nxb-1
             tetmn=tetxb(jw)
             tetpl=tetxb(jw+1)
             romn=roxb(jw)
             ropl=roxb(jw+1)
           if(tetp.le.tetpl .AnD. tetp.ge.tetmn) go to 347
            enddo

 347       continue

           ro0(j)=(ropl*(tetp-tetmn)+romn*(tetpl-tetp))/(tetpl-tetmn)

 210    CONTINUE


           return
           end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_b(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!          include'comrec.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                  ui(nip,njp),g(nip,njp),
     +                  ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

           real*8 ron(ntp),ronm(ntp)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)

                !    	write(6,*) 'grid_b:eter'
                !        write(6,*) 'gridb:igdf,nstep',igdf,nstep
                !        call wrb
                !        pause 'gridb:enter:wrb'

            if(nstep.eq.0) then

         if(igdf.eq.0)then

	        do i=1,iplas
           psia(i)=(iplas-i)/(iplas-1.d0)
              enddo

	        do i=1,iplas1
	     dpsda(i)=-1.d0
              enddo

         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3) then

	        do i=1,iplas
           psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
              enddo

	        do i=1,iplas1
           dpsda(i)=(1-2*i)/(iplas-1.d0)
              enddo

         endif
         endif

              rm=xm    !
              zm=ym    !  position of magn. axis
                       !

              psim=um-up
              psip=0.d0

                dtet=2.d0*pi/(nt-2)
                teta(1)=-dtet

            do 30 j=2,nt
             teta(j)=teta(j-1)+dtet
 30         continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

            psia(1)=1.d0

          nsh=1

          xs(nsh)=x(imax)
          ys(nsh)=y(jmax)
         fun(nsh)=u(imax,jmax)

        do k=-1,1
          ii= imax+k
        do l=-1,1
          jj= jmax+l

          if(ii.ne.imax .OR. jj.ne.jmax) then

          nsh=nsh+1
          xs(nsh)=x(ii)
          ys(nsh)=y(jj)
         fun(nsh)=u(ii,jj)

	    endif

        enddo
        enddo

          call deriv5(xs,ys,fun,nsh,5,dp)

           !write(6,*) dp

	     stpx=(x(imax+1)-x(imax))
	     stpy=(y(jmax+1)-y(jmax))*0.75d0

            do i=2,iplas

             u0=psia(i)
           ur0=up+u0*(um-up)

             do j=1,nt

         ttj=teta(j)

        ro2j=(ur0-um)/( 0.5d0*dp(3)*cos(ttj)**2
     +                +     dp(4)*cos(ttj)*sin(ttj)
     +                + 0.5d0*dp(5)*sin(ttj)**2 )

          ron(j)=sqrt(ro2j)
          ronm(j)=ron(j)    

             r(i,j)=rm+ron(j)*cos(teta(j))
             z(i,j)=zm+ron(j)*sin(teta(j))

             ro(i,j)=ron(j)
             psin(i,j)=u0

             enddo  ! j-loop

	     ibeg=i+1

C	 write(6,*) 'ro(1)',ro(i,1),ro(i,2),stpx,stpy

	   if( ron(2) .GT. dmax1(stpx,stpy) ) go to 2799

             enddo   !i-loop

 2799     continue

                ! pause 'pausew'

            do 100 i=iplas,ibeg,-1

             u0=psia(i)

            call loopL_b(teta,nt,ron,u0)

C          write(6,*) 'loopL i=',i,u0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             do j=1,nt

              if(ron(j) .LT. ronm(j)) then
               ron(j)=(ro(i+1,j)+ronm(j)*(i-ibeg+1.d0))/(i-ibeg+2.d0)
              endif

             enddo						  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do 110 j=1,nt

             r(i,j)=rm+ron(j)*cos(teta(j))
             z(i,j)=zm+ron(j)*sin(teta(j))

             ro(i,j)=ron(j)
             psin(i,j)=u0

 110        continue
 100        continue

            do 120 j=1,nt

             r(1,j)=rm
             z(1,j)=zm
             ro(1,j)=0.d0
             psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             ro(i,nt)=ro(i,2)

 333    continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo

          return
          end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_b_ef(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!          include'comrec.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                  ui(nip,njp),g(nip,njp),
     +                  ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

           real*8 ron(ntp),ronm(ntp),xstepp,rggr(iplas)
           real*8 psitemp(iplas),uuuu(9),xxx(9,2)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)
           real*8 utemp(nip,njp),um_t,up_t
           real*8 xm_t,ym_t,up_new
           real*8 ucaz(500)
	integer imax_t,jmax_t
	
            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)

                !    	write(6,*) 'grid_b:eter'
                !        write(6,*) 'gridb:igdf,nstep',igdf,nstep
                !        call wrb
                !        pause 'gridb:enter:wrb'

!            if(nstep.eq.0) then

!	write(*,*) igdf
	
!         if(igdf.eq.0)then

!	        do i=1,iplas
!           psia(i)=(iplas-i)/(iplas-1.d0)
!              enddo

!	        do i=1,iplas1
!	     dpsda(i)=-1.d0
!              enddo

!         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3) then

	        do i=1,iplas
           psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
              enddo

	        do i=1,iplas1
           dpsda(i)=(1-2*i)/(iplas-1.d0)
              enddo
	utemp=u
	up_t=up
	xm_t=xm
	ym_t=ym
	um_t=um


!	write(9754,*) xm,ym,um,xm_t,ym_t,um_t
!	write(9754,*) imax,jmax,imax_t,jmax_t
!	write(9754,*) up,up_new,up_t
!	write(*,*) up_t
              rm=xm_t    !
              zm=ym_t    !  position of magn. axis
                       !




              psim=um_t-up_t
              psip=0.d0

                dtet=2.d0*pi/(nt-2)
                teta(1)=-dtet

            do 30 j=2,nt
             teta(j)=teta(j-1)+dtet
 30         continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

            psia(1)=1.d0


! first find boundary
	i=iplas
	xstepp=0.1*sqrt((x(2)-x(1))**2.+(y(2)-y(1))**2.)
	do j=1,nt
		do k=1,nip*njp
						ronzz=xstepp*k
             dumx=rm+ronzz*cos(teta(j))
             dumy=zm+ronzz*sin(teta(j))
			j1=((dumx-x(1))/(x(2)-x(1)))+1	
			j2=((dumy-y(1))/(y(2)-y(1)))+1	
		uuuu(1)=utemp(j1-1,j2-1)
		uuuu(2)=utemp(j1,j2-1)
		uuuu(3)=utemp(j1-1,j2)
		uuuu(4)=utemp(j1,j2)
		uuuu(5)=utemp(j1-1,j2+1)
		uuuu(6)=utemp(j1,j2+1)
		uuuu(7)=utemp(j1+1,j2+1)
		uuuu(8)=utemp(j1+1,j2)
		uuuu(9)=utemp(j1+1,j2-1)

	call biquadratic_interp_spidef(x(j1-1),x(j1),x(j1+1),
     & y(j2-1),y(j2),y(j2+1),uuuu,
     &  dumx,dumy,
     &  u002)
		 ucaz(k)=u002
		if (ucaz(k).eq.up_t) then
			ro(i,j)=ronzz
			goto 5561
		endif
		if (ucaz(k).lt.up_t) then
			ronzz2=xstepp*(k-1)
			ro(i,j)=ronzz2+(up_t-ucaz(k-1))/(ucaz(k-1)-ucaz(k))*(ronzz2-ronzz)
			goto 5561
		endif
		enddo

5561	continue

	psitemp(1)=1.
	psitemp(iplas)=0.
	rggr(iplas)=ro(iplas,j)
	rggr(1)=0.
!	write(332,'(6444E25.11)') rm,zm,um_t,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)
	do i=2,iplas-1
		rggr(i)=ro(iplas,j)*(i-1.)/(iplas-1.)
		dumx=rm+rggr(i)*cos(teta(j))
		dumy=zm+rggr(i)*sin(teta(j))
		j1=((dumx-x(1))/(x(2)-x(1)))+1	
		j2=((dumy-y(1))/(y(2)-y(1)))+1	

		xxx(1,1)=x(j1-1)
		xxx(1,2)=y(j2-1)
		uuuu(1)=utemp(j1-1,j2-1)
		xxx(2,1)=x(j1)
		xxx(2,2)=y(j2-1)
		uuuu(2)=utemp(j1,j2-1)
		xxx(3,1)=x(j1-1)
		xxx(3,2)=y(j2)
		uuuu(3)=utemp(j1-1,j2)
		xxx(4,1)=x(j1)
		xxx(4,2)=y(j2)
		uuuu(4)=utemp(j1,j2)
		xxx(5,1)=x(j1-1)
		xxx(5,2)=y(j2+1)
		uuuu(5)=utemp(j1-1,j2+1)
		xxx(6,1)=x(j1)
		xxx(6,2)=y(j2+1)
		uuuu(6)=utemp(j1,j2+1)
		xxx(7,1)=x(j1+1)
		xxx(7,2)=y(j2+1)
		uuuu(7)=utemp(j1+1,j2+1)
		xxx(8,1)=x(j1+1)
		xxx(8,2)=y(j2)
		uuuu(8)=utemp(j1+1,j2)
		xxx(9,1)=x(j1+1)
		xxx(9,2)=y(j2-1)
		uuuu(9)=utemp(j1+1,j2-1)

	call biquadratic_interp_spidef(x(j1-1),x(j1),x(j1+1),
     & y(j2-1),y(j2),y(j2+1),uuuu,
     &  dumx,dumy,
     &  u002)
!	write(*,*) i,j,u002
!		u002=bqlin_(dumx,dumy,xxx,uuuu)
!		u002=blin_(dumx,dumy,xxx(4,1),xxx(7,1),
!     & xxx(4,2),xxx(7,2),uuuu(4),uuuu(8),
!     & uuuu(7),uuuu(6))
			psitemp(i)=(u002-up_t)/(um_t-up_t)
!	write(332,'(6444E25.11)') dumx,dumy,u002,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)
	enddo
!	write(332,'(6444E25.11)') rm+rggr(iplas)*cos(teta(j)),
!     & zm+rggr(iplas)*sin(teta(j)),up_t,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)


!	write(986,*) j
		call reinterp_back_quad((1.-psitemp(1:iplas))**0.5,rggr(1:iplas),
     & iplas, 
     &   (1.-psia(1:iplas))**0.5,ro(1:iplas,j),iplas)

!	write(7489,*) j,(1.-psitemp(1:iplas))**0.5,rggr(1:iplas),
!     & (1.-psia(1:iplas))**0.5,ro(1:iplas,j)

!	if (j.eq.1) then
!			write(*,*) 'j',j,psitemp,rggr,psia(1:iplas),ro(1:iplas,j)
!	endif
	
	enddo !j cycle



            do i=2,iplas

             u0=psia(i)
           ur0=up_t+u0*(um_t-up_t)

             do j=1,nt

             r(i,j)=rm+ro(i,j)*cos(teta(j))
             z(i,j)=zm+ro(i,j)*sin(teta(j))

             psin(i,j)=u0

             enddo  ! j-loop

             enddo   !i-loop

            do  j=1,nt

             r(1,j)=rm
             z(1,j)=zm
             ro(1,j)=0.d0
             psin(1,j)=1.d0

						enddo

        do  i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             ro(i,nt)=ro(i,2)

				enddo

             do j=1,nt
               ronor(iplas,j)=1.d0
             enddo				  

             do j=1,nt 
             do i=1,iplas

               ronor(i,j)=ro(i,j)/ro(iplas,j)

               psin(i,j)=psia(i)

             enddo				  
             enddo				  

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo


!	open(32,file='fort.889')
!	write(32,'(6444E25.11)') r(1:iplas,1:nt),
!     & z(1:iplas,1:nt),psi(1:iplas,1:nt),
!     & x(1:65),y(1:65),utemp(1:65,1:65)	
!	close(32)!


          return
          end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_b_ef_old(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!          include'comrec.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                  ui(nip,njp),g(nip,njp),
     +                  ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

           real*8 ron(ntp),ronm(ntp),xstepp,rggr(iplas)
           real*8 psitemp(iplas),uuuu(9),xxx(9,2),clsq(6)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)
           real*8 utemp(nip,njp),um_t,up_t
           real*8 xm_t,ym_t,up_new
	integer imax_t,jmax_t
	
            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            sqrt(xx)=dsqrt(xx)

                !    	write(6,*) 'grid_b:eter'
                !        write(6,*) 'gridb:igdf,nstep',igdf,nstep
                !        call wrb
                !        pause 'gridb:enter:wrb'

!            if(nstep.eq.0) then

	utemp=u
!	write(*,*) igdf
	
!         if(igdf.eq.0)then

!	        do i=1,iplas
!           psia(i)=(iplas-i)/(iplas-1.d0)
!              enddo

!	        do i=1,iplas1
!	     dpsda(i)=-1.d0
!              enddo

!         elseif(igdf.eq.1 .OR. igdf.eq.2 .OR. igdf.eq.3) then

	        do i=1,iplas
           psia(i)=1.d0-((i-1)/(iplas-1.d0))**2
              enddo

	        do i=1,iplas1
           dpsda(i)=(1-2*i)/(iplas-1.d0)
              enddo

!         endif
!         endif

	imax_t=imax
	jmax_t=jmax
	
!find new max
	do k=-3,3
		i= imax+k
		do l=-3,3
			j= jmax+l
			if (u(i,j).gt.u(imax_t,jmax_t)) then
				imax_t=i
				jmax_t=j	
			endif
	enddo
	enddo
!find new mag axis
c---definition of Um - poloidal flux function value at magn. axis.

          nsh=1
          xs(nsh)=x(imax_t)
          ys(nsh)=y(jmax_t)
         fun(nsh)=u(imax_t,jmax_t)

        do 200 k=-1,1

          i= imax_t+k

        do 210 l=-1,1

          j= jmax_t+l

          if(i.eq.imax_t .AND. j.eq.jmax_t) go to 210
          nsh=nsh+1
          xs(nsh)=x(i)
          ys(nsh)=y(j)
         fun(nsh)=u(i,j)

 210    continue
 200    continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call lsq_sur6(xs,ys,fun,nsh,clsq,rrmm,zzmm,uumm,dp)
          xm_t=rrmm
          ym_t=zzmm
          um_t=uumm

	call right0_ef(up_new)

	utemp=u
	up_t=up_new


!	write(9754,*) xm,ym,um,xm_t,ym_t,um_t
!	write(9754,*) imax,jmax,imax_t,jmax_t
!	write(9754,*) up,up_new,up_t

              rm=xm_t    !
              zm=ym_t    !  position of magn. axis
                       !




              psim=um_t-up_t
              psip=0.d0

                dtet=2.d0*pi/(nt-2)
                teta(1)=-dtet

            do 30 j=2,nt
             teta(j)=teta(j-1)+dtet
 30         continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

            psia(1)=1.d0


! first find boundary
	i=iplas
	xstepp=0.25*sqrt((x(2)-x(1))**2.+(y(2)-y(1))**2.)
	do j=1,nt
		do k=1,nip*njp
						ronzz=xstepp*k
             dumx=rm+ronzz*cos(teta(j))
             dumy=zm+ronzz*sin(teta(j))
			j1=((dumx-x(1))/(x(2)-x(1)))+1	
			j2=((dumy-y(1))/(y(2)-y(1)))+1	
			x00=x(j1)
			y00=y(j2)
			x002=x(j1+1)
			y002=y(j2+1)
			u00=blin_(dumx,dumy,x00,x002,y00,y002,utemp(j1,j2)
     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1))
		if (u00.eq.up_t) then
			ro(i,j)=ronzz
			goto 5561
		endif
		if (u00.lt.up_t) then
			ronzz2=xstepp*(k-1)
 	   dumx=rm+ronzz2*cos(teta(j))
 	   dumy=zm+ronzz2*sin(teta(j))
			j1=((dumx-x(1))/(x(2)-x(1)))+1	
			j2=((dumy-y(1))/(y(2)-y(1)))+1	
			x00=x(j1)
			y00=y(j2)
			x002=x(j1+1)
			y002=y(j2+1)
			u002=blin_(dumx,dumy,x00,x002,y00,y002,utemp(j1,j2)
     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1))
			ro(i,j)=ronzz2+(up_t-u002)/(u002-u00)*(ronzz2-ronzz)
			goto 5561
		endif
		enddo

5561	continue

	psitemp(1)=1.
	psitemp(iplas)=0.
	rggr(iplas)=ro(iplas,j)
	rggr(1)=0.
!	write(332,'(6444E25.11)') rm,zm,um_t,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)
	do i=2,iplas-1
		rggr(i)=ro(iplas,j)*(i-1.)/(iplas-1.)
		dumx=rm+rggr(i)*cos(teta(j))
		dumy=zm+rggr(i)*sin(teta(j))
		j1=((dumx-x(1))/(x(2)-x(1)))+1	
		j2=((dumy-y(1))/(y(2)-y(1)))+1	

		xxx(1,1)=x(j1-1)
		xxx(1,2)=y(j2-1)
		uuuu(1)=utemp(j1-1,j2-1)
		xxx(2,1)=x(j1)
		xxx(2,2)=y(j2-1)
		uuuu(2)=utemp(j1,j2-1)
		xxx(3,1)=x(j1-1)
		xxx(3,2)=y(j2)
		uuuu(3)=utemp(j1-1,j2)
		xxx(4,1)=x(j1)
		xxx(4,2)=y(j2)
		uuuu(4)=utemp(j1,j2)
		xxx(5,1)=x(j1-1)
		xxx(5,2)=y(j2+1)
		uuuu(5)=utemp(j1-1,j2+1)
		xxx(6,1)=x(j1)
		xxx(6,2)=y(j2+1)
		uuuu(6)=utemp(j1,j2+1)
		xxx(7,1)=x(j1+1)
		xxx(7,2)=y(j2+1)
		uuuu(7)=utemp(j1+1,j2+1)
		xxx(8,1)=x(j1+1)
		xxx(8,2)=y(j2)
		uuuu(8)=utemp(j1+1,j2)
		xxx(9,1)=x(j1+1)
		xxx(9,2)=y(j2-1)
		uuuu(9)=utemp(j1+1,j2-1)

	call biquadratic_interp_spidef(x(j1-1),x(j1),x(j1+1),
     & y(j2-1),y(j2),y(j2+1),uuuu,
     &  dumx,dumy,
     &  u002)
!	write(*,*) i,j,u002
!		u002=bqlin_(dumx,dumy,xxx,uuuu)
!		u002=blin_(dumx,dumy,xxx(4,1),xxx(7,1),
!     & xxx(4,2),xxx(7,2),uuuu(4),uuuu(8),
!     & uuuu(7),uuuu(6))
			psitemp(i)=(u002-up_t)/(um_t-up_t)
!	write(332,'(6444E25.11)') dumx,dumy,u002,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)
	enddo
!	write(332,'(6444E25.11)') rm+rggr(iplas)*cos(teta(j)),
!     & zm+rggr(iplas)*sin(teta(j)),up_t,
!     & x00,y00,x002,y002,utemp(j1,j2)
!     & ,utemp(j1+1,j2),utemp(j1+1,j2+1),utemp(j1,j2+1),psitemp(i)


!	write(986,*) j
		call reinterp_back_quad((1.-psitemp(1:iplas))**0.5,rggr(1:iplas),
     & iplas, 
     &   (1.-psia(1:iplas))**0.5,ro(1:iplas,j),iplas)

!	write(7489,*) j,(1.-psitemp(1:iplas))**0.5,rggr(1:iplas),
!     & (1.-psia(1:iplas))**0.5,ro(1:iplas,j)

!	if (j.eq.1) then
!			write(*,*) 'j',j,psitemp,rggr,psia(1:iplas),ro(1:iplas,j)
!	endif
	
	enddo !j cycle



            do i=2,iplas

             u0=psia(i)
           ur0=up_t+u0*(um_t-up_t)

             do j=1,nt

             r(i,j)=rm+ro(i,j)*cos(teta(j))
             z(i,j)=zm+ro(i,j)*sin(teta(j))

             psin(i,j)=u0

             enddo  ! j-loop

             enddo   !i-loop

            do  j=1,nt

             r(1,j)=rm
             z(1,j)=zm
             ro(1,j)=0.d0
             psin(1,j)=1.d0

						enddo

        do  i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)
             ro(i,1)=ro(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)
             ro(i,nt)=ro(i,2)

				enddo

             do j=1,nt
               ronor(iplas,j)=1.d0
             enddo				  

             do j=1,nt 
             do i=1,iplas

               ronor(i,j)=ro(i,j)/ro(iplas,j)

               psin(i,j)=psia(i)

             enddo				  
             enddo				  

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

         do i=1,nr
         do j=1,nt

          psi(i,j)=psip+psin(i,j)*(psim-psip)

         enddo
         enddo


!	open(32,file='fort.889')
!	write(32,'(6444E25.11)') r(1:iplas,1:nt),
!     & z(1:iplas,1:nt),psi(1:iplas,1:nt),
!     & x(1:65),y(1:65),utemp(1:65,1:65)	
!	close(32)!


          return
          end

Cquadratic inerpolation

	subroutine biquadratic_interp_spidef(x1,x2,x3,
     & z1,z2,z3,
     & uuuu,r2,thet2,
     &  y2)
!x(j1-1),x(j1),x(j1+1),
!     & y(j2-1),y(j2),y(j2+1),uuuu,
!     &  dumx,dumy,
!     &  u002
	implicit none
	integer i,j,k
	double precision x1,x2,x3,z1,z2,z3
	double precision r2,y2,A,B,C,yy(3),uuuu(9)
	double precision r1(3,3),thet1(3),y1(3,3)
	double precision r0(3),p0(3),t1,thet2,t2,t3

	r1(1,1)=x1
	r1(2,1)=x2
	r1(3,1)=x3
	r1(:,2)=r1(:,1)
	r1(:,3)=r1(:,1)	
	thet1(1)=z1
	thet1(2)=z2
	thet1(3)=z3
	y1(1,1)=uuuu(1)
	y1(2,1)=uuuu(2)
	y1(3,1)=uuuu(9)
	y1(1,2)=uuuu(3)
	y1(2,2)=uuuu(4)
	y1(3,2)=uuuu(8)
	y1(1,3)=uuuu(5)
	y1(2,3)=uuuu(6)
	y1(3,3)=uuuu(7)


	do j=1,3
		r0=r1(:,j)
		p0=y1(:,j)
		t1=p0(1)
		t2=p0(2)
		t3=p0(3)
		z1=r0(1)
		z2=r0(2)
		z3=r0(3)
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2))
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2)
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		yy(j)=A*(r2**2.0)+B*r2+C
	enddo

		t1=yy(1)
		t2=yy(2)
		t3=yy(3)
		z1=thet1(1)
		z2=thet1(2)
		z3=thet1(3)
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2))
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2)
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		y2=A*(thet2**2.0)+B*thet2+C

!	write(*,*) x1,x2,x3,z1,z2,z3,r1,thet1,y1,y2
	end

C----------------------------------------------------------------------|
		subroutine griddata_spidef(x1,y1,Nx1, 
     &   x2,y2,Nx2)
Cend call ! works for random data x1,y1
C
C----------------------------------------------------------------------|
	implicit none

	integer Nx1,Nx2,i,j,k,jdone

	double precision x1(Nx1)
	double precision y1(Nx1)
	double precision x2(Nx2)
	double precision y2(Nx2),A,B,C
	double precision z1,z2,z3,z4,z5,z11
	double precision t1,t2,t3,t4,t5,t11,t10
	double precision z6,z7,z8,z9,t6,t7,t8,t9


	do i=1,Nx2
		jdone=0
		
		t4=x2(i)

		do j=2,Nx1-1
		z1=x1(j-1)
		t1=y1(j-1)
		z2=x1(j)
		t2=y1(j)
		z3=x1(j+1)
		t3=y1(j+1)


		
	if (t4.eq.z1 .and. jdone.eq.0) then
		y2(i)=t1
		jdone=1
!	write(986,*) i,j,z1,z2,z3,t1,t2,t3,t4,y2(i)
	endif

	if (t4.eq.z2 .and. jdone.eq.0) then
		y2(i)=t2
		jdone=1
!	write(986,*) i,j,z1,z2,z3,t1,t2,t3,t4,y2(i)
	endif

	if (t4.eq.z3 .and. jdone.eq.0) then
		y2(i)=t3
		jdone=1
!	write(986,*) i,j,z1,z2,z3,t1,t2,t3,t4,y2(i)
	endif
	
	if (t4.gt.z1 .and. t4.lt.z3 .and. jdone.eq.0) then
		t1=y1(j-1)
		t2=y1(j)
		t3=y1(j+1)
		
				A=(t3-t2-(z3-z2)*(t1-t2)/(z1-z2))
     &     /((z3-z2)*(z3-z1))
				B=(t1-t2)/(z1-z2)
     &     -A*(z1+z2)
				C=t2-A*(z2**2.0)-B*z2

		y2(i)=A*(t4**2.0)+B*t4+C


		jdone=1
!	write(986,*) i,j,z1,z2,z3,t1,t2,t3,t4,y2(i)
	endif

		enddo

	enddo
C	write(55,*) Nx1,x1,y1
C	write(55,*) Nx2,x2,y2
C	pause
	end 
C end of qinterp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grid_b_ef2(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!          include'comrec.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                  ui(nip,njp),g(nip,njp),
     +                  ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

           real*8 ron(ntp),ronm(ntp),xstepp,rggr(iplas)
           real*8 psitemp(iplas)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

 
!	write(*,*) iplas,nt
!	open(32,file='fort.1888')
!	write(32,'(8444E25.11)') r(1:iplas,1:nt),
!     & z(1:iplas,1:nt),psi(1:iplas,1:nt),
!     & x(1:65),y(1:65),u(1:65,1:65),
!     & psia(1:iplas)	
!	close(32)


          return
          end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine grid_b1_ef(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!           include 'comrec.inc'

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                   ui(nip,njp),g(nip,njp),
     +                   ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

            real*8 ron(ntp),ronm(ntp)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
!            sqrt(xx)=dsqrt(xx)

                !    	write(6,*) 'grid_b:enter'
                !        write(6,*) 'gridb:igdf,nstep',igdf,nstep
                !        call wrb
                !        pause 'gridb:enter:wrb'

              rm=xm    !
              zm=ym    !  position of magn. axis
                       !
              !psim=um-up
              psip=0.d0

                dtet=2.d0*pi/(nt-2)
                teta(1)=-dtet

            do 30 j=2,nt
             teta(j)=teta(j-1)+dtet
 30         continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

            psia(1)=1.d0
             u0=0.d0

            call loopL_b(teta,nt,ron,u0)

C          write(6,*) 'loopL i=',i,u0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             do j=1,nt
               ronor(iplas,j)=1.d0
             enddo				  

             do j=1,nt 
             do i=1,iplas

               ro(i,j)=ronor(i,j)*ron(j)

               r(i,j)=rm+ro(i,j)*cos(teta(j))
               z(i,j)=zm+ro(i,j)*sin(teta(j))

               psin(i,j)=psia(i)

             enddo				  
             enddo				  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do 120 j=1,nt

             r(1,j)=rm
             z(1,j)=zm

             ro(1,j)=0.d0
             ronor(1,j)=0.d0
             psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)

             ro(i,1)=ro(i,nt1)
             ronor(i,1)=ronor(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)

             ro(i,nt)=ro(i,2)
             ronor(i,nt)=ronor(i,2)

 333    continue

         do i=1,nr
         do j=1,nt
          psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
         enddo

	!do i=1,iplas
	!dpdpsi(i)=tabp(psia(i))
	!dfdpsi(i)=tabf(psia(i))
	!enddo

           !   call wrb
           !     pause 'grid_b1:wrb'
           !     stop
           !    	write(6,*) 'grid_b:exit'

          return
          end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine grid_b1(igdf,nstep)

           include 'double.inc'
           include 'dim.inc'
           include 'parrc1.inc'
          parameter(nshp=10)
           include 'compol.inc'
!           include 'comrec.inc'

         common /comrz/ x(nip),y(njp),dx(nip),dy(njp),
     +                  dxi(nip),dyj(njp),x12(nip)

        common /compot/ u(nip,njp),ue(nip,njp),un(nip,njp),
     +                   ui(nip,njp),g(nip,njp),
     +                   ux0,ux1,ux2,up,um,xm,ym,
     +                  xx0,yx0,xx1,yx1,xx2,yx2,imax,jmax,
     +                  ix1,jx1,ix2,jx2,
     +                  xx10,yx10,xx20,yx20,xm0,ym0,
     +                  psi_bon

            real*8 ron(ntp),ronm(ntp)
           real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
!            sqrt(xx)=dsqrt(xx)

                !    	write(6,*) 'grid_b:enter'
                !        write(6,*) 'gridb:igdf,nstep',igdf,nstep
                !        call wrb
                !        pause 'gridb:enter:wrb'

              rm=xm    !
              zm=ym    !  position of magn. axis
                       !
              !psim=um-up
              psip=0.d0

                dtet=2.d0*pi/(nt-2)
                teta(1)=-dtet

            do 30 j=2,nt
             teta(j)=teta(j-1)+dtet
 30         continue

            teta(1)=teta(nt1)-2.d0*pi
            teta(nt)=teta(2)+2.d0*pi

            psia(1)=1.d0
             u0=0.d0

            call loopL_b(teta,nt,ron,u0)

C          write(6,*) 'loopL i=',i,u0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             do j=1,nt
               ronor(iplas,j)=1.d0
             enddo				  

             do j=1,nt 
             do i=1,iplas

               ro(i,j)=ronor(i,j)*ron(j)

               r(i,j)=rm+ro(i,j)*cos(teta(j))
               z(i,j)=zm+ro(i,j)*sin(teta(j))

               psin(i,j)=psia(i)

             enddo				  
             enddo				  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do 120 j=1,nt

             r(1,j)=rm
             z(1,j)=zm

             ro(1,j)=0.d0
             ronor(1,j)=0.d0
             psin(1,j)=1.d0

 120        continue

        do 333 i=1,nr

             r(i,1)=r(i,nt1)
             z(i,1)=z(i,nt1)

             ro(i,1)=ro(i,nt1)
             ronor(i,1)=ronor(i,nt1)

             r(i,nt)=r(i,2)
             z(i,nt)=z(i,2)

             ro(i,nt)=ro(i,2)
             ronor(i,nt)=ronor(i,2)

 333    continue

         do i=1,nr
         do j=1,nt
          psi(i,j)=psip+psin(i,j)*(psim-psip)
         enddo
         enddo

	!do i=1,iplas
	!dpdpsi(i)=tabp(psia(i))
	!dfdpsi(i)=tabf(psia(i))
	!enddo

           !   call wrb
           !     pause 'grid_b1:wrb'
           !     stop
           !    	write(6,*) 'grid_b:exit'
!	open(32,file='fort.888')
!	write(32,'(6444E25.11)') r(1:iplas,1:nt),
!     & z(1:iplas,1:nt),psi(1:iplas,1:nt),
!     & x(1:65),y(1:65),u(1:65,1:65)	
!	close(32)

          return
          end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !subroutine arc_x_bnd(nteta,nbtab,rbtab,zbtab)
         subroutine arc_x_bnd(nteta)

         use bnd_modul       

          include   'double.inc'
          include 'dim.inc'
          !parameter(nbtabp=1000)
          !common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
          !dimension rbtab(*),zbtab(*)


          parameter(ntz=1000, nsz=6*ntz+16)
          real*8 uk1(ntz), vk1(ntz) 
          real*8 ukw(ntz), vkw(ntz) 
          real*8 t1(ntz), t1n(ntz) 
	    real*8 w1(nsz), w2(nsz), w3(nsz)
cx

	common
     *  /c_kpr/kpr

           do ib=2,nbtab
             dr=dabs(rbtab(ib)-rbtab(1))
             dz=dabs(zbtab(ib)-zbtab(1))
            if(dr.lt.1.d-6 .AND. dz.lt.1.d-6) then
             nbnd=ib-1
		   go to 1947
            endif
	     enddo
             nbnd=nbtab

 1947   continue

cx
        m=nbnd+1
        do j=1,nbnd
          ukw(j)=rbtab(j)
          vkw(j)=zbtab(j)
        enddo
        ukw(m)=ukw(1)
        vkw(m)=vkw(1)
cx        min(z) as a starting
c         zmin=vkw(1)
c         jvmin=1
c         do j=2,m
c           if(vkw(j).lt.zmin) then
c             zmin=vkw(j)
c             jvmin=j
c           endif
c         enddo
cx      x-point as starting
cx        internal point as middle
        rmin=ukw(1)
        rmax=ukw(1)
        zmin=vkw(1)
        zmax=vkw(1)
        do j=2,m
          if(ukw(j).lt.rmin) then
            rmin=ukw(j)
          endif
          if(ukw(j).gt.rmax) then
            rmax=ukw(j)
          endif
          if(vkw(j).lt.zmin) then
            zmin=vkw(j)
          endif
          if(vkw(j).gt.zmax) then
            zmax=vkw(j)
          endif
        enddo
        rint=0.5d0*(rmin+rmax)
        zint=0.5d0*(zmin+zmax)
cx        orientation
        clock=(ukw(1)-rint)*(vkw(2)-zint)-(vkw(1)-zint)*(ukw(2)-rint)
        j=1
        cmin=((ukw(j)-ukw(m-1))*(ukw(j+1)-ukw(j))
     &       +(vkw(j)-vkw(m-1))*(vkw(j+1)-vkw(j)))
     &      /dsqrt(((ukw(j)-ukw(m-1))**2+(vkw(j)-vkw(m-1))**2)
     &            *((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
        jvmin=1
        do j=2,m-1
          ccur=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))
     &         +(vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j)))
     &        /dsqrt(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2)
     &              *((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
!!aai!!           if(clock.lt.0.d0) ccur=-ccur
          if(ccur.lt.cmin) then
            cmin=ccur
            jvmin=j
          endif
        enddo
cx
        do j=jvmin,m-1
          uk1(j-jvmin+1)=ukw(j)
          vk1(j-jvmin+1)=vkw(j)
        enddo
        do j=1,jvmin-1
          uk1(m-jvmin+j)=ukw(j)
          vk1(m-jvmin+j)=vkw(j)
        enddo
        uk1(m)=uk1(1)
        vk1(m)=vk1(1)

          if(kpr.eq.1) then
        write(*,*) ' bound points ',m
        write(*,*) ' bound start point ',uk1(1),vk1(1)
        write(*,*) ' min angle ',
     &(1.d0-dacos(cmin)/(4.d0*datan(1.d0)))*180.d0,' degree'
          endif

        if((1.d0-dacos(cmin)/(4.d0*datan(1.d0)))*180.d0.lt.100.d0) then
          uk1(2)=0.5d0*(uk1(1)+uk1(3))
          uk1(m-1)=0.5d0*(uk1(m)+uk1(m-2))
          vk1(2)=0.5d0*(vk1(1)+vk1(3))
          vk1(m-1)=0.5d0*(vk1(m)+vk1(m-2))
          if(kpr.eq.1) then
          write(*,*) ' linear interpolated near x-point '
	    endif
	  endif
cx
        do j=1,m
	    ukw(j)=uk1(j)
	    vkw(j)=vk1(j)
	  enddo
cx      second x-point search
        j=2
        cmin=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))
     &       +(vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j)))
     &      /dsqrt(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2)
     &            *((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
        jvmin2=2
        do j=2,m-1
          ccur=((ukw(j)-ukw(j-1))*(ukw(j+1)-ukw(j))
     &         +(vkw(j)-vkw(j-1))*(vkw(j+1)-vkw(j)))
     &        /dsqrt(((ukw(j)-ukw(j-1))**2+(vkw(j)-vkw(j-1))**2)
     &              *((ukw(j+1)-ukw(j))**2+(vkw(j+1)-vkw(j))**2))
!!aai!!   if(clock.lt.0.d0) ccur=-ccur
          if(ccur.lt.cmin) then
            cmin=ccur
            jvmin2=j
          endif
        enddo
cx      check the x-point angle
      if((1.d0-dacos(cmin)/(4.d0*datan(1.d0)))*180.d0.lt.100.d0) then !!!!! 
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!

          uk1(jvmin2+1)=0.5d0*(uk1(jvmin2)+uk1(jvmin2+2))
          uk1(jvmin2-1)=0.5d0*(uk1(jvmin2)+uk1(jvmin2-2))
          vk1(jvmin2+1)=0.5d0*(vk1(jvmin2)+vk1(jvmin2+2))
          vk1(jvmin2-1)=0.5d0*(vk1(jvmin2)+vk1(jvmin2-2))

          if(kpr.eq.1) then
          write(*,*) ' min angle ',
     &     (1.d0-dacos(cmin)/(4.d0*datan(1.d0)))*180.d0,' degree'
          write(*,*) ' linear interpolated near the second x-point '
          endif

cx        arclength from the first to the second x-point
          x2len=0.d0
          tolen=0.d0
	    do j=2,m
	      culen=dsqrt((uk1(j)-uk1(j-1))**2+(vk1(j)-vk1(j-1))**2)
	      tolen=tolen+culen
	      if(j.le.jvmin2) then
	        x2len=x2len+culen
	      endif
	    enddo
cx        distribute the point over the two branches
cx        first branch
        m11=nteta-1
        m1x=(m11-1)*(x2len/tolen)+1
	  m2x=m11-m1x+1
!        write(*,*) ' arclength step ratio on sx branches ',
!     &             (x2len/(m1x-1)) / ((tolen-x2len)/(m2x-1))
        mx=jvmin2
        hbnd=1.d0/(m1x-1)
        t1(1)=0.d0
        do j=2,m1x-1
          t1(j)=t1(j-1)+hbnd
        enddo
        t1(m1x)=1.d0
        CALL ARCM(MX,UK1,VK1,T1N)
        CALL SPLNA1(NTZ,MX,T1N,UK1,M1X,T1,NSZ,W1,W2,W3)
        CALL SPLNA1(NTZ,MX,T1N,VK1,M1X,T1,NSZ,W1,W2,W3)
        
cx
        ISBO=20
        EPSBO=1.d-10
        DO ISB=1,ISBO
          ISBW=ISB
          CALL ARCM(M1X ,UK1,VK1,T1N)
          CALL SPLNA1(NTZ,M1X ,T1N,UK1,M1X,T1,NSZ,W1,W2,W3)
          CALL SPLNA1(NTZ,M1X ,T1N,VK1,M1X,T1,NSZ,W1,W2,W3)
          ABO=0.d0
          DO J=1,M1X
            ABO=DMAX1(ABO,DABS(T1N(J)-T1(J)))
          ENDDO
          IF(ABO.LE.EPSBO) GO TO 7776
        ENDDO
 7776   CONTINUE
!        WRITE(*,*) ' '
!        WRITE(*,*) '     ITERATIONS TO MATCH ARCL. DISTR. ',
!     &                   'AT THE BOUNDARY ',ISBW
cx        second branch
        mx=m-jvmin2+1
        do j=1,mx
	    uk1(j+m1x-1)=ukw(j+jvmin2-1)
	    vk1(j+m1x-1)=vkw(j+jvmin2-1)
	  enddo
        hbnd=1.d0/(m2x-1)
        t1(1)=0.d0
        do j=2,m2x-1
          t1(j)=t1(j-1)+hbnd
        enddo
        t1(m2x)=1.d0
        CALL ARCM(MX,UK1(m1x),VK1(m1x),T1N)
        CALL SPLNA1(NTZ,MX,T1N,UK1(m1x),M2X,T1,NSZ,W1,W2,W3)
        CALL SPLNA1(NTZ,MX,T1N,VK1(m1x),M2X,T1,NSZ,W1,W2,W3)
        
cx
        ISBO=20
        EPSBO=1.d-10
        DO ISB=1,ISBO
          ISBW=ISB
          CALL ARCM(M2X ,UK1(m1x),VK1(m1x),T1N)
          CALL SPLNA1(NTZ,M2X ,T1N,UK1(m1x),M2X,T1,NSZ,W1,W2,W3)
          CALL SPLNA1(NTZ,M2X ,T1N,VK1(m1x),M2X,T1,NSZ,W1,W2,W3)
          ABO=0.d0
          DO J=1,M2X
            ABO=DMAX1(ABO,DABS(T1N(J)-T1(J)))
          ENDDO
          IF(ABO.LE.EPSBO) GO TO 7777
        ENDDO
 7777   CONTINUE
        WRITE(*,*) ' '
        WRITE(*,*) ' 2x    ITERATIONS TO MATCH ARCL. DISTR. ',
     &                   'AT THE BOUNDARY ',ISBW
	  else
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
        !WRITE(*,*) 'no second x-point '
cx        no second x-point
cx      respline to equal arclength
        m11=nteta-1
        hbnd=1.d0/(m11-1)
        t1(1)=0.d0
        do j=2,m11-1
          t1(j)=t1(j-1)+hbnd
        enddo
        t1(m11)=1.d0
        CALL ARCM(M ,UK1,VK1,T1N)
        CALL SPLNA1(NTZ,M ,T1N,UK1,M11,T1,NSZ,W1,W2,W3)
        !WRITE(*,*) 'SPLNA1 for uk1'
        CALL SPLNA1(NTZ,M ,T1N,VK1,M11,T1,NSZ,W1,W2,W3)
        !WRITE(*,*) 'SPLNA1 for wk1'
        
cx
C...    ARCLENGTH MESH FOR NEW BO. (ISBO ITERATIONS)
        ISBO=20
        EPSBO=1.d-10
        DO ISB=1,ISBO
          ISBW=ISB
          CALL ARCM(M11 ,UK1,VK1,T1N)
          CALL SPLNA1(NTZ,M11 ,T1N,UK1,M11,T1,NSZ,W1,W2,W3)
          CALL SPLNA1(NTZ,M11 ,T1N,VK1,M11,T1,NSZ,W1,W2,W3)
          ABO=0.d0
          DO J=1,M11
            ABO=DMAX1(ABO,DABS(T1N(J)-T1(J)))
          ENDDO
          IF(ABO.LE.EPSBO) GO TO 7778
        ENDDO
 7778   CONTINUE
        !WRITE(*,*) ' '
        !WRITE(*,*) ' 1x   ITERATIONS TO MATCH ARCL. DISTR. ',
!     &                   'AT THE BOUNDARY ',ISBW
	  endif
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!

        nbtab=m11
      deallocate( rbtab, zbtab )
      allocate( rbtab(nbtab), zbtab(nbtab) )

cx
        do j=1,m11
c aai
         if(clock.gt.0.) then 
          rbtab(j)=uk1(j)
          zbtab(j)=vk1(j)
         else
          rbtab(j)=uk1(m11-j+1)
          zbtab(j)=vk1(m11-j+1)
         endif
c aai
        enddo
        !nbtab=m11

       return
	 end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
          SUBROUTINE ARCM(M,US,VS,AW)
C
C...      ARCLENGTH MESH FOR GIVEN CURVE
C
          implicit real*8(a-h,o-z)
C
          DIMENSION US(1:M),VS(1:M),AW(1:M)
C...      A.L. COMPUTATION
          AW(1)=0.d0
          DO J=2,M
            AWJ=DSQRT( (US(J)-US(J-1))**2+(VS(J)-VS(J-1))**2 )
            AW(J)=AW(J-1)+AWJ
          ENDDO
C
          RWM=1.d0/AW(M)
          DO J=1,M
            AW(J)=AW(J)*RWM
          ENDDO
C
          RETURN
          END
C       
C...    SPLINE RECOSTRUCTION FOR FU GIVEN NEW MESH: SK(NK)
C
        SUBROUTINE SPLNA1(NAZ,N,S,FU,NK,SK,NSZ,W,WW,WWW)
        implicit real*8(a-h,o-z)
C
        DIMENSION S(1:NAZ)
        DIMENSION FU(1:NAZ)
        DIMENSION SK(1:NAZ)
        DIMENSION W(1:NSZ),WW(1:NSZ),WWW(1:NSZ)
        DIMENSION CW(1:4)
C
C       ... SPLINE INTERPOLATION
C
        N4=N+4
C
C...    MATCH ENDS
        SK(1)=S(1)
        SK(NK)=S(N)
C
        CALL E01BAF(N,S,FU,W,WW,NSZ,WWW,NSZ,IFAIL)
        DO I=1,NK
C         CALL E02BBF(N4,W,WW,SK(I),CW,IFAIL)
          CALL E02BCF(N4,W,WW,SK(I),0,CW,IFAIL)
          FU(I)=CW(1)
        ENDDO
C
        RETURN
        end

         subroutine bnd_mov

           include 'double.inc'
           include 'dim.inc'
          parameter(nbtabp=1000)
          parameter(nb4=nbtabp+4,nb6=nb4*6)
           include 'compol.inc'
           common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
	common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh

           real*8 robn(nbtabp),tetbn(nbtabp)
           real*8 RRK(nb4),CCK(nb4),WRK(nb6),CWK(4)

!            cos(xx)=dcos(xx)
!            acos(xx)=dacos(xx)
!            sin(xx)=dsin(xx)
!            sqrt(xx)=dsqrt(xx)

             do ib=1,nbtab
              rbtab(ib)=rbtab(ib)+0.01
             enddo

        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine furgrid(mdim,mg,af_r,bf_r,af_z,bf_z)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'
            !parameter(mdim=50)           
            dimension af_r(nrp,0:mdim),bf_r(nrp,0:mdim)
            dimension af_z(nrp,0:mdim),bf_z(nrp,0:mdim)

            dimension x(nrp,ntp),y(nrp,ntp)
            
            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
!            sqrt(xx)=dsqrt(xx)


!!! r(i,j)=SUMMA  (af_r(i,m)*cos(m*teta)+bf_r(i,m)*sin(m*teta))
!!! z(i,j)=SUMMA  (af_z(i,m)*cos(m*teta)+bf_z(i,m)*sin(m*teta))
!!!        m=0,mg


             !mg=6

            do i=1,iplas
             do m=0,mg
              af_r(i,m)=0.d0
              bf_r(i,m)=0.d0
              af_z(i,m)=0.d0
              bf_z(i,m)=0.d0
             enddo
            enddo

            do i=1,iplas
            
              do j=2,nt-1
               rc=0.5d0*(r(i,j)+r(i,j+1))              
               zc=0.5d0*(z(i,j)+z(i,j+1))              
               dtet=0.5d0*(teta(j+1)-teta(j))/pi              
               af_r(i,0)=af_r(i,0)+rc*dtet
               af_z(i,0)=af_z(i,0)+zc*dtet
              enddo
           
           if(i.gt.1) then           
             do m=1,mg
              do j=2,nt-1
               tet0=teta(j)*m              
               tet1=teta(j+1)*m              
               r0=r(i,j)
               r1=r(i,j+1) 
               z0=z(i,j)
               z1=z(i,j+1) 
               rc=0.5d0*(r(i,j)+r(i,j+1))   
               zc=0.5d0*(z(i,j)+z(i,j+1)) 
               tetc=0.5d0*(teta(j+1)+teta(j))*m             
               dtet=(teta(j+1)-teta(j))/pi
               
               drc= (r1-r0)/(teta(j+1)-teta(j))
               dzc= (z1-z0)/(teta(j+1)-teta(j))
               sint0=sin(tet0)
               sint1=sin(tet1)
               cost0=cos(tet0)
               cost1=cos(tet1)
                             
               af_r1= rc*(sint1-sint0)/m/pi
               af_z1= zc*(sint1-sint0)/m/pi
               bf_r1=-rc*(cost1-cost0)/m/pi
               bf_z1=-zc*(cost1-cost0)/m/pi
               
           af_r2=drc*( 0.5d0*(sint0+sint1)*dtet+(cost1-cost0)/pi/m )/m
           af_z2=dzc*( 0.5d0*(sint0+sint1)*dtet+(cost1-cost0)/pi/m )/m
           bf_r2=drc*(-0.5d0*(cost0+cost1)*dtet+(sint1-sint0)/pi/m )/m
           bf_z2=dzc*(-0.5d0*(cost0+cost1)*dtet+(sint1-sint0)/pi/m )/m
             
             af_r(i,m)=af_r(i,m)+af_r1+af_r2
             af_z(i,m)=af_z(i,m)+af_z1+af_z2
             bf_r(i,m)=bf_r(i,m)+bf_r1+bf_r2
             bf_z(i,m)=bf_z(i,m)+bf_z1+bf_z2
                            
!               af_r1=af_r(i,m)+rc*dtet*cos(tetc)
!               af_z1=af_z(i,m)+zc*dtet*cos(tetc)
!               bf_r1=bf_r(i,m)+rc*dtet*sin(tetc)
!               bf_z1=bf_z(i,m)+zc*dtet*sin(tetc)
!             af_r2=af_r(i,m)+0.5d0*(r0*cos(tet0)+r1*cos(tet1))*dtet
!             af_z2=af_z(i,m)+0.5d0*(z0*cos(tet0)+z1*cos(tet1))*dtet
!             bf_r2=bf_r(i,m)+0.5d0*(r0*sin(tet0)+r1*sin(tet1))*dtet
!             bf_z2=bf_z(i,m)+0.5d0*(z0*sin(tet0)+z1*sin(tet1))*dtet
!             af_r(i,m)=0.5d0*(af_r1+af_r2)
!             af_z(i,m)=0.5d0*(af_z1+af_z2)
!             bf_r(i,m)=0.5d0*(bf_r1+bf_r2)
!             bf_z(i,m)=0.5d0*(bf_z1+bf_z2)
              enddo
             enddo
           endif          

            enddo

              go to 100

!!test!!!!!!!!!!!!

            do i=1,iplas
             do j=1,nt
              x(i,j)=0.d0
              y(i,j)=0.d0
             enddo
            enddo

            do i=1,iplas
             do j=2,nt-1
              x(i,j)=af_r(i,0)
              y(i,j)=af_z(i,0)
              tetk=teta(j)             
              do m=1,mg
         x(i,j)=x(i,j) + af_r(i,m)*cos(m*tetk)+ bf_r(i,m)*sin(m*tetk)
         y(i,j)=y(i,j) + af_z(i,m)*cos(m*tetk)+ bf_z(i,m)*sin(m*tetk)
              enddo
             enddo
            enddo

            do i=1,iplas
             x(i,1)=x(i,nt-1)
             y(i,1)=y(i,nt-1)
             x(i,nt)=x(i,2)
             y(i,nt)=y(i,2)
            enddo


            do i=1,iplas
             do j=1,nt
              r(i,j)=x(i,j)
              z(i,j)=y(i,j)
             enddo
            enddo

             call wrb
             
!!test!!!!!!!!!!!!

 100       continue


         return
         end        
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine furvsl(mdim,mg,af_r,bf_r,af_z,bf_z)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'
         common /comvsl/ tetvs(nrp,ntp),dtetc(nrp,ntp)
            !parameter(mdim=50)           
            dimension af_r(nrp,0:mdim),bf_r(nrp,0:mdim)
            dimension af_z(nrp,0:mdim),bf_z(nrp,0:mdim)

            dimension x(nrp,ntp),y(nrp,ntp)
            
            cos(xx)=dcos(xx)
            sin(xx)=dsin(xx)
            !sqrt(xx)=dsqrt(xx)


!!! r(i,j)=SUMMA  (af_r(i,m)*cos(m*teta)+bf_r(i,m)*sin(m*teta))
!!! z(i,j)=SUMMA  (af_z(i,m)*cos(m*teta)+bf_z(i,m)*sin(m*teta))
!!!        m=0,mg


             !mg=6

            do i=1,iplas
             do m=0,mg
              af_r(i,m)=0.d0
              bf_r(i,m)=0.d0
              af_z(i,m)=0.d0
              bf_z(i,m)=0.d0
             enddo
            enddo

            do i=1,iplas
            
              do j=2,nt-1
               rc=0.5d0*(r(i,j)+r(i,j+1))              
               zc=0.5d0*(z(i,j)+z(i,j+1))              
               dtet=0.5d0*(tetvs(i,j+1)-tetvs(i,j))/pi              
               af_r(i,0)=af_r(i,0)+rc*dtet
               af_z(i,0)=af_z(i,0)+zc*dtet
              enddo
           
           if(i.gt.1) then           
             do m=1,mg
              do j=2,nt-1
               tet0=tetvs(i,j)*m              
               tet1=tetvs(i,j+1)*m              
               r0=r(i,j)
               r1=r(i,j+1) 
               z0=z(i,j)
               z1=z(i,j+1) 
               rc=0.5d0*(r(i,j)+r(i,j+1))   
               zc=0.5d0*(z(i,j)+z(i,j+1))
               tetc=0.5d0*(tetvs(i,j+1)+tetvs(i,j))*m             
               dtet=(tetvs(i,j+1)-tetvs(i,j))/pi
               drc= (r1-r0)/(tetvs(i,j+1)-tetvs(i,j))
               dzc= (z1-z0)/(tetvs(i,j+1)-tetvs(i,j))
               sint0=sin(tet0)
               sint1=sin(tet1)
               cost0=cos(tet0)
               cost1=cos(tet1)
                             
               af_r1= rc*(sint1-sint0)/m/pi
               af_z1= zc*(sint1-sint0)/m/pi
               bf_r1=-rc*(cost1-cost0)/m/pi
               bf_z1=-zc*(cost1-cost0)/m/pi
               
           af_r2=drc*( 0.5d0*(sint0+sint1)*dtet+(cost1-cost0)/pi/m )/m
           af_z2=dzc*( 0.5d0*(sint0+sint1)*dtet+(cost1-cost0)/pi/m )/m
           bf_r2=drc*(-0.5d0*(cost0+cost1)*dtet+(sint1-sint0)/pi/m )/m
           bf_z2=dzc*(-0.5d0*(cost0+cost1)*dtet+(sint1-sint0)/pi/m )/m
             
             af_r(i,m)=af_r(i,m)+af_r1+af_r2
             af_z(i,m)=af_z(i,m)+af_z1+af_z2
             bf_r(i,m)=bf_r(i,m)+bf_r1+bf_r2
             bf_z(i,m)=bf_z(i,m)+bf_z1+bf_z2
             
               !af_r1=af_r(i,m)+rc*dtet*cos(tetc)
               !af_z1=af_z(i,m)+zc*dtet*cos(tetc)
               !bf_r1=bf_r(i,m)+rc*dtet*sin(tetc)
               !bf_z1=bf_z(i,m)+zc*dtet*sin(tetc)
             !af_r2=af_r(i,m)+0.5d0*(r0*cos(tet0)+r1*cos(tet1))*dtet
             !af_z2=af_z(i,m)+0.5d0*(z0*cos(tet0)+z1*cos(tet1))*dtet
             !bf_r2=bf_r(i,m)+0.5d0*(r0*sin(tet0)+r1*sin(tet1))*dtet
             !bf_z2=bf_z(i,m)+0.5d0*(z0*sin(tet0)+z1*sin(tet1))*dtet
             !af_r(i,m)=0.5d0*(af_r1+af_r2)
             !af_z(i,m)=0.5d0*(af_z1+af_z2)
             !bf_r(i,m)=0.5d0*(bf_r1+bf_r2)
             !bf_z(i,m)=0.5d0*(bf_z1+bf_z2)
             
              enddo
             enddo
           endif          

            enddo

              go to 100

!!test!!!!!!!!!!!!

            do i=1,iplas
             do j=1,nt
              x(i,j)=0.d0
              y(i,j)=0.d0
             enddo
            enddo

            do i=1,iplas
             do j=2,nt-1
              x(i,j)=af_r(i,0)
              y(i,j)=af_z(i,0)
              tetk=tetvs(i,j)             
              do m=1,mg
         x(i,j)=x(i,j) + af_r(i,m)*cos(m*tetk)+ bf_r(i,m)*sin(m*tetk)
         y(i,j)=y(i,j) + af_z(i,m)*cos(m*tetk)+ bf_z(i,m)*sin(m*tetk)
              enddo
             enddo
            enddo

            do i=1,iplas
             x(i,1)=x(i,nt-1)
             y(i,1)=y(i,nt-1)
             x(i,nt)=x(i,2)
             y(i,nt)=y(i,2)
            enddo


            do i=1,iplas
             do j=1,nt
              r(i,j)=x(i,j)
              z(i,j)=y(i,j)
             enddo
            enddo

             call wrb
             
!!test!!!!!!!!!!!!

 100       continue


         return
         end        
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine vsl_teta
           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'
           
         common /com_mf/ bpol(nrp,ntp),btor(nrp,ntp),
     &                   br(nrp,ntp),bz(nrp,ntp),
     &                   rc(nrp,ntp),zc(nrp,ntp)
           
         common /comvsl/ tetvs(nrp,ntp),dtetc(nrp,ntp)
            
            !sqrt(xx)=dsqrt(xx)
            atan(xx)=datan(xx)

!!!!!find theta zero
          jc0=0
         do j=2,nt-1
           dzm=z(iplas,j)-zm
           dzp=z(iplas,j+1)-zm
           drm=r(iplas,j)-rm
           drp=r(iplas,j+1)-rm
          if(dzm.le.0.d0 .AnD. dzp.ge.0.d0 .AnD. drm.gt.0.d0 ) then
           jc0=j
           tetpl=atan(dzp/drp)
           tetmn=atan(dzm/drm)
           exit
          endif       
         enddo 

!!!!!vsl d_theta in (i+0.5,j)

        do i=1,iplas-1
        
           tetmax = 0.d0
         do j=2,nt-1
          dl=0.5d0*(dlt(i,j)+dlt(i+1,j))
          Bp=bpol(i,j)
          Bf=btor(i,j)
          R05=rc(i,j)
          dtetc(i,j) = Bf*dl/Bp/R05
          tetmax = tetmax+Bf*dl/Bp/R05
         enddo 
  
         do j=2,nt-1
          dtetc(i,j) = dtetc(i,j)*2.d0*pi/tetmax
         enddo 
         
        enddo 
           
!!!!!vsl d_theta in (i,j)
           
         do j=2,nt-1
          dtetc(iplas,j) = 1.5d0*dtetc(iplas-1,j)-0.5d0*dtetc(iplas-2,j)
          !dtetc(iplas,j) = dtetc(iplas-1,j)
         enddo 
         
        do i=iplas-1,2,-1
         do j=2,nt-1
          dtetc(i,j) = 0.5d0*(dtetc(i,j)+dtetc(i-1,j))
         enddo 
        enddo 
           
!!!!!vsl theta in (i,j)
           
        do i=2,iplas
          tetvs(i,2)=0.d0
         do j=3,nt-1
          tetvs(i,j)=tetvs(i,j-1)+dtetc(i,j-1)
         enddo 
        enddo 
          
        do i=2,iplas
        
         tetvspl=tetvs(i,jc0+1)
         tetvsmn=tetvs(i,jc0)
         tetvs0=(tetvsmn*tetpl-tetvspl*tetmn)/(tetpl-tetmn)
         
         do j=2,nt-1
          tetvs(i,j)=tetvs(i,j)-tetvs0
         enddo 
         
        enddo 
         
        !do i=2,iplas
         !do j=2,nt-1
          !if(tetvs(i,j).lt.0.d0) tetvs(i,j)=tetvs(i,j)+2.d0*pi
         !enddo 
        !enddo 
          
          
        do i=2,iplas
          tetvs(i,1)=tetvs(i,nt-1)-2.d0*pi !
          tetvs(i,nt)=tetvs(i,2)+2.d0*pi !
        enddo 
          
         do j=1,nt
          tetvs(1,j)=tetvs(2,j)
         enddo 
          
         
          
!            do i=1,iplas
!             do j=1,nt
!              cur(i,j)=tetvs(i,j)
!             enddo
!            enddo

!             call wrb
           
          
          
         return
         end                   
          
