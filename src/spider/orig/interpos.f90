!
MODULE interpos
  !
  ! Module with an interpolation routine using cubic splines with tension
  ! It also allows various boundary conditions:
  !    on the second derivative at end points or anywhere in the input interval
  !    on the first derivative
  !    on the function itself as well (thus can impose the value y(x) at some x position)
  !
  ! The minimization uses the standard chi squarred, included the contribution of the
  ! integral of the 2nd derivative weighted by the value of the tension. Therefore
  ! a tension=0 gives the standard cubic splines. The tension can be given as an array
  ! which will be interpreted as "error bars" for each y_i(x_i), normalised to the smallest error bar.
  !
  ! a test program shows use of cbsplgen0 and splibnd
  !
  ! typical call of cbsplgen0:
  !   call cbsplgen0(xin,yin,nin,xout,yout,youtp,youtpp, &
  !      &  nout,ioptder,iextrapo,tension,nbc,ybc,iflag)
  !       xin, yin (nin): input array
  !       xout,yout,youtp,youtpp(nout): output arrays of values of new spline on xout, 1st der (youtp) and 2nd der (youtpp)
  !       ioptder: 0, 1 or 2 if needs to compute yout only or youtp as well or all 3 (func, 1st der and 2nd der)
  !       iextrapo: option for extrapolation
  !                 0 : stop if out of bound of xin
  !                 1 : linear extrapolation (safest extrapolation)
  !                 2 : quadratic extrapolation (relatively good)
  !                 3 : cubic extrapolation
  !                 21: use quadratic within deltax and linear further
  !                 31: use cubic within deltax and linear    further
  !                 32: use cubic within deltax and quadratic further
  !
  !       tension: value of tension on each xin point
  !               (usually same value on each point, but allows closer match to yin in some regions, like error bars)
  !       nbc(2), ybc(2): left/right boundary conditions. Usually nbc=0, ybc=0 to have 2nd der=0 on both ends
  !           nbc=0, 1 or 2: 2nd der., 1st der. or function is specified by ybc
  !           useful also: nbc(1)=1, ybc(1)=0. to specify 1st der=0 at left boundary
  !           other options available: see comments in cbsplgen0.f90 and cbfitbnd.f90
  !       iflag: if not 0, there was a problem (not well implemented yet)
  !
  ! If one wants to compute first spline and then interpolate many times on various xout meshes:
  !   call cbsplgen0(xin,yin,nin,xin,yinspl,dummy,yinsplpp,nout, &
  !      &  2,iextrapo,tension,nbc,ybc,iflag)
  !   then use xin, yinspl, and yinsplpp anywhere in the program like:
  !   do i=1,nout
  !     call splibnd(xin,yinspl,yinsplpp,nin,xout(i),yout(i), youtp(i), youtpp(i),iextrapo)
  !   end do
  !
  ! This module is provided by CRPP - EPFL. Comments/questions to: olivier.sauter@epfl.ch
  ! It is based on the paper: HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280
  !
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: cbsplgen0                       ! Cubic spline with tension, head call
  PUBLIC :: splibnd                         ! Compute interpolated values given 2nd der of spline

CONTAINS

   include 'cbsplgenpkg.f90'

END MODULE interpos
