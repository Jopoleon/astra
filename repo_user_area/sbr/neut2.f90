	subroutine	NEUT2
!C---------------------------------------------------------April 2021---|
!C	Input:	ABC,NA1,NA,AMJ,NAB,NNCX
!C		AMAIN(j),ZMAIN(j),TE(j),TI(j),NE(j),NI(j),ZEF(j),SNNBM(j)
!C		ENCL,ENWM
!C		NNCL,NNWM
!C	Warning:	none
!C	Output:	NN,	TN,    ALBPL
!C-------------------------Christian Schuster chrischu@ipp.mpg.de-------|
	use parameter_inc
	use const_inc
	use status_inc
	use outcmn_inc
	use iso_fortran_env
!C	implicit none
!	include	'for/parameter.inc'
!	include 'for/const.inc'
!C	double precision ABC,AMJ,ENCL,ENWM,NNCL,NNWM,ALBPL
!C	integer NA1,NA,NAB,NNCX
!	include 'for/status.inc'
!	INCLUDE 'for/outcmn.inc'
	
	
	integer :: iCl, iWm, hiInd, loInd, lastSetNnInd, lastSetIntX
	integer, parameter :: nV = 50
	real(real64) :: vVals(nV), eVals(nV), y(nV), intCoeff(nV), yInit(nV)
	real(real64) :: oldNn(2*NA1), newNn(2*NA1), newTn(2*NA1)
	real(real64) :: sCxNi(NA1), sIonNe(NA1), xVals(NA1)
	real(real64) :: lastNn, lastTn, lastCalculatedNn
!C----------------------------------------------------------------------|
!	
	call calculate()
	
	
	contains
!C======================================================================|

	subroutine calculate()
		implicit none
		integer :: j, ipar(2), idid, i1, i2
		real(real64) :: work(8*nV+22),rpar(2),rtol(nV),atol(nV),x,xe
		integer :: iwork(22)
		
		call init()
		rtol = 1.0D-6
		atol = 1.0D-18
		
		oldNn = 0
		do j=1,3
			
			!print *,sIonNe
			call prepareIter()
			!call rhs(nV, RHO(NA1) * (1-float(2*j) / 1000), y, y, 1., 1)
			work = 0
			iwork = 0
			x = xVals(NA1)
			xe = 0
			!print *,size(work), size(iwork)
			call DOPRI5(nV,rhs,x,y,xe,rtol,atol,1, &    
     &		 solOut,1,work,size(work),iwork,size(iwork),rpar,ipar,idid)
		 	!print *,lastTn
			work = 0
		 	iwork = 0
		 	x = 0
		 	xe = -xVals(NA1)
		 	!print *,'start second integration'
		 	call DOPRI5(nV,rhs,x,y,xe,rtol,atol,1,    &
     &		 solOut,1,work,size(work),iwork,size(iwork),rpar,ipar,idid)
!			newNn = 1e-3
			oldNn = newNn
			!print *,"ran",j,idid
			!print *,newNn
			!print *,sIonNe
!			do while (.true.)
			 
!			end do
			
		enddo
		do j=1,NA1
			i1 = NA1-j+1
			i2 = NA1+j
			nn(j) = newNn(i1) + newNn(i2)
			tn(j) = (newNn(i1)*newTn(i1) + newNn(i2)*newTn(i2)) / nn(j)
			!print *,j,nn(j),tn(j),newTn(i1),i1
		end do
		nn = nn / newNn(1)! nn isn't actually the neutral 
		! density, it is just the relative one
		!print *,nn
		
		
		! last step: calculate the albedo
		!print *,"*****"
		ALBPL = sum(y * vVals * intCoeff)/sum(yInit * vVals * intCoeff)
		!print *,newNn(1),newNn(2*NA1),ALBPL
		!print *,yInit/yInit
		!print *,"*******"
		!print *,xVals
		!do while (.true.)
		
		!end do
	end

	subroutine init()
		implicit none
		integer :: j
		
		do j=1,nV
			eVals(j) = 10 ** ((j-1.) * 4 / (nV-1.0)) * 0.001
		end do
		
!C	ok, we now have to find the closest energy values to ENCL and ENWM
!C 	find the indices for both energies
		iCl = minloc(dabs(ENCL - eVals), dim = 1)
		iWm = minloc(dabs(ENWM - eVals), dim = 1)
!C	Adjust the corresponding energies a little such that 
!C	the users will is fulfilled as accurate as possible.

		if (iCl == iWm) then
			print *,'******************************'
			print *,'ENCL and ENWM are too similar!'
			print *,'undefined behavior in NEUT2!'
			print *,'******************************'
		end if

		eVals(iCl) = ENCL
		eVals(iWm) = ENWM
		
		vVals = dsqrt(2*eVals*1.602e-19*1000/(AMJ*1.6605e-27))
		
!C	calculate the coefficients for quadrature with trapezoidal rule
!C	usage: sum(intCoeff * y) ! gives 1-directional nn
		intCoeff(1) = 0.5*(vVals(2) - vVals(1))
		intCoeff(nV) = 0.5*(vVals(nV) - vVals(nV-1))
		do j=2,nV-1
			intCoeff(j) = 0.5*(vVals(j+1) - vVals(j-1))
		end do
		
		xVals = RHO(1:NA1)
		call setRates()
	end
	
	subroutine setRates()
		integer :: j
		do j=1,NA1
			include 'fml/svcx'
			sCxNi(j)=SVCX*NI(j)
			include 'fml/svie'
			sIonNe(j)=SVIE*NE(j)
		end do
	end
	
	subroutine prepareIter()
		implicit none
		real(real64) :: yTmp(nV)
		integer :: j
		yTmp = 0.
		yTmp(iCl) = 1.
		yTmp = yTmp / sum(yTmp * intCoeff) * NNCL
		
		y = 0.
		y(iWm) = 1.
		y = y / sum(y * intCoeff) * NNWM
		
		y = y + yTmp
		yInit = y
		
		hiInd = NA1
		loInd = NA1-1
		
		lastSetNnInd = 1
		newNn(1) = NNCL + NNWM
		lastNn = newNn(1)
		!print *,'lastNn',lastNn
		lastTn = sum(y * eVals * intCoeff) * 2/lastNn
		newTn(1) = lastTn
		!print *,'lastTn',lastTn
		lastSetIntX = 1
		
		! print *,"int test"
		! print *,sum(exp(-vVals/vVals(1)) * intCoeff)
		! print *,(exp(-vVals(nV/vVals(1)))-exp(-1.))*vVals(1)
!		print *,RHO(1:NA1),AMETR(1:NA1)
	end

	subroutine rhs(n, x, y, f, RPAR, IPAR)
		implicit none
		integer, intent(in) :: n, IPAR
		real(real64), intent(in) :: x, y(n), RPAR
		real(real64), intent(out) :: f(n)
		
		real(real64) :: maxwellian(n)
		real(real64) :: xAbs, tiLoc, sCxNiLoc, sIonNeLoc
		real(real64) :: nLoc, nMirror, d, d1, d2, yInt
		integer      :: jLo, jHi, lLo, lHi
		
		!print *,"rhs",x
		xAbs = abs(x)
		
!C		Determine interval for interpolation
!C		Most likely we are very close to the last position, therefore we
!C		start form there. Otherwise a binary search would be smarter.
		do while(xVals(loInd) > xAbs .and. loInd > 1)
			hiInd = hiInd - 1
			loInd = loInd - 1
		end do
		
		do while(xVals(hiInd) < xAbs .and. hiInd < NA1)
			hiInd = hiInd + 1
			loInd = loInd + 1
		end do
		
!C		Interpolate the values we need. For F(x)=F(-x) quantities:
		d = xVals(hiInd) - xVals(loInd)
		d1 = xVals(hiInd) - xAbs
		d2 = xAbs - xVals(loInd)
		if (xAbs < xVals(1)) then
			tiLoc = ti(1)
			sCxNiLoc = sCxNi(1)
			sIonNeLoc = sIonNe(1)
		else
			tiLoc = (ti(loInd)*d1 + ti(hiInd)*d2) / d
			sCxNiLoc = (sCxNi(loInd)*d1 + sCxNi(hiInd)*d2) / d
			sIonNeLoc = (sIonNe(loInd)*d1 + sIonNe(hiInd)*d2) / d
		end if
!C		For Nn it is a little more complicated
!C		First: outgoing nn from last iteration
		
		yInt = sum(y*intCoeff)
		!print *,'here'
		
		if (x >= 0) then
			jLo = -loInd + NA1 + 1
			jHi = -hiInd + NA1 + 1
			lLo = 2*NA1 + 1 - jLo
			lHi = 2*NA1 + 1 - jHi
			nLoc = yInt
			lastCalculatedNn = yInt
			!print *,"indices",jLo,jHi
			nMirror = (oldNn(lLo)*d1 + oldNn(lHi)*d2) / d
			!print *,nMirror
			!print *,y
			!print *,'****************'
		else
			jLo = loInd + NA1
			jHi = hiInd + NA1
			lLo = 2*NA1 + 1 - jLo
			lHi = 2*NA1 + 1 - jHi
			nLoc = (oldNn(jLo)*d1 + oldNn(jHi)*d2) / d
		
!C			now we have to see what the inward nn is
!C			Close to the axis it can happen that we don't yet have the 
!C			data we need. We then extrapolate with the last value.
			if (xAbs <= xVals(1)) then
				!print *,lastSetNnInd,newNn(lastSetNnInd)
				!nMirror = newNn(lastSetIntX)
				nMirror = lastCalculatedNn
			else
				!print *,x
				
				nMirror = (newNn(lLo)*d1 + newNn(lHi)*d2) / d
				!print *,nMirror
				!if (dabs(nMirror) > 1.0D-8 .and. xAbs < 0.02) then 
					!print *,nMirror
					!print *,x,xVals(loInd),loInd
					!print *,newNn(jLo), newNn(jHi)
					!print *,jLo,jHi,loInd,hiInd,NA1
				!end if
			end if
			
		end if
		
!C		We distribute CX neutrals according to a Maxwellian
!C		(thermalized ions). *2 because eVals only contains positives.
		maxwellian = exp(-eVals / tiLoc)
		maxwellian = maxwellian / (sum(maxwellian * intCoeff)*2)
		
		f = sIonNeLoc*y ! ionization losses
		f = f + sCxNiLoc*nLoc*y/yInt ! CX losses
		f = f - sCxNiLoc*(nLoc+nMirror)*maxwellian ! CX gains
		!if (abs(x) < 1e-2) then
		!print *,x,"alive",f(1),nLoc,nMirror,jLo,jHi!,newNn(lastSetNnInd)
		!print *,lastCalculatedNn,nMirror,nLoc
		!print *,x,sIonNeLoc*y(1),sCxNiLoc*nLoc*y(1)/yInt,
!     &		sCxNiLoc*(nLoc+nMirror)*maxwellian(1)
		!end if
		! if (x > 1.6) then
		! 	print *,"nn",nLoc,nMirror
		! 	print *,x
		! 	print *,lLo,lHi,loInd,hiInd
		! end if
		!print *,"alive",y
		f = f / vVals 
		!print *,nLoc,nMirror,x,maxval(abs(f)),sCxNiLoc,sIonNeLoc
		!print *,sIonNe
	end

	subroutine solOut(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
		implicit none
		integer, intent(in) :: nr, n, nd, ipar, icomp(nd)
		real(real64), intent(in) :: xOld, x, rpar
		real(real64), intent(in) :: y(n),con(5*nd)
		integer :: irtrn
		
		integer :: indToSet, indToSetAbs
		real(real64) :: xSet, nnNow, tnNow, xGrid
	
		indToSet = lastSetIntX + 1
		if (indToSet > 2*NA1) then
			return
		else if (indToSet > NA1) then
			indToSetAbs = indToSet - NA1
			!print *,indToSet,NA1,indToSetAbs
			xGrid = -xVals(indToSetAbs)
		else
			indToSetAbs = NA1 - indToSet + 1
			xGrid = xVals(indToSetAbs)
		end if
		
		!print *,"alive", xGrid, x, indToSet, indToSetAbs
		
		nnNow = sum(y * intCoeff)
		!print *,nnNow
		tnNow = sum(y * eVals * intCoeff) * 2/nnNow
		do while(xGrid >= x .and. indToSetAbs <= NA1)
			!print *,indToSetAbs
			newNn(indToSet) = (lastNn*(x-xGrid) + nnNow*(xGrid-xOld))  & 
     &			/(x-xOld)
			newTn(indToSet) = (lastTn*(x-xGrid) + tnNow*(xGrid-xOld))  &
     &			/(x-xOld)
!	 		print *,"------------"
!			print *,x,indToSet,newNn(indToSet)!,lastNn,nnNow
!	   	 	print *,"------------"
			!print *,lastNn, nnNow, indToSet
			!print *,newNn(indToSet)
			!print *,y
			!print *,'**********'
			!if (newTn(1) .gt. 0.1) then
		!		print *,'*****BEEP*****',indToSet,x
		!	end if
			lastSetIntX = indToSet
			lastSetNnInd = indToSetAbs
			!print *,"set",lastSetNnInd
			!print *,"************",xGrid
			!print *,y
			
			indToSet = indToSet + 1
			if (indToSet > 2*NA1) then
				return
			else if (indToSet > NA1) then
				indToSetAbs = indToSet - NA1
				!print *,indToSet,NA1,indToSetAbs
				xGrid = -xVals(indToSetAbs)
			else
				indToSetAbs = NA1 - indToSet + 1
				xGrid = xVals(indToSetAbs)
			end if
		end do
		lastNn = nnNow
		lastTn = tnNow
	return
	end

	end
