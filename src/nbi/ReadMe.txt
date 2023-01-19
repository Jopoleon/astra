(1) Neutronics
!I prpopose to introduce 4 new variables in ASTRA status.inc a few arrays  which describe the neutron and tritium sources intensities [1/s] due to fast d or t NBI: 
		stnbdp(j)
		sdnbdp1(j)
		sdnbdp2(j)
		sdnbtp(j)


!How to use to calculate particle source density:

! t-nbi - d-plasmas: source [10^19/m3/s] of 14.1MeV n:Sn14=Ndeut(j)*stnbdp(j)
! d-nbi - d-plasmas: source [10^19/m3/s] of 2.45MeV n:Sn2450=(Ndeut(j)-Nibm(j))*sdnbdp1(j)
! d-nbi - d-plasmas: source [10^19/m3/s] of 1 MeV t:St1008=(Ndeut(j)-Nibm(j))*sdnbdp2(j)
! d-nbi - t-plasmas: source [10^19/m3/s] of 14.1MeV n:Sn14=(Ndeut(j)-Nibm(j))*sdnbtp(j)

======================================

! in present version I borrowed the arrays which have different meaning (nbiinj.f)
! thus in the after introduction of special arrays in ASTRA it would be necessary to comment in nbinj.f the following lines
C BPEMEHHO vvvvvvvvvvvvvvvvvvvvvvvvvvvvv		
	do j=1,na1
		NNBM2(j)	=SNIBM2(j)
		NNBM3(j)	=stnbdp(j)
		SNIBM1(j)	=sdnbdp1(j)
		SNIBM2(j)	=sdnbdp2(j)
		SNIBM3(j)	=sdnbtp(j)
	enddo
C BPEMEHHO^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! when you comment the lines then 
	SNIBM1(j)/SNIBM2(j)/SNIBM3(j)
! will be the H/D/T ion source [10^19/m3/s] from H/D/T NBI which in general will be shifted in respect to the electron source SNEBM(j)  

======================================


!there are a few options of the control of the sources foreseen to speed-up simulations in the case when the sources are negligible:

C CNBI4	= control of beam-plasma fusion 	            (Def=1)
C      	!  1/0/2 finite Te,Ti/off/Ti=0,finite Te 

======================================
! all neutron/t sources due to beams are calculated for relaxed fast particle profile (i.e. correspond to CNB4	=1)
! in genereral it could be extended to timedependent disrtibution function? but it could bec0me much more time consuming


===============================================================================================================
(2) Switch to TORIC+SSFPQL:
!
C CNB4	=0/1/2 fast ion source only (for TORIC+SSFPQL)/ steady F-P /time dependent Fokker-Planck Solver	(Def=1)

! If CNB4=0 then source for toric is saved in the file 
	dat/toric.nbi 
! on the transport radial mesh NA1-1 (but 1 point less)
! All other sources caused by presence of the NBI are not calculated (to be calculated by TORIC+SSFPQL)

! the last point jR=NA1 is missed because the source there is 0 by definition. 
! If that is a problem for consistency in TORIC (less points (NA1-1) than in transport grid) 
! I can fill it by "0" in ASTRA output?

================================================================================================================
(3) Output of the distribution function
! In the case of the time dependent FP solver (CNB4	=2) 
! it is possible to save the distribution function in ASCII to the file 
	dat/nbion.dat
! to provide the output you have to put 1 in the file 
	dat/nbi.dat
================================================================================================================
(4) time-dependent FP timestep
! Accuracy of linerisation drops with increase of timestep ~ < TAUs/10. 
! For ITER case the desirable maximal timesteps TAU ~ 0.1 - 0.05 s

