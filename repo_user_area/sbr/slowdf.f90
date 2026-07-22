      subroutine SLOWDF(pnbi_powers,slow_time)

!EFable 2025, compute fast ions slowing down for beam ions
			
	use status
        use scalars
        use read_input, only: AWD, nml_file

			implicit none

	integer J,i,ishot
				
! beam energy in keV

	double precision beam_energy, beam_mass,beam_charge
	double precision beam_energy_fractions(3),v0
	
	double precision slow_time(NRD),niffo(nrd),coulg(nrd)
	double precision coulgi(nrd,5),vcrit(nrd),dumb1(nrd),tausso(nrd)
double precision :: part_mix(3, 8),pnbi_powers(nrd), blippo

integer n_nbi, tim_prev
double precision a_beam(8),z_beam(8),einj(8)
data tim_prev/0/	

character(len=120) :: as_nml, pinj_file, pinj_file2, limiter_file, table_path

namelist / partmix / part_mix
namelist / nbi_par / n_nbi, a_beam, z_beam, einj, pinj_file


save tim_prev, n_nbi, einj, part_mix, a_beam, z_beam


if (tim_prev == 0) then  ! init
    tim_prev=1
    as_nml = TRIM(AWD) // TRIM(nml_file)

    open(53, FILE=TRIM(as_nml), delim='apostrophe')
    read(53, nml=partmix)
    read(53, nml=nbi_par)
    close(53)

endif


slow_time=1.
niffo=0.
blippo=0.
 COULGI=0.

 COULG = 15.9 - 0.5*LOG(NE) + log(TE)
!3 ion species: D, T, impurity
 COULGi(1:nrd,1) =  18.48-3.*LOG(1.) -0.5*LOG(NDEUT)+1.5*LOG(TI)
 COULGi(1:nrd,2) =  18.48-3.*LOG(1.) -0.5*LOG(NTRIT)+1.5*LOG(TI)
if (NIZ1(1)>0.) COULGi(1:nrd,3) =  18.48-3.*LOG(ZIM1) -0.5*LOG(NIZ1)+1.5*LOG(TI)
if (NIZ2(1)>0.) COULGi(1:nrd,4) =  18.48-3.*LOG(ZIM2) -0.5*LOG(NIZ2)+1.5*LOG(TI)
if (NIZ3(1)>0.) COULGi(1:nrd,5) =  18.48-3.*LOG(ZIM3) -0.5*LOG(NIZ3)+1.5*LOG(TI)

dumb1=NDEUT*COULGI(1:nrd,1)/2.+NTRIT*COULGI(1:nrd,2)/3.
if (AIM1>0.) dumb1=dumb1+NIZ1*COULGI(1:nrd,3)/AIM1*ZIM1**2.
if (AIM2>0.) dumb1=dumb1+NIZ2*COULGI(1:nrd,4)/AIM2*ZIM2**2.
if (AIM3>0.) dumb1=dumb1+NIZ3*COULGI(1:nrd,5)/AIM3*ZIM3**2.

vcrit=(1.51e14*(TE*1.e3)**(1.5)*dumb1/NE/COULG)**(1./3.)

do j=1,n_nbi

beam_energy_fractions(1)=part_mix(1,j) ! full energy
beam_energy_fractions(2)=part_mix(2,j) ! 2/3
beam_energy_fractions(3)=part_mix(3,j) ! 1/3
beam_mass=a_beam(j)
beam_charge=z_beam(j)
beam_energy=einj(j)/1.e3

tausso=6.27e8*beam_mass/beam_charge**2./COULG*(TE*1.e3)**(1.5)/(NE*1.e13)

do i=1,3
 if (i==1) v0=sqrt(2.*beam_energy*1.e3*1.602e-19/beam_mass/1.67262e-27)
 if (i==2) v0=sqrt(2.*beam_energy*2./3.*1.e3*1.602e-19/beam_mass/1.67262e-27)
 if (i==3) v0=sqrt(2.*beam_energy*1./3.*1.e3*1.602e-19/beam_mass/1.67262e-27)
 niffo=niffo+pnbi_powers(j)*tausso/3.*log((v0**3.+vcrit**3.)/vcrit**3.)*beam_energy_fractions(i)
enddo
 blippo=blippo+pnbi_powers(j)

enddo !launchers

if (blippo>0.) slow_time=blippo/niffo

!write(*,*) 'niffo:',n_nbi, beam_energy,beam_mass,slow_time(1),niffo(1),tausso(1),vcrit(1),v0,dumb1(1)
		
end




