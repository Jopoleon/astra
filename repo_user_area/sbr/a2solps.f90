subroutine A2SOLPS(irad, i_mod)

! Write the input file for SOLPS to do its own calculation of 1 time step of astra  
! equ-file call:
!     A2SOLPS(nint(CBND2*ROC/HRO), int(CSOL3))>::::;

use scalars, only: NA1, ROC, AIM1, AIM2, AIM3
use status, only: CAR45, CAR46, CAR47, QN, &
    QE, QI, QF1, QF2, QF3, DF1, DF2, DF3, VF1, VF2, VF3, &
    NMAIN, NIZ1, NIZ2, NIZ3, ZMAIN, ZIM1, ZIM2, ZIM3, &
    XI, HE, DN, CN, TE, TI, SLAT, RHO
use read_input, only: awd, equ_file, exp_file

implicit none

integer, intent(in) :: irad, i_mod
integer :: switch, lang, jrho, iostatus, n_imp, i_part, remain, ios
double precision, dimension(NA1) :: qmain 
double precision :: dna1, dna2, dna3, dna4, hce, hci, vsa
character(len=160) :: sicas_dir, date_str, f_switch='astra_switch.txt'

data i_part/0/
save i_part

write(*, *) 'A2SOLPS'
sicas_dir = TRIM(awd) // '/sicas/' // TRIM(exp_file) // '/' // TRIM(equ_file)
call execute_command_line('mkdir -p ' // TRIM(sicas_dir))
call chdir(TRIM(sicas_dir))

if (i_part == 0) then
    i_part = i_part + 1
    return
endif

NIZ1(1:irad) = CAR45(1:irad)
NIZ2(1:irad) = CAR46(1:irad)
NIZ3(1:irad) = CAR47(1:irad)
remain = MOD(i_part, i_mod)
i_part = i_part + 1
if (remain /= 0) return

lang = NA1 - irad + 1

do jrho=1, NA1
    qmain(jrho) = (qn(jrho)  - zim1(jrho)*qf1(jrho) - &
        zim2(jrho)*qf2(jrho) - zim3(jrho)*qf3(jrho))/zmain(jrho)
enddo

n_imp = 0
write(1001, *) ' to solps fluxes in (MW,10^19 part/s)/m^2:'
write(1001, *) lang
write(1001, *) QE(irad)
write(1001, *) QI(irad)
write(1001, *) QN(irad)
write(1001, *) QMAIN(irad)
write(1001, *) TE(irad)*1000
write(1001, *) TI(irad)*1000
write(1001, *) (DN(jrho), jrho=irad, NA1)
write(1001, *) (HE(jrho), jrho=irad, NA1)
write(1001, *) (XI(jrho), jrho=irad, NA1)
write(1001,*) SLAT(irad)
if (AIM1 >= 0.5) then
    write(1001, *) (DF1(jrho), jrho=irad, NA1)
endif
if (AIM2 >= 0.5) then
    write(1001, *) (DF2(jrho), jrho=irad, NA1)
endif
if (AIM3 >= 0.5) then
    write(1001, *) (DF3(jrho), jrho=irad, NA1)
endif
if (AIM1 >= 0.5) then
    write(1001, *) QF1(irad)
    n_imp = n_imp + 1
endif
if (AIM2 >= 0.5) then
    write(1001, *) QF2(irad)
    n_imp = n_imp + 1
endif
if (AIM3 >= 0.5) then
    write(1001, *) QF3(irad)
    n_imp = n_imp + 1
endif
write(1001,*) (RHO(jrho)/ROC, jrho=irad, NA1)
write(1001,*) NMAIN(irad)
if (AIM1 >= 0.5) then
    write(1001, *) NIZ1(irad)
endif
if (AIM2 >= 0.5) then
    write(1001, *) NIZ2(irad)
endif
if (AIM3 >= 0.5) then
    write(1001, *) NIZ3(irad)
endif
if (AIM1 >= 0.5) then
    write(1001, *) nint(ZIM1(irad))
endif
if (AIM2 >= 0.5) then
    write(1001, *) nint(ZIM2(irad))
endif
if (AIM3 >= 0.5) then
    write(1001, *) nint(ZIM3(irad))
endif
write(1001, *) (CN(jrho), jrho=irad, NA1)
if (AIM1 >= 0.5) then
    write(1001, *) (VF1(jrho), jrho=irad, NA1)
endif
if (AIM2 >= 0.5) then
    write(1001, *) (VF2(jrho), jrho=irad, NA1)
endif
if (AIM3 >= 0.5) then
    write(1001, *) (VF3(jrho), jrho=irad, NA1)
endif
write(1001, *) n_imp
close(1001)

write(*, *) 'Before SOLPS switch'

switch = 0
do while (switch == 0)
    open(10001, file=TRIM(f_switch))
    read(10001, *, IOSTAT=iostatus) switch
    close(10001)
    if (switch == 0) then
        open(unit=872, file='solps_switch.txt', status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            print *,  'Error opening file. IOSTAT =',  ios
            stop
        endif
        write(872, '(I0)') 1
        close(872)
    endif
    call SLEEP (30)
enddo
call date_and_time(date=date_str)
write(*, *) 'SOLPS pass'
write(*, *) date_str
open(10001, file=TRIM(f_switch))
write(10001, '(I0)') 0
close(10001)

call chdir(TRIM(awd))

end subroutine A2SOLPS
