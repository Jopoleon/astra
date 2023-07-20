subroutine rrsudg(time_ext, dt_smlk)

use fs_coupling_variables, only: fs_dt_smlk, fs_ipl_in
use const_inc, only: BTOR, CPEL1, CIMP3, CBND3, &
        ZRD70, ZRD71, ZRD73, ZRD84, ZRD93, &
        CSCL1, CDYM3, CDVM7, CDWM5, CDWM6, &
        CDMJ1, CDMJ2, CDMJ3, CDMJ4, CDJM5, CDJM6, CDJM7, CDJM8, &
        CSOL1, &
        CHE1, CHE3, &
        CDHJ1, CDHJ2, CDHJ3, CDHJ4, CDHJ5, &
        CV3, CV4, CV6, CV13
use status_inc, only: CAR32, CAR33
use outcmn_inc, only: machine, vcoil

implicit none

real*8, intent(out) :: time_ext, dt_smlk
real*8 :: vcoiltmp(10), gvcoil(15)

save vcoiltmp

if (trim(MACHINE) == 'demo') then
    call shmr( &
        CPEL1, CIMP3, CV4,  CBND3, &
        CSCL1, ZRD71, ZRD73, CDYM3, &
        CDWM5, CDWM6, &
        CDJM5, CDJM6, CDJM7,  CDJM8, ZRD70, &
        CV13, CSOL1, &
        CHE1, CDHJ1, CDHJ2, &
        CHE3, CDHJ3, CDHJ4, &
        CV3, CDHJ5, &
        CDMJ1, CDMJ2, CDMJ3, CDMJ4, vcoil, &
        time_ext, CV6, fs_dt_smlk)
        CV13 = MAX(1., CV13)  ! finite pump speed to avoid NaN
        dt_smlk = fs_dt_smlk  ! simulink tau defined in equ log
elseif (trim(MACHINE) == 'iter') then
    call shmr( &
        CPEL1, CIMP3, CV4,  CBND3, &
        CSCL1, ZRD71, ZRD73, CDYM3, &
        CDWM5, CDWM6, &
        CDJM5, CDJM6, CDJM7,  CDJM8, ZRD70, &
        CV13, CSOL1, &
        CHE1, CDHJ1, CDHJ2, &
        CHE3, CDHJ3, CDHJ4, &
        CV3, CDHJ5, &
        CDMJ1, CDMJ2, CDMJ3, CDMJ4, gvcoil, &
        time_ext, CV6, fs_dt_smlk)
        CV13 = MAX(1., CV13)  ! finite pump speed to avoid NaN
        dt_smlk = fs_dt_smlk  ! simulink tau defined in equ log
else if (trim(MACHINE) == 'aug') then
    call shmr( &
        CPEL1, CV13, CDMJ1, CDMJ2, CDMJ3, CDMJ4, &
        ZRD84, CAR32(1: 8), CAR32(9: 16), CAR32(17: 24), CAR32(25: 26), &
        vcoiltmp(1: 10), CAR33(1: 24), BTOR, &
        time_ext, CV6, fs_dt_smlk)

    BTOR = abs(BTOR) !Btor defined here absolute value. sign has to be given separatly
    dt_smlk = fs_dt_smlk ! simulink tau defined in equ log

    vcoil(1) = vcoiltmp(1) - vcoiltmp(2)
    vcoil(2) = vcoiltmp(2) - vcoiltmp(3)
    vcoil( 3: 10) = vcoiltmp(3: 10)
    vcoil(11: 12) = 0.

    ZRD93 = time_ext + dt_smlk
endif  

return
end subroutine rrsudg
