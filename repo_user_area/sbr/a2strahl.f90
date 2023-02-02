subroutine A2STRAHL(tau_start, zneocl, dzneocl, dimpsol, ydimp, yvimp, &
    y_zcharge, ymimp0, ynimp, yzeff, ydneo, yvneo, rrates_in, &
    prate1, prate2, shot_in)

use parameter_inc, only: NRD
use strahl, only: astra2strahl

implicit none

double precision, intent(in) :: tau_start, zneocl, dzneocl, dimpsol, &
    rrates_in, prate1, prate2, shot_in
double precision, dimension(NRD), intent(in) :: ydneo, yvneo
double precision, intent(out) :: ymimp0
double precision, dimension(NRD), intent(out) :: ydimp, yvimp, y_zcharge, ynimp, yzeff

call ASTRA2STRAHL(tau_start, zneocl, dzneocl, dimpsol, ydimp, yvimp, &
    y_zcharge, ymimp0, ynimp, yzeff, ydneo, yvneo, rrates_in, &
    prate1, prate2, shot_in)

return
end subroutine a2strahl
