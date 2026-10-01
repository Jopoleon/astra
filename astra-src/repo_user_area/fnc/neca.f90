! NECH [10#19/m#3]: Chord average density (r) [m]
!    Integral {0,r} ( NE ) dl / a
!    (Yushmanov 11-MAY-87)
double precision function NECAR()

use scalars, only: NA1
use status, only: NE

implicit none

necar = sum(NE(1:NA1))/NA1

end function NECAR
