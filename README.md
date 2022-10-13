
Repository of the ASTRA code (G. Pereverzev, P. N. Yushmanov)

Supported platforms:

  IPP tok
  IPP lxts
  IPP-cz
  gateway
  iter-sdcc

Server installation: execute
  ./install.sh

Or if you prefer
  module load intel
  make

To have the final installation, edit the 4-liners install.sh once for all and execute it. Adapt then the astra_rc file on the client side to see all the lib paths (astra_core, equil, nbi, rabbit, torbeam, tglf, qlkz).
