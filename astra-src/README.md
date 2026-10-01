Repository of the ASTRA code

## Maintainer

Current maintainer:
- Giovanni Tardini <giovanni.tardini@ipp.mpg.de>
- Emiliano Fable <emiliano.fable@ipp.mpg.de>

## License

This project is licensed under the GNU Lesser General Public License v2.1
(or later). See the LICENSE file for details.

Clone:
```
  git clone git@gitlab.mpcdf.mpg.de:git/astra.git a8
```

Install:
```
  cd a8
  chmod u+x install.sh
  ./install.sh
```

Compile or execute:
```
  cd a8
  exe/as_exe
```

Supported platforms:
  IPP tok
  IPP hz-ld-prod
  IPP-cz
  gateway
  iter-sdcc
  GA-iris
  GA-omega
  Perlmutter
  mit.edu
  puhti (VTT)
  rat2 (Padua)
  freia (ukaea)

The supported platforms are automatically recognised. Check with
```
  cd a8
  ./get_platform
```
