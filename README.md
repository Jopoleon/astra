Repository of the ASTRA code

## Maintainer

Current maintainers:
- Giovanni Tardini <giovanni.tardini@ipp.mpg.de>
- Emiliano Fable <emiliano.fable@ipp.mpg.de>

## License

This project is licensed under the GNU Lesser General Public License v2.1
(or later). See the LICENSE file for details.

Clone:
```
  git clone https://gitlab.mpcdf.mpg.de/git/astra a8
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
  Eurofusion gateway
  iter-sdcc
  GA-iris
  GA-omega
  Perlmutter
  mit.edu
  puhti, tohtori (VTT)
  rat2 (Padua)
  freia (ukaea)
  Columbia university
  Sevilla university

The supported platforms are automatically recognised. Check with
```
  cd a8
  ./get_platform
```
