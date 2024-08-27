Repository of the ASTRA code (G. Pereverzev, P. N. Yushmanov)

Clone:
```
  git clone https://gitlab.mpcdf.mpg.de/git/astra.git a8
```

Install (after clone or pull):
```
  cd a8
  chmod 744 install.sh
  ./install.sh
```

Compile or execute:
```
  cd a8
  exe/as_exe
```

Supported platforms:
  IPP tok
  IPP ldaug
  IPP lddsk
  IPP cs-ld-prod
  IPP-cz
  gateway
  iter-sdcc
  GA-iris
  GA-omega
  Perlmutter
  mit.edu
  rat2 (Padua)

The supported platforms are automatically recognised. If not, execute
```
  cd a8
  ./get_platform
```
