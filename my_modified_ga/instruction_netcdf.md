# Recompile NetCDF for intel fortran

- Download [netcdf-fortran](https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.0/netcdf-fortran-4.6.0.tar.gz)
- untar
- `cd` into the source
- `./configure --prefix="$(realpath ../netcdf)" FC=ifort FCFLAGS=-O3`
- `make -j10`
- `make install`
