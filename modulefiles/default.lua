
-- -*- lua -*-

-- Name of the application
local AppName = "astra"

-- Version number
local Version = "4.1.1"

-- Sets module help message
help(
[[
    ASTRA is a suite of transport codes
]])


whatis("Name        : " .. AppName)
whatis("Version     : " .. Version)

conflict("env/gcc8.x-pgf20.11", "env/gcc8.x", "env/pgf20.11")
if(not isloaded("env/intel2020"))
then
    load("env/intel2020")
end
if(not isloaded("python/3.11"))
then
    load("python/3.11")
end

setenv("ASTRA_EXT","/fusion/projects/codes/astra/c8/8.2/intel2020/ASTRA_LIBRARIES_EXT")
setenv("MKL_LIBDIR", os.getenv("MKLROOT") .. "/lib/intel64")
setenv("COMP_LIBDIR", "/fusion/usc/opt/intel2020/compilers_and_libraries/linux/lib/intel64")
setenv("NETCDF_LIBDIR", os.getenv("NETCDF_DIR") .. "/lib")
setenv("NETCDF_LIB",  "-L" ..  os.getenv("NETCDF_LIBDIR") .. " -lnetcdf -lnetcdff -Wl,-rpath," ..  os.getenv("NETCDF_LIBDIR"))
setenv("NETCDF_INC",  os.getenv("NETCDF_DIR") .. "/include")

setenv("AFC", "ifort")
setenv("ACC",  "cc -O")
setenv("NETCDF_INC",  os.getenv("NETCDF_DIR") .. "/include")
setenv("XLBR",  "-L/usr/lib64 -lX11 -Wl,--start-group " .. os.getenv("COMP_LIBDIR") .. "/libintlc.so.5 " .. os.getenv("MKL_LIBDIR") .. "/libmkl_intel_lp64.a " .. os.getenv("MKL_LIBDIR") .. "/libmkl_core.a " .. os.getenv("MKL_LIBDIR") .. "/libmkl_sequential.a " .. os.getenv("MKL_LIBDIR") .. "/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group")

set_alias("ae", "exe/as_exe")
