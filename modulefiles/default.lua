
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

setenv("AFC", "ifort")
setenv("ACC",  "cc -O")

set_alias("ae", "exe/as_exe")
