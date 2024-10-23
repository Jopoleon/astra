
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
