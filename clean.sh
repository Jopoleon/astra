#!/bin/bash -f

rootdir=`dirname $0`       # may be relative path
AWD=`cd $rootdir && pwd`  # ensure absolute path

cd $AWD

make -f exe/Makefile clean
