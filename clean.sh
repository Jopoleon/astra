#!/bin/bash

rootdir=`dirname $0`       # may be relative path
export AWD=`cd $rootdir && pwd`  # ensure absolute path

cd $AWD

make -f exe/Makefile clean
