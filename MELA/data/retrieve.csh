#!/bin/tcsh -f
# Script to retrieve the MCFM library from the link specified in download.url

cd `dirname $0`/$1

set LIB="lib$2.so"
if (! -e "$LIB") then
  set filedir=`cat download.url`
  set MCFMfilename="$filedir$LIB"
  wget --no-check-certificate -q "$MCFMfilename"
endif

set MADLIB="libMG_SMEFTsim_v$3.so"
if (! -e "$MADLIB") then
  set filedir=`cat download.url`
  set Madgraphfilename="$filedir$MADLIB"
  wget --no-check-certificate -q "$Madgraphfilename"
endif
