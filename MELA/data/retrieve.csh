#!/bin/tcsh -f
# Script to retrieve the MCFM library from the link specified in download.url

cd `dirname $0`/$1

set LIB="lib$2.so"
if (! -e "$LIB") then
  set filedir=`cat download.url`
  set MCFMfilename="$filedir$LIB"
  set Madgraphfilename="$filedir"libSMEFTsim.so
  wget --no-check-certificate -q "$filename" "$Madgraphfilename"
endif
