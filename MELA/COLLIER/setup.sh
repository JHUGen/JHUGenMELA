#!/bin/bash

{

set -euo pipefail

scriptdir=$(dirname $0)
curdir=$(pwd)

cd $scriptdir

pkgname="collier-1.2.0"
pkgdir="COLLIER-1.2"
tarname=$pkgname".tar.gz"
tarweb="https://collier.hepforge.org/downloads/"$tarname
libname="libcollier.so"
tmpdir="colliertmp"

if [[ $# > 0 ]] && [[ "$1" == *"clean"* ]]; then

  rm -f *.so
  rm -f *.o
  rm -f *.mod
  rm -f *.f
  rm -f *.F
  rm -f *.F90
  for f in $(ls ./);do
    if [[ -d $f ]]; then
      rm -rf $f
    fi
  done

  rm -f ${MELA_LIB_PATH}/../*/$libname

else

  if [[ ! -f "${MELA_LIB_PATH}/$libname" ]]; then
    rm -rf $tmpdir
    rm -f $tarname
    wget --no-check-certificate $tarweb
    mkdir $tmpdir
    tar -xvzf $tarname -C $tmpdir
    rm -f $tarname
    mv $tmpdir"/"$pkgdir"/src/"* ./
    rm -rf $tmpdir

    make $@
  fi

fi

cd $curdir

}
