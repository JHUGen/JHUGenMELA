#!/bin/bash

(
set -euo pipefail


getMELAARCH(){
  GCCVERSION=$(gcc -dumpversion)
  if [[ "$GCCVERSION" == "4.3"* ]] || [[ "$GCCVERSION" == "4.4"* ]] || [[ "$GCCVERSION" == "4.5"* ]]; then # v1 of MCFM library
    echo slc5_amd64_gcc434
  elif [[ "$GCCVERSION" == "4"* ]] || [[ "$GCCVERSION" == "5"* ]] || [[ "$GCCVERSION" == "6"* ]]; then # v2 of MCFM library
    echo slc6_amd64_gcc630
  elif [[ "$GCCVERSION" == "7"* ]]; then # v3 of MCFM library
    echo slc7_amd64_gcc700
  elif [[ "$GCCVERSION" == "8.0"* ]] || [[ "$GCCVERSION" == "8.1"* ]] || [[ "$GCCVERSION" == "8.2"* ]]; then # v4 of MCFM library
    echo slc7_amd64_gcc820
  elif [[ "$GCCVERSION" == "8"* ]]; then
    echo slc7_amd64_gcc830
  elif [[ "$GCCVERSION" == "12"* ]]; then
    echo el9_amd64_gcc12
  else
  #elif [[ "$GCCVERSION" == "9"* ]]; then
    echo slc7_amd64_gcc920
  #else
  fi
}

checkPYBIND11_INSTALL(){
  python3 -c "import pybind11" > /dev/null 2>&1 
  local status=$?
  echo $status
}

pyBIND11_STATUS=$(checkPYBIND11_INSTALL)


cd $(dirname ${BASH_SOURCE[0]})

MELADIR="$(readlink -f .)"
MCFMVERSION=mcfm_709
declare -i doDeps=0
declare -i doPrintEnv=0
declare -i doPrintEnvInstr=0
declare -i needROOFITSYS_ROOTSYS=0
declare -a setupArgs=()

for farg in "$@"; do
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "deps" ]]; then
    doDeps=1
  elif [[ "$fargl" == "env" ]]; then
    doPrintEnv=1
  elif [[ "$fargl" == "envinstr" ]]; then
    doPrintEnvInstr=1
  else
    setupArgs+=( "$farg" ) 
  fi
done
declare -i nSetupArgs
nSetupArgs=${#setupArgs[@]}

# Determine the MELA_ARCH string
# If CMS SCRAM_ARCH is set as an environment variable and it is present in data/,
# use SCRAM_ARCH as MELA_ARCH.
mela_arch=$(getMELAARCH)
mela_lib_path="${MELADIR}/data/${mela_arch}"

if [ "$pyBIND11_STATUS" != 0 ]; then
  echo "Cannot identify the python3 package PYBIND11. Please install the package or enter an area where it is installed."
  exit 1
fi

if [[ -z "${ROOFITSYS+x}" ]] && [[ $doDeps -eq 0 ]]; then
  if [[ $(ls ${ROOTSYS}/lib | grep -e libRooFitCore) != "" ]]; then
    needROOFITSYS_ROOTSYS=1
  else
    echo "Cannot identify ROOFITSYS. Please set this environment variable properly."
    exit 1
  fi
fi

printenv () {
  if [[ -z "${MELA_ARCH+x}" ]] || [[ "${MELA_ARCH}" != "${mela_arch}" ]]; then
    echo "export MELA_ARCH=${mela_arch}"
  fi

  if [[ -z "${MELA_LIB_PATH+x}" ]] || [[ "${MELA_LIB_PATH}" != "${mela_lib_path}" ]]; then
    echo "export MELA_LIB_PATH=${mela_lib_path}"
  fi

  ldlibappend="${mela_lib_path}"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$ldlibappend"* ]]; then
    echo "export LD_LIBRARY_PATH=${ldlibappend}${end}"
  fi

  pythonappend="${MELADIR}/python:${mela_lib_path}"
  end=""
  if [[ ! -z "${PYTHONPATH+x}" ]]; then
    end=":${PYTHONPATH}"
  fi
  if [[ "${end}" != *"$pythonappend"* ]]; then
    echo "export PYTHONPATH=${pythonappend}${end}"
  fi

  if [[ $needROOFITSYS_ROOTSYS -eq 1 ]]; then
    echo "export ROOFITSYS=${ROOTSYS}"
  fi
}
doenv () {
  if [[ -z "${MELA_ARCH+x}" ]] || [[ "${MELA_ARCH}" != "${mela_arch}" ]]; then
    export MELA_ARCH="${mela_arch}"
    echo "Temporarily using MELA_ARCH as ${MELA_ARCH}"
  fi

  if [[ -z "${MELA_LIB_PATH+x}" ]] || [[ "${MELA_LIB_PATH}" != "${mela_lib_path}" ]]; then
    export MELA_LIB_PATH="${mela_lib_path}"
    echo "Temporarily using MELA_LIB_PATH as ${MELA_LIB_PATH}"
  fi

  ldlibappend="${mela_lib_path}"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$ldlibappend"* ]]; then
    export LD_LIBRARY_PATH="${ldlibappend}${end}"
    echo "Temporarily using LD_LIBRARY_PATH as ${LD_LIBRARY_PATH}"
  fi

  pythonappend="${MELADIR}/python:${mela_lib_path}"
  end=""
  if [[ ! -z "${PYTHONPATH+x}" ]]; then
    end=":${PYTHONPATH}"
  fi
  if [[ "${end}" != *"$pythonappend"* ]]; then
    export PYTHONPATH="${pythonappend}${end}"
    echo "Temporarily using PYTHONPATH as ${PYTHONPATH}"
  fi

  if [[ $needROOFITSYS_ROOTSYS -eq 1 ]]; then
    export ROOFITSYS=${ROOTSYS}
    echo "Temporarily using ROOFITSYS as ${ROOTSYS}"
  fi
}
dodeps () {
  ${MELADIR}/COLLIER/setup.sh "${setupArgs[@]}"
  tcsh ${MELADIR}/data/retrieve.csh ${MELA_ARCH} $MCFMVERSION
  ${MELADIR}/downloadNNPDF.sh
}
printenvinstr () {
  echo
  echo "remember to do"
  echo
  echo 'eval $('${BASH_SOURCE[0]}' env)'
  echo "or"
  echo 'eval `'${BASH_SOURCE[0]}' env`'
  echo
  echo "if you are using a bash-related shell, or you can do"
  echo
  echo ${BASH_SOURCE[0]}' env'
  echo
  echo "and change the commands according to your shell in order to do something equivalent to set up the environment variables."
  echo
}

if [[ $doPrintEnv -eq 1 ]]; then
    printenv
    exit
elif [[ $doPrintEnvInstr -eq 1 ]]; then
    printenvinstr
    exit
fi

if [[ $nSetupArgs -eq 0 ]]; then
    setupArgs+=( -j 1 )
    nSetupArgs=2
fi


doenv

if [[ $doDeps -eq 1 ]]; then
    : ok
elif [[ "$nSetupArgs" -eq 1 ]] && [[ "${setupArgs[0]}" == *"clean"* ]]; then
    #echo "Cleaning C++"
    make clean

    #echo "Cleaning FORTRAN"
    pushd ${MELADIR}/fortran &> /dev/null
    make clean
    rm -f ../data/${MELA_ARCH}/libjhugenmela.so
    popd &> /dev/null

    #echo "Cleaning COLLIER"
    ${MELADIR}/COLLIER/setup.sh "${setupArgs[@]}"

    exit
elif [[ "$nSetupArgs" -ge 1 ]] && [[ "$nSetupArgs" -le 2 ]] && [[ "${setupArgs[0]}" == *"-j"* ]]; then
    : ok
else
    echo "Unknown arguments:"
    echo "  ${setupArgs[@]}"
    echo "Should be nothing, env, or clean"
    exit 1
fi


dodeps
if [[ $doDeps -eq 1 ]]; then
    exit
fi

pushd ${MELADIR}/fortran &> /dev/null
make "${setupArgs[@]}"
popd &> /dev/null

make "${setupArgs[@]}"
printenvinstr

)
