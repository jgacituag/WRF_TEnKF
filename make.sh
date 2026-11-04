#!/bin/bash
#=======================================================================
# Build letkf.exe, update_wrf_time.exe and obs_ope.exe.
#
# Compiler settings for FUJITSU_FUGAKU, INTEL_FUGAKU and INTEL_HYDRA are
# copied from wrfutil_hydra/letkf_code/make_all.sh (commit 586cc96), and
# the compilation order from wrfutil_hydra/letkf_code/src/wrf/letkf/make.sh.
# GNU_LOCAL is only meant for syntax checks on a workstation.
# On Fugaku the NetCDF paths point to /vol0004/ra000007/data/jruiz, which
# is readable by the group; set NETCDF to override them.
#
# Usage: ./make.sh [FUJITSU_FUGAKU|INTEL_FUGAKU|INTEL_HYDRA|GNU_LOCAL]
# Output: dist/letkf_<COMPILATION_NAME>.tar (copy it to $LETKFDIR/letkf.tar)
#=======================================================================
set -e
export COMPILATION_NAME=${1:-FUJITSU_FUGAKU}
DEBUG=${DEBUG:-FALSE}
[ "$DEBUG" == "TRUE" ] && DEBUG_FLAGS=" -g -traceback "

case $COMPILATION_NAME in
FUJITSU_FUGAKU)
  PARF90="mpifrtpx"
  FFLAGS="-O3 -Nalloc_assign ${DEBUG_FLAGS}"
  OMP="-Kopenmp"
  NETCDF=${NETCDF:-"/vol0004/ra000007/data/jruiz/comp_libs_fujitsu/netcdf/"}   #same library as wrfutil machine.conf
  LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5  -lfjprofmpi -lmpi_cxx"
  INC_NETCDF="-I$NETCDF/include/"
  LBLAS="-SSL2BLAMP"
  ;;
INTEL_FUGAKU)
  . /opt/intel/oneapi/setvars.sh intel64
  PARF90="mpiifort"
  FFLAGS="-O3 -xHost -convert big_endian -FR ${DEBUG_FLAGS}"
  OMP="-qopenmp"
  NETCDF=${NETCDF:-"/vol0004/ra000007/data/jruiz/comp_libs/netcdf4/"}
  LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5"
  INC_NETCDF="-I$NETCDF/include/"
  LBLAS="-qmkl"
  ;;
INTEL_HYDRA)
  source /opt/load-libs.sh 1
  PARF90="mpif90"
  FFLAGS="-O3 -FR ${DEBUG_FLAGS}"
  OMP="-qopenmp"
  NETCDF="/opt/netcdf/netcdf-4/intel/2021.4.0/"
  HDF5="/opt/hdf5/hdf5-1.10.5/intel/2021.4.0/"
  LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf  -L$HDF5/lib -lhdf5_fortran -lhdf5_hl -lhdf5"
  INC_NETCDF="-I$NETCDF/include/"
  LBLAS="-qmkl"
  ;;
GNU_LOCAL)
  PARF90="mpif90"
  FFLAGS="-O2 -ffree-line-length-none -fallow-argument-mismatch -fconvert=big-endian ${DEBUG_FLAGS}"
  OMP="-fopenmp"
  INC_NETCDF="$(nf-config --fflags)"
  LIB_NETCDF="$(nf-config --flibs)"
  LBLAS="-llapack -lblas"
  ;;
*)
  echo "Unknown COMPILATION_NAME: $COMPILATION_NAME"; exit 1 ;;
esac

ROOT=$(cd "$(dirname "$0")" && pwd)
BUILD=$ROOT/build/$COMPILATION_NAME
rm -rf $BUILD && mkdir -p $BUILD $ROOT/dist
cd $BUILD
ln -s $ROOT/src/common/*.f90 $ROOT/src/common/*.f .
ln -s $ROOT/src/wrf_common/*.f90 .
ln -s $ROOT/src/letkf/*.f90 .

F="$PARF90 $OMP $FFLAGS"

# letkf.exe
$F -c SFMT.f90
$F -c common.f90
$F -c common_mpi.f90
$F -c common_mtx.f90
$PARF90 $OMP -O3 -c netlib.f
$F -c common_letkf.f90
$F -c module_map_utils.f90
$F $INC_NETCDF -c common_wrf.f90
$F $INC_NETCDF -c common_namelist.f90
$F -c common_smooth2d.f90
$F $INC_NETCDF -c common_obs_wrf.f90
$F $INC_NETCDF -c common_mpi_wrf.f90
$F $INC_NETCDF -c letkf_obs.f90
$F $INC_NETCDF -c letkf_tools.f90
$F $INC_NETCDF -c letkf.f90
$F $INC_NETCDF -o letkf.exe SFMT.o common.o common_mpi.o common_mtx.o netlib.o \
   common_letkf.o module_map_utils.o common_wrf.o common_namelist.o common_smooth2d.o \
   common_obs_wrf.o common_mpi_wrf.o letkf_obs.o letkf_tools.o letkf.o $LBLAS $LIB_NETCDF

# update_wrf_time.exe
$F $INC_NETCDF -c update_wrf_time.f90
$F $INC_NETCDF -o update_wrf_time.exe SFMT.o common.o module_map_utils.o common_wrf.o \
   update_wrf_time.o $LBLAS $LIB_NETCDF

# obs_ope.exe
$F $INC_NETCDF -c obs_ope.f90
$F $INC_NETCDF -o obs_ope.exe SFMT.o common.o common_mpi.o module_map_utils.o common_wrf.o \
   common_namelist.o common_smooth2d.o common_obs_wrf.o common_mpi_wrf.o letkf_obs.o \
   obs_ope.o $LBLAS $LIB_NETCDF

tar -cf $ROOT/dist/letkf_${COMPILATION_NAME}.tar letkf.exe update_wrf_time.exe obs_ope.exe
echo "Built $ROOT/dist/letkf_${COMPILATION_NAME}.tar"
