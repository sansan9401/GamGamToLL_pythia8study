# PYTHIA configuration file.
# Generated on Mon Mar 22 12:05:55 CET 2021 with the user supplied options:
# --prefix=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/tmp/BUILDROOT/b182f211d06440c9dce5558bbfddbc88/opt/cmssw/slc7_amd64_gcc900/external/pythia8/303-ljfedo3
# --with-hepmc2=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/hepmc/2.06.10-ljfedo
# --with-lhapdf6=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/lhapdf/6.2.3-ljfedo2

# Install directory prefixes.
PREFIX_BIN=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/tmp/BUILDROOT/b182f211d06440c9dce5558bbfddbc88/opt/cmssw/slc7_amd64_gcc900/external/pythia8/303-ljfedo3/bin
PREFIX_INCLUDE=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/tmp/BUILDROOT/b182f211d06440c9dce5558bbfddbc88/opt/cmssw/slc7_amd64_gcc900/external/pythia8/303-ljfedo3/include
PREFIX_LIB=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/tmp/BUILDROOT/b182f211d06440c9dce5558bbfddbc88/opt/cmssw/slc7_amd64_gcc900/external/pythia8/303-ljfedo3/lib
PREFIX_SHARE=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/tmp/BUILDROOT/b182f211d06440c9dce5558bbfddbc88/opt/cmssw/slc7_amd64_gcc900/external/pythia8/303-ljfedo3/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/hepmc/2.06.10-ljfedo/
HEPMC2_INCLUDE=-I/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/hepmc/2.06.10-ljfedo/include
HEPMC2_LIB=-L/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/hepmc/2.06.10-ljfedo/lib -Wl,-rpath,/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/hepmc/2.06.10-ljfedo/lib -lHepMC

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/lhapdf/6.2.3-ljfedo2/bin/
LHAPDF6_INCLUDE=-I/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/lhapdf/6.2.3-ljfedo2/include
LHAPDF6_LIB=-L/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/lhapdf/6.2.3-ljfedo2/lib -Wl,-rpath,/data/cmsbld/jenkins/workspace/auto-builds/CMSSW_11_3_0_pre5-slc7_amd64_gcc900/build/CMSSW_11_3_0_pre5-build/slc7_amd64_gcc900/external/lhapdf/6.2.3-ljfedo2/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=lhapdf-config
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=lhapdf-config
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=/cvmfs/cms.cern.ch/slc7_amd64_gcc900/cms/cmssw/CMSSW_11_3_0/external/slc7_amd64_gcc900/bin/../../../../../../../slc7_amd64_gcc900/lcg/root/6.22.08-ljfedo/lib

YODA_USE=false
YODA_CONFIG=lhapdf-config
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=lhapdf-config
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=lhapdf-config
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=lhapdf-config
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=lhapdf-config
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=
