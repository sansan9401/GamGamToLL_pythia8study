source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc900/cms/cmssw/CMSSW_11_3_0/src/
eval `scram runtime -sh`
cd -
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PWD/Plotter
