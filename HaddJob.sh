#!/bin/bash
echo "========================================="
fout=$1
fin=$2

wd="/user/yvekemans"
rel="MasterThesis/ZHTollSS_analysis"
reldir=$wd/$rel
cd $reldir
echo ""$reldir

#export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
#eval `scramv1 runtime -sh`
#source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "CMSSW intialized to "$CMSSW_BASE
source ./setup.sh

hadd -f $fout $fin/*.root
