#!/bin/bash
echo "========================================="
cfg=$1
cuts=$2
output=$3
flow=$4
fhigh=$5
dname=$6
year=$7

echo "First arg: $cfg"
echo "Second arg: $cuts"
echo "Third arg: $output"
echo "Fourth arg: $flow"
echo "Fifth arg: $fhigh"
echo "Sixth arg: $dname"

wd="/user/sdansana"
rel="CMSSW_10_6_27/src/HToSS_analysis"
reldir=$wd/$rel
cd $reldir
#cd /user/sdansana/CMSSW_10_6_27/src/HToSS_analysis
echo ""$reldir


#export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
#eval `scramv1 runtime -sh`
#source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "CMSSW intialized to "$CMSSW_BASE
source ./setup.sh

python src/RDataFrames_analyzer.py -c $cfg --cuts $cuts -o $output -y $year --flow $flow --fhigh $fhigh --dname $dname
