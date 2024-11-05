# Point LHAPDF to pdf sets
# export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-slc6-gcc8-opt/share/LHAPDF:/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.2.1/
#export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/share/LHAPDF:/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.5.1/ 
#export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc11-opt/share/LHAPDF:/cvmfs/cms.cern.ch/lhapdf/pdfsets/6.5.3/ 
# above taken from lcg releases website

#export TQZ_TOOLS_PATH='.'

# Enable newrt gcc and python
#source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_95 x86_64-slc6-gcc8-opt 
# source cvmfs from centos7 since HTCondor is set up for that -> Check and see if it works after jobs get over
# modify maybe the condor script to reflect and select same architecture
# source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh 
#source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh 
source /cvmfs/sft.cern.ch/lcg/views/LCG_105b/x86_64-centos7-gcc11-opt/setup.sh 
#### LCG_105 needed for RDataFrames full functionality - ROOT 6.30
#source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh 
