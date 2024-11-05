# ZHTollSS_analysis
This repository will be used for the extension of the HToSS_analysis to look at the ZH production mode. Currently, only a structure of the analysis is present - More information/parts to be imported from the HToSS_analysis as and when needed. This contains .py and .yaml files to run the analysis using RDataframes in pyROOT. 

## Setup the repo
1. Get the main branch:
```
git clone https://github.com/HToSS-Search/ZHTollSS_analysis.git
```
2. Checkout relevant branch (XXX = Yolan_dev / Soumya_dev): ``` git checkout -b XXX ```
3. For copying changes in local directory to remore repo, use:
```
git add .
git commit -a
git push
```

## Run the analysis
1. Source environment: ``` source setup.sh ```
2. For Z-mass peak using alternate (simplified) analysis (which uses JetMET study data):

For MC:   
```
python3 src/RDataFrames_Zpeak_alt.py -o plots/UL2017/MC.root -y UL2017 -c configs/2017/datasets/DYJetsToLL_M/DYJetsToLL_M-50_alt.yaml
```

For Data:   
```
python3 src/RDataFrames_Zpeak_alt.py -o plots/UL2017/Data.root -y UL2017 -c configs/2017/datasets/Data/mumuRun2017F_Zpeak.yaml
```

This produces a ROOT file. Use the right config file even though at the moment it is not used properly.

## Glossary
1. ``` configs ``` - Contains dataset locations + cuts yaml files used as two separate config files used to run the analysis
2. ``` plots ``` - Used typically to store final histogram root files
3. ``` scripts ``` - Can be used to store the scripts needed to run the analysis OR make plots (empty at the moment)
4. ``` scale_factors ``` - Contains files to scale the MC with Data/MC corrections
5. ``` pileup ``` - Pileup weights (currently not applied)
6. ``` src ``` - Contains the main analysis files
