#This is the RDF Analyzer for the ZH analysis
import numpy as np
from array import array
import os,argparse
import re
import yaml
import ROOT
import time
import sys

from ZH_analysis_helper import gInterpreter_PFIsolation, gInterpreter_lv, gInterpreter_getIndices, gInterpreter_getKinematics, gInterpreter_pairselection, gInterpreter_matching, gInterpreter_diObjectLxy, gInterpreter_std_map
from sf_helper import get_pileup, gInterpreter_SF

ROOT.gROOT.SetBatch(True)

cpu_count = 4 # give the same for request_cpus on condor script
ROOT.ROOT.EnableImplicitMT(cpu_count) #comment it when testing

def maxfilenumber(path):
    list_of_files = os.listdir(path)
    n = [int(re.findall("\d+", str(file))[0]) for file in list_of_files]
    return min(n),max(n)

start_time = time.time()
parser = argparse.ArgumentParser(description='Plot stacked histogram')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("--cuts", dest="cuts",   help="Enter cuts file to process", type=str)
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
parser.add_argument("-y","--year", dest="year", help="Year for processing", type=str)
parser.add_argument("--flow", dest="flow", help="start file number",default=1, type=int)
parser.add_argument("--fhigh", dest="fhigh", help="end file number", default=100,type=int)
parser.add_argument("--onlyweights", dest="onlyweights", help="just store MCweights weighted with cs*lumi", action='store_true')
parser.add_argument("--dname", dest="dname", help="stores the name of dataset - esp. needed for LLP reweighting", type=str)
# parser.add_argument("--total", dest="onlyweights", help="just store MCweights weighted with cs*lumi", action='store_true')

args = parser.parse_args()
data_name = args.dname

directories2 = []

fin = open(args.config,'r')
# f = open("params.txt",'w')
conf_pars = yaml.safe_load(fin)
# data_name = conf_pars['name']

lumi_scale = {'UL2016_APV': 19500, 'UL2016': 16800,'UL2017': 41480,'UL2018': 59830 } #in pb-1
lumi_factor = lumi_scale[args.year]


data_loc = conf_pars['locations']
cross_section = 1 if 'Run' in args.config else conf_pars['cross_section']
sum_wts = 1 if 'Run' in args.config else conf_pars['sum_weights']
# lumi = 1 if 'Run' in args.config else 41474 #2017 for now
# lumi = 1 if 'Run' in args.config else lumi_factor
# lumi = 1 if 'Run' in args.config else 4247.682053046 #2017D for now
isData = True if 'Run' in args.dname else False



###################### LOADING ALL FILES FOR PROCESSING ############################
isOldNtuple = False
if 'almorton' in data_loc:
    isOldNtuple = True
# data_loc = data_loc[0]
if data_loc[-1] != '/':
    data_loc = data_loc+'/'
sys.stderr.write(data_loc)
if isOldNtuple:
    directories = [data_loc] #below only for old ntuples
else:
    directories = [data_loc+d+"/" for d in os.listdir(data_loc) if os.path.isdir(os.path.join(data_loc, d))]

#sys.stderr.write(directories)
#sys.stderr.write(directories2)
# sys.stderr.write(directories+directories2)
directories=directories+directories2

treeName = "makeTopologyNtupleMiniAOD/tree"
list_of_files = []
for dirtmp in directories:
#   sys.stderr.write(dirtmp)
    flow, fhigh = maxfilenumber(dirtmp)
    if ('ZHTollSS' not in args.config):
        if flow > args.fhigh or fhigh < args.flow:
            continue
        if args.flow >= flow: 
            flow = args.flow
        if args.fhigh <= fhigh:
            fhigh = args.fhigh
    for i in range(flow, fhigh+1):
        fno = str(i)
        fistr = dirtmp+"output_"+fno+".root"
        if not os.path.exists(fistr):
            continue
        try:
            root_file = ROOT.TFile.Open(fistr)
            if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
                raise Exception(f"Error opening file: {fistr}")
            # else:
            #   continue
        except Exception as e:
            sys.stderr.write(f"Error processing file {fistr}: {e}")
            continue
        list_of_files.append(fistr)

#sys.stderr.write(str(list_of_files))

#Hardcoded MC Files for the sum of weights calculation
#list_of_files = ['/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_13.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_23.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_22.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_28.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_9.root']
 
###################### CALCULATION OF SUM OF WEIGHTS ####################
if args.onlyweights:
    # sum_wts calculated here
    sys.stderr.write("enters sum weights calculation")
    #sys.stderr.write(list_of_files)
    if not isData:
        print('weights...')
        if not isinstance(list_of_files, list):
            file = ROOT.TFile(list_of_files)
            weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
            weightPlot.SetDirectory(0)
            file.Close()
        else:
            print(list_of_files[0])
            file = ROOT.TFile(list_of_files[0])
            weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
            weightPlot.SetDirectory(0)
            file.Close()
            for i,fistr in enumerate(list_of_files):
                if i==0:
                    continue
                if not os.path.exists(fistr):
                    continue
                try:
                    # Open the ROOT file
                    file = ROOT.TFile.Open(fistr)
                    print(fistr)
                    # Check if the file was opened successfully
                    if not file or file.IsZombie() or file.TestBit(ROOT.TFile.kRecovered):
                        raise Exception(f"Error opening file: {fistr}")

                    # Process the file
                    tmpPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
                    # total_entries=total_entries+tmpPlot.GetEntries()
                    weightPlot.Add(tmpPlot)
                    # Close the file
                    file.Close()
                    # continue
                except Exception as e:
                    sys.stderr.write(f"Error processing file {fistr}: {e}")
                    continue  # Continue to the next file in case of an error
        totalEvents_ = weightPlot.GetBinContent(2) - weightPlot.GetBinContent(3) # bins filled from 1, but bins available from 0
        sum_wts = totalEvents_
    sys.stderr.write("\n\nsum of weights:"+str(sum_wts)+"\n")
    fout = ROOT.TFile(args.output,"RECREATE")
    weightPlot.Write()
    fout.Close()
    sys.stderr.write("Time taken: --- %s seconds ---" % (time.time() - start_time)+'\n')
    quit()


if not isinstance(list_of_files, list):
    file = ROOT.TFile(list_of_files)
    cutPlot = file.Get("makeTopologyNtupleMiniAOD/eventFilterAND").Clone()
    weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
    eventPlot = file.Get("makeTopologyNtupleMiniAOD/eventcount").Clone()
    weightPlot.SetDirectory(0)
    cutPlot.SetDirectory(0)
    eventPlot.SetDirectory(0)
    file.Close()
else:
    file = ROOT.TFile(list_of_files[0])
    cutPlot = file.Get("makeTopologyNtupleMiniAOD/eventFilterAND").Clone()
    weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
    eventPlot = file.Get("makeTopologyNtupleMiniAOD/eventcount").Clone()
    weightPlot.SetDirectory(0)
    cutPlot.SetDirectory(0)
    eventPlot.SetDirectory(0)
    file.Close()
    for i,fistr in enumerate(list_of_files):
        if i==0:
            continue
        if not os.path.exists(fistr):
            continue
        try:
            # Open the ROOT file
            file = ROOT.TFile.Open(fistr)

            # Check if the file was opened successfully
            if not file or file.IsZombie() or file.TestBit(ROOT.TFile.kRecovered):
                raise Exception(f"Error opening file: {fistr}")

            # Process the file
            tmpPlot = file.Get("makeTopologyNtupleMiniAOD/eventFilterAND").Clone()
            # total_entries=total_entries+tmpPlot.GetEntries()
            cutPlot.Add(tmpPlot)
            tmpPlot2 = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
            weightPlot.Add(tmpPlot2)
            tmpPlot3 = file.Get("makeTopologyNtupleMiniAOD/eventcount").Clone()
            eventPlot.Add(tmpPlot3)
            # Close the file
            file.Close()
            # continue
        except Exception as e:
            sys.stderr.write(f"Error processing file {fistr}: {e}")
            continue  # Continue to the next file in case of an error



###################### LOADING YAML FOR CUTS ############################
fcuts = open(args.cuts,'r')
# f = open("params.txt",'w')
fcut_pars = yaml.safe_load(fcuts)

cut_pars = fcut_pars['cuts']

if fcut_pars['hadronType'] == "kaon":
    chsMass_ = 0.493677 # 0.13957061;//pion mass
else:
    chsMass_ = 0.13957061

muonMass_ = 0.105 #GeV

#ADD TRIGGER CONDITIONS HERE WHEN AVAILABLE
#----

#UPDATE WITH RECOSCALAR ONCE IMPLEMENTED
mu_cuts = cut_pars['muons']
Z_cuts = cut_pars['recoZ']
ch_cuts = cut_pars['hadrons']
dihadron_cuts = cut_pars['dihadron']
recohiggs_cuts = cut_pars['recohiggs']
#recoscalar_cuts =cut_pars['recoscalar']

Zmass_low = str(Z_cuts['massLower'])
Zmass_high = str(Z_cuts['massUpper'])

diChPt_ = str(dihadron_cuts['pt'])
diChdR_ = str(dihadron_cuts['dR'])
diChFlag_ = dihadron_cuts['flag']


s1s2_mass_low = recohiggs_cuts['massLow']
s1s2_mass_high = recohiggs_cuts['massHigh']
unblind = recohiggs_cuts['unblind']
sideband = recohiggs_cuts['sideband']

hmass_low, hmass_high = 110, 140

#Higgs mass cut is always there; when flag is true, use tight mass conditions else use loose
#This is the peak region [122.5, 127.5]
higgs_mass_peak = 'higgs_mass >='+str(s1s2_mass_low)+' && '+'higgs_mass <='+str(s1s2_mass_high)

#This is the offpeak region ]-inf, 122.5[  U  ]127.5, +inf[
higgs_mass_peak_blinded = 'higgs_mass <'+str(s1s2_mass_low)+' || '+'higgs_mass >'+str(s1s2_mass_high)

#This is loosest peak region ]110, 140[
higgs_mass_peak_loose = f'higgs_mass >{hmass_low}'+' && '+f'higgs_mass <{hmass_high}'

#This is loose peak region ]120,130[
higgs_mass_peak_tight = f'higgs_mass > {s1s2_mass_low-sideband}'+' && '+f'higgs_mass < {s1s2_mass_high+sideband}'

#Higgs_mass_peak is tighter than higgs_mass_peak_tight
higgs_mass_cuts = higgs_mass_peak_tight
higgs_mass_cuts_loose = higgs_mass_peak_loose


#### Splitting into regions #####
s1lxysign,s2lxysign = fcut_pars['regions']['s1lxysign'],fcut_pars['regions']['s2lxysign']
R1 = f's1_lxysign_tmp <= {s2lxysign} && s2_lxysign_tmp <= {s1lxysign}'
R2h1h2 = f's1_lxysign_tmp > {s1lxysign} && s2_lxysign_tmp <= {s2lxysign}'
R2h3h4 = f's1_lxysign_tmp <= {s1lxysign} && s2_lxysign_tmp > {s2lxysign}'
R3 = f's1_lxysign_tmp > {s1lxysign} && s2_lxysign_tmp > {s2lxysign}'



###################### THE RDF ANALYZER ############################
gInterpreter_lv()
gInterpreter_getIndices()
gInterpreter_getKinematics() 
gInterpreter_pairselection()
gInterpreter_matching()
gInterpreter_diObjectLxy()
gInterpreter_PFIsolation()

#TESTING
#list_of_files_bis = ['/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_13.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_23.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_22.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_28.root','/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/ZH_HToSSTo4Hadrons_ZToLL_MH125_MS1p2_ctauS0_TuneCP2_13TeV-powheg-pythia8/RunIIUL17/241208/new/output_9.root']
#df = ROOT.RDataFrame(treeName, '/pnfs/iihe/cms/store/user/sdansana/HToSS/MC/nTuples/SingleMuon/Run2017F-UL2017_MiniAODv2-v1_DataUL2017F_ZH_production/250317_225506/0000/output_1.root') 
#df = ROOT.RDataFrame(treeName, list_of_files_bis)

#ORIGINAL
df = ROOT.RDataFrame(treeName, list_of_files)
sys.stderr.write('\nTree loaded in succesfully')

totalEntries = df.Count().GetValue()
sys.stderr.write(f"\nTotal Entries : {totalEntries}\n")


#Set up the weights
dataset_weight = cross_section/sum_wts

if "PUup" in args.output:
    shift="up"
elif "PUdown" in args.output:
    shift="down"
else:
    shift="nominal"

if isData:
    sys.stderr.write('\nLooking at Data.......')
    df = df.Define('weight', '1')
else:
    sys.stderr.write('\nLooking at MC.......')

    gInterpreter_std_map() # Accesses c++ map
    pileup_ratio,pileup_edges = get_pileup(args.year,shift)
    pileupMap = {e: r for e, r in zip(pileup_edges[:-1], pileup_ratio)}
    pileupMap_c = ROOT.GetTheMap()
    for key in pileupMap:
        pileupMap_c[key]=pileupMap[key]
    #Define the dataset weights MCweight * cross_section/sum_wts
    df = df.Define('weight_', 'processMCWeight')
    df = df.Define('weightOnlyDataset','weight_*{}'.format(dataset_weight))
    #Rescale dataset weight with pileup weights
    df = df.Define('PUReweight_sf','GetTheMap()[floor(numVert)]')
    df = df.Define('weight_tmp','weightOnlyDataset*PUReweight_sf')
    #CURRENTLY HARDCODE WEIGHT TO BE DATASETWEIGHT*PILEUP UNTIL FIXED
    df = df.Define('weight', 'weight_tmp')

#LATER ON IN CODE WILL UPDATE WEIGHTS WITH SCALE FACTORS WHERE NEEDED

##########################################################################################
#####################################Generated Objects####################################
##########################################################################################

if not isData:

########## Generated Muons ##########


    #Define muon pt of Z muons and all muons
    isMuZ = 'genParId == 13 && genParMotherId == 23'
    isAntiMuZ = 'genParId == -13 && genParMotherId == 23'
    isMuZComb = 'abs(genParId) == 13 && genParMotherId == 23'


    df_genMuZ = df.Define('genNegMuZ_pT', f'genParPt[{isMuZ}]')\
            .Define('genPosMuZ_pT', f'genParPt[{isAntiMuZ}]')
    #Set up kinematics for muonZ sorted according to descending pT
    df_genMuZ = df_genMuZ.Define('genMuZ_pT_unsorted',f'genParPt[{isMuZComb}]')\
            .Define('sort_indices', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genMuZ_pT_unsorted))')\
            .Define('genMuZ_pT', 'Take(genMuZ_pT_unsorted, sort_indices)')\
            .Define('genMuZ_sizeZ', 'genMuZ_pT.size()')\
            .Define('genMuZ_eta', f'Take(genParEta[{isMuZComb}], sort_indices)')\
            .Define('genMuZ_phi', f'Take(genParPhi[{isMuZComb}], sort_indices)')

    #Define leading and subleading muon Z pt
    df_genMuZ = df_genMuZ.Define('genMuZ_pT_leading', 'genNegMuZ_pT[0] > genPosMuZ_pT[0] ? genNegMuZ_pT[0]:genPosMuZ_pT[0]')\
            .Define('genMuZ_pT_subleading', 'genNegMuZ_pT[0] > genPosMuZ_pT[0] ? genPosMuZ_pT[0]:genNegMuZ_pT[0]')

    #Filter out the events with 2 or more genMuons from Z
    sys.stderr.write('Total events: ' + str(df.Count().GetValue()) + '\n')
    sys.stderr.write('Events with minimum 1 muons and muons are from Z: ' + str( df_genMuZ.Filter('genMuZ_sizeZ >= 1').Count().GetValue()) + '\n')
    sys.stderr.write('Events with minimum 2 muons and muons are from Z: ' + str(df_genMuZ.Filter('genMuZ_sizeZ >= 2').Count().GetValue()) + '\n')
    sys.stderr.write('Events with exactly 2 muons from Z: ' + str(df_genMuZ.Filter('genMuZ_sizeZ == 2').Count().GetValue()) + '\n')

    df_genMuZ = df_genMuZ.Filter('genMuZ_sizeZ >= 2')


    #Set up the invariant Z mass
    df_genMuZ = df_genMuZ.Define('mu_mass', str(muonMass_))\
            .Define('genZ_invmass','myDiLep(genMuZ_pT, genMuZ_eta, genMuZ_phi, mu_mass).M()')\
            .Define('genZ_pT', 'myDiLep(genMuZ_pT, genMuZ_eta, genMuZ_phi, mu_mass).Pt()')\
            .Define('genZ_deltaR', 'getDeltaR(genMuZ_eta, genMuZ_phi)')

    #Cut on the Z mass
    df_genMuZ_invcut = df_genMuZ.Filter('genZ_invmass < ' + Zmass_high + ' && genZ_invmass > ' + Zmass_low)


########## Generated Scalar to hh ##########

    sys.stderr.write('\n\nGen Muons\n')
    isScalar = 'abs(genParId) == 9000006'
    isK = 'abs(genParId) == 321'
    isKfromS1 = 'abs(genParId) == 321 && genParMotherId == 9000006'
    isKfromS2 = 'abs(genParId) == 321 && genParMotherId == -9000006'

    #Check genId
    df_genScalar = df_genMuZ.Define('genKId', f'genParId[{isK}]').Define('genK_size', 'genKId.size()')\
            .Define('genKS1Id', f'genParId[{isKfromS1}]')\
            .Define('genS1', f'genParMotherId[{isKfromS1}]')\
            .Define('genKS2Id', f'genParId[{isKfromS2}]')\
            .Define('genS2', f'genParMotherId[{isKfromS2}]')\
            .Define('genKS1_size', 'genKS1Id.size()').Define('genKS2_size', 'genKS2Id.size()')

    #Setup gen Kinematics
    df_genScalar = df_genScalar.Define('genKS1_pT_unsorted', f'genParPt[{isKfromS1}]').Define('genKS2_pT_unsorted', f'genParPt[{isKfromS2}]')\
            .Define('indicesKS1', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genKS1_pT_unsorted))')\
            .Define('indicesKS2', 'ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(genKS2_pT_unsorted))')\
            .Define('genKS1_pT', 'Take(genKS1_pT_unsorted, indicesKS1)').Define('genKS2_pT', 'Take(genKS2_pT_unsorted, indicesKS2)')\
            .Define('genKS1_eta', f'Take(genParEta[{isKfromS1}], indicesKS1)').Define('genKS2_eta', f'Take(genParEta[{isKfromS2}], indicesKS2)')\
            .Define('genKS1_phi', f'Take(genParPhi[{isKfromS1}], indicesKS1)').Define('genKS2_phi', f'Take(genParPhi[{isKfromS2}], indicesKS2)')\
            .Define('genKS1_pT_leading', 'genKS1_pT[0]').Define('genKS2_pT_leading', 'genKS2_pT[0]')\
            .Define('genKS1_pT_subleading', 'genKS1_pT[1]').Define('genKS2_pT_subleading', 'genKS2_pT[1]')

    #Filter out the events which have at least 2 hadrons from both scalars
    df_genScalar.Filter('genKS1_size >= 2 && genKS2_size >= 2')
    sys.stderr.write('scalar to hh events and Z to mumu in event: ' + str(df_genScalar.Count().GetValue()) + '\n')

    df_genScalar = df_genScalar.Define('genKS1_dPhi', 'ROOT::VecOps::DeltaPhi(genKS1_phi[0], genKS1_phi[1])')\
            .Define('genKS2_dPhi', 'ROOT::VecOps::DeltaPhi(genKS2_phi[0], genKS2_phi[1])')\
            .Define('genKS1_dR', 'ROOT::VecOps::DeltaR(genKS1_eta[0], genKS1_eta[1], genKS1_phi[0], genKS1_phi[1])')\
            .Define('genKS2_dR', 'ROOT::VecOps::DeltaR(genKS2_eta[0], genKS2_eta[1], genKS2_phi[0], genKS2_phi[1])')


    #Setup invariant mass
    df_genScalar = df_genScalar.Define('chs_mass', str(chsMass_))\
            .Define('genK1_lv', f'ROOT::Math::PtEtaPhiMVector(genKS1_pT[0], genKS1_eta[0], genKS1_phi[0], {chsMass_})')\
            .Define('genK2_lv', f'ROOT::Math::PtEtaPhiMVector(genKS1_pT[1], genKS1_eta[1], genKS1_phi[1], {chsMass_})')\
            .Define('genK3_lv', f'ROOT::Math::PtEtaPhiMVector(genKS2_pT[0], genKS2_eta[0], genKS2_phi[0], {chsMass_})')\
            .Define('genK4_lv', f'ROOT::Math::PtEtaPhiMVector(genKS2_pT[1], genKS2_eta[1], genKS2_phi[1], {chsMass_})')\
            .Define('genS1_lv', 'genK1_lv + genK2_lv')\
            .Define('genS2_lv', 'genK3_lv + genK4_lv')\
            .Define('genS1_invmass1', 'genS1_lv.M()')\
            .Define('genS2_invmass1', 'genS2_lv.M()')\
            .Define('genS1_invmass', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, chs_mass).M()')\
            .Define('genS1_pT', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, chs_mass).Pt()')\
            .Define('genS1_eta', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, chs_mass).Eta()')\
            .Define('genS1_phi', 'myDiLep(genKS1_pT, genKS1_eta, genKS1_phi, chs_mass).Phi()')\
            .Define('genS2_invmass', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, chs_mass).M()')\
            .Define('genS2_pT', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, chs_mass).Pt()')\
            .Define('genS2_eta', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, chs_mass).Eta()')\
            .Define('genS2_phi', 'myDiLep(genKS2_pT, genKS2_eta, genKS2_phi, chs_mass).Phi()')\
            .Define('genS12_dPhi', 'ROOT::VecOps::DeltaPhi(genS1_phi, genS2_phi)')\
            .Define('genS12_dR', 'ROOT::VecOps::DeltaR(genS1_eta, genS2_eta, genS1_phi, genS2_phi)')



##########################################################################################
###################################Reconstructed Objects##################################
##########################################################################################



##### Reconstructed Muons #####

sys.stderr.write('\n\nReconstructed Muons\n')
if not isData:
    sys.stderr.write('Events from gen Mu: ' + str(df_genMuZ.Count().GetValue()) + '\n')
    df_recoMu = df_genMuZ.Filter('numMuonPF2PAT >= 2', '2 or more recoMu')
else:
    sys.stderr.write('Total Events to start from before muon selection: ' + str(df.Count().GetValue()) + '\n')
    df_recoMu = df.Filter('numMuonPF2PAT >= 2', '2 or more recoMu')

sys.stderr.write('Events after numMuonPF2PAT >= 2 cut: ' + str(df_recoMu.Count().GetValue()) + '\n')
recoMu_qualitycut = 'abs(muonPF2PATEta) < ' + str(mu_cuts['eta']) + ' && muonPF2PATLooseCutId &&  muonPF2PATPt > ' + str(mu_cuts['pt'])
recoMu_qualitycutEta = 'abs(muonPF2PATEta) < ' + str(mu_cuts['eta'])
recoMu_qualitycutId = 'muonPF2PATLooseCutId'
recoMu_qualitycutPt = 'muonPF2PATPt > ' + str(mu_cuts['pt'])

#Take first two entries which correspond to highest pT muons, and some quality cuts on the muons
df_recoMu = df_recoMu.Define("indices", "ROOT::RVec<int>{0, 1}")\
        .Define('recoMu_pT', f'Take(muonPF2PATPt[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_eta', f'Take(muonPF2PATEta[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_phi', f'Take(muonPF2PATPhi[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_ch', f'Take(muonPF2PATCharge[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_RelIsoLoose', f'Take(muonPF2PATPfIsoLoose[{recoMu_qualitycut}], indices)')\
        .Define('recoMu_RelIsoTight', f'Take(muonPF2PATPfIsoTight[{recoMu_qualitycut}], indices)')

sys.stderr.write('PfRelIsoLoose pass: ' + str(df_recoMu.Filter('recoMu_RelIsoLoose[0] == 1 && recoMu_RelIsoLoose[1] == 1').Count().GetValue()) + '\n')
sys.stderr.write('PfRelIsoTight pass: ' + str(df_recoMu.Filter('recoMu_RelIsoTight[0] == 1 && recoMu_RelIsoTight[1] == 1').Count().GetValue()) + '\n')

df_recoMu = df_recoMu.Filter('Min(recoMu_pT) > ' + str(mu_cuts['pt']), 'Quality filter')
df_recoMu = df_recoMu.Define('recoMu_pT_leading', 'recoMu_pT[0]').Define('recoMu_pT_subleading', 'recoMu_pT[1]')
sys.stderr.write('Events after quality cuts: ' + str(df_recoMu.Count().GetValue()) + '\n')

#Only keep events with leading pT above cut and opposite charged muons
recoMu_cut = 'Max(recoMu_pT) > ' + str(mu_cuts['ptLeading']) + ' && recoMu_ch[0]*recoMu_ch[1] < 0'
df_recoMu_cut = df_recoMu.Filter(recoMu_cut, 'Muon cuts')
sys.stderr.write('Events after leadingPt, opposite charge cut: ' + str(df_recoMu_cut.Count().GetValue()) + '\n')

#Set up the invariant Z mass
df_recoMu_cut = df_recoMu_cut.Define('mu_mass_', str(muonMass_))\
        .Define('Z_mass', 'myDiLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass_).M()')\
        .Define('Z_pT', 'myDiLep(recoMu_pT, recoMu_eta, recoMu_phi, mu_mass_).Pt()')\
        .Define('Z_deltaR', 'getDeltaR(recoMu_eta, recoMu_phi)')

invMassCut = 'Z_mass < ' + Zmass_high + ' && Z_mass > ' + Zmass_low
df_recoMu_cut_invcut = df_recoMu_cut.Filter(invMassCut, 'Invariant mass cuts')
sys.stderr.write('Events after Zmass cut: ' + str(df_recoMu_cut_invcut.Count().GetValue()) + '\n')

#sys.stderr.write('CUTFLOWREPORT MUONS BELOW\n')
#sys.stderr.write(df_recoMu_cut_invcut.Report().Print())



##### Reconstructed hadrons ######

sys.stderr.write('\n\nReconstructed hadrons')

isRecoCh = 'abs(packedCandsPdgId) == 211 && packedCandsCharge!=0 && packedCandsHasTrackDetails==1' #reconstructed charged hadrons 

#cuts from muons
df_recoCh = df_recoMu_cut_invcut.Define('pCandChId', f'packedCandsPdgId[{isRecoCh}]')\
        .Define('recoCh_charge', f'packedCandsCharge[{isRecoCh}]')\
        .Define('recoCh_trkDetail', f'packedCandsHasTrackDetails[{isRecoCh}]')\
        .Define('chs_mass', str(chsMass_))\
        .Define('recoCh_px', f'packedCandsPx[{isRecoCh}]')\
        .Define('recoCh_py', f'packedCandsPy[{isRecoCh}]')\
        .Define('recoCh_pz', f'packedCandsPz[{isRecoCh}]')
        
#Check how many of the 3000 events have well reconstructed ch, keeping in mind all 3000 events have s->KK at gen level
df_recoCh_noMu = df.Define('pCandChIdNoMu', f'packedCandsPdgId[{isRecoCh}]')\
        .Define('recoChNoMu_charge', f'packedCandsCharge[{isRecoCh}]')

#Cuts on charged hadrons
ch_cuts =  'lv_trk_pt > ' + str(ch_cuts['pt']) + ' && ' + 'abs(lv_trk_eta) < ' + str(ch_cuts['eta'])

df_recoCh_cut = df_recoCh.Define('LVs', f'makeLVs(pCandChId, recoCh_charge, recoCh_trkDetail, recoCh_px, recoCh_py, recoCh_pz, {chsMass_}, true)')\
        .Define('lv_trk_pt','getKinematics(LVs, "pt")')\
        .Define('lv_trk_eta','getKinematics(LVs, "eta")')\
        .Define('ch_trk_pt',f'lv_trk_pt[{ch_cuts}]')\
        .Define('ch_trk_eta',f'lv_trk_eta[{ch_cuts}]')\
        .Define('LVs_sel',f'LVs[{ch_cuts}]')\
        .Define('ch_trk_sel','LVs_sel.size()>=4')

#Filter on the events which have 4 or more Charged Hadrons
df_recoCh_sel = df_recoCh_cut.Filter('ch_trk_sel')

#Set up the 2 hadron pairs, related to the two scalars
df_diCh = df_recoCh_sel.Define('cand_globalidx_hadronsonly',f'getIndices(numPackedCands)[{isRecoCh}]')\
        .Define('ch_globalidx1',f'cand_globalidx_hadronsonly[{ch_cuts}]')\
        .Define('ch1_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch2_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx1',f'''getDileptonCand1(ch_globalidx1, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch1_lv, ch2_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_globalidx2', 'ch_globalidx1[ch_globalidx1 != ch_pair_idx1[0] && ch_globalidx1 != ch_pair_idx1[1]]')\
        .Define('ch3_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch4_lv','ROOT::Math::PxPyPzMVector(0.,0.,0.,0.)')\
        .Define('ch_pair_idx2',f'''getDileptonCand1(ch_globalidx2, packedCandsCharge, packedCandsPx, packedCandsPy, packedCandsPz, {chsMass_}, numChsTrackPairs,chsTkPairIndex1,chsTkPairIndex2,chsTkPairTk1Px,chsTkPairTk1Py,chsTkPairTk1Pz,chsTkPairTk2Px,chsTkPairTk2Py,chsTkPairTk2Pz, "hadron", packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, ch3_lv, ch4_lv, {diChPt_}, {diChdR_})''')\
        .Define('ch_pair_check','ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0 && ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0')\
        .Define('ch_pair_checksingle','(ch_pair_idx1[0]>=0 && ch_pair_idx1[1]>=0) || (ch_pair_idx2[0]>=0 && ch_pair_idx2[1]>=0)')
       
df_diCh_check_single = df_diCh.Filter('ch_pair_checksingle')


#Seperate the kinematic variables for the hadron pairs
df_diCh_check = df_diCh.Filter('ch_pair_check')\
        .Define('ch1_pT_leading', 'ch1_lv.Pt()').Define('ch2_pT_subleading', 'ch2_lv.Pt()')\
        .Define('ch1_eta', 'ch1_lv.Eta()').Define('ch2_eta', 'ch2_lv.Eta()')\
        .Define('ch1_phi', 'ch1_lv.Phi()').Define('ch2_phi', 'ch2_lv.Phi()')\
        .Define('ch3_pT_leading', 'ch3_lv.Pt()').Define('ch4_pT_subleading', 'ch4_lv.Pt()')\
        .Define('ch3_eta', 'ch3_lv.Eta()').Define('ch4_eta', 'ch4_lv.Eta()')\
        .Define('ch3_phi', 'ch3_lv.Phi()').Define('ch4_phi', 'ch4_lv.Phi()')\
        .Define('ch1_px', 'ch1_lv.Px()').Define('ch1_py', 'ch1_lv.Py()').Define('ch1_pz', 'ch1_lv.Pz()')\
        .Define('ch2_px', 'ch2_lv.Px()').Define('ch2_py', 'ch2_lv.Py()').Define('ch2_pz', 'ch2_lv.Pz()')\
        .Define('ch3_px', 'ch3_lv.Px()').Define('ch3_py', 'ch3_lv.Py()').Define('ch3_pz', 'ch3_lv.Pz()')\
        .Define('ch4_px', 'ch4_lv.Px()').Define('ch4_py', 'ch4_lv.Py()').Define('ch4_pz', 'ch4_lv.Pz()')\
        .Define('ch1_charge', 'packedCandsCharge[ch_pair_idx1[0]]').Define('ch2_charge', 'packedCandsCharge[ch_pair_idx1[1]]')\
        .Define('ch3_charge', 'packedCandsCharge[ch_pair_idx2[0]]').Define('ch4_charge', 'packedCandsCharge[ch_pair_idx2[1]]')\
        .Define('ch1_nhits_tmp','packedCandsPseudoTrkNumberOfHits[ch_pair_idx1[0]]')\
        .Define('ch2_nhits_tmp','packedCandsPseudoTrkNumberOfHits[ch_pair_idx1[1]]')\
        .Define('ch3_nhits_tmp','packedCandsPseudoTrkNumberOfHits[ch_pair_idx2[0]]')\
        .Define('ch4_nhits_tmp','packedCandsPseudoTrkNumberOfHits[ch_pair_idx2[1]]')

#Define the angular differences
df_diCh_check = df_diCh_check.Define('ch12_dEta', 'ch1_eta - ch2_eta').Define('ch34_dEta', 'ch3_eta - ch4_eta')\
        .Define('ch12_dPhi', 'ROOT::VecOps::DeltaPhi(ch1_phi, ch2_phi)').Define('ch34_dPhi', 'ROOT::VecOps::DeltaPhi(ch3_phi, ch4_phi)')\
        .Define('ch12_dR', 'ROOT::VecOps::DeltaR(ch1_eta, ch2_eta, ch1_phi, ch2_phi)').Define('ch34_dR', 'ROOT::VecOps::DeltaR(ch3_eta, ch4_eta, ch3_phi, ch4_phi)')


sys.stderr.write('\nAmount of events passing muon cuts (no invmass cut): ' + str(df_recoMu_cut.Count().GetValue()))
#sys.stderr.write('\nAmount of events with ch reco larger than 4 per events: ' + str(df_recoCh_noMu.Filter('recoChNoMu_charge.size() >= 4').Count().GetValue()))
sys.stderr.write('\nAmount of events with >= 4 ch, events passed muon selection: ' + str(df_recoCh.Filter('recoCh_charge.size() >= 4').Count().GetValue()))
sys.stderr.write('\nAmount of events with >= 4 ch after hadron cuts (pT & Eta): ' + str(df_recoCh_sel.Count().GetValue()))
sys.stderr.write('\nAmount of events where single pair found: ' + str(df_diCh_check_single.Count().GetValue()))
sys.stderr.write('\nAmount of events where two pairs found: ' +  str(df_diCh_check.Count().GetValue()))

sys.stderr.write('\nnum of Ch events after quality cut  with 1 or more svVertexChi2: ' + str(df_diCh.Filter('svVertexChi2.size() >= 1').Count().GetValue()))
sys.stderr.write('\nnum of Ch events after quality cut with 2 or more svVertexChi2: ' +str(df_diCh.Filter('svVertexChi2.size() >= 2').Count().GetValue()))


if not isData:

############ GEN MATCHING ON THE PAIRS #############

    sys.stderr.write('\nGen matching of the reco chs pairs...')
    #Check if pairs can be found as gen matched pairs
    isnKfromS1 = 'genParId == -321 && genParMotherId == 9000006'
    ispKfromS1 = 'genParId == 321 && genParMotherId == 9000006'
    isnKfromS2 = 'genParId == -321 && genParMotherId == -9000006'
    ispKfromS2 = 'genParId == 321 && genParMotherId == -9000006'

    df_diCh_genmatch = df_diCh_check.Define('genParId_nKS1', f'genParId[{isnKfromS1}]').Define('genParId_pKS1', f'genParId[{ispKfromS1}]')\
            .Define('genParCharge_nKS1', f'genParCharge[{isnKfromS1}]').Define('genParCharge_pKS1', f'genParCharge[{ispKfromS1}]')\
            .Define('genParId_nKS2', f'genParId[{isnKfromS2}]').Define('genParId_pKS2', f'genParId[{ispKfromS2}]')\
            .Define('genParCharge_nKS2', f'genParCharge[{isnKfromS2}]').Define('genParCharge_pKS2', f'genParCharge[{ispKfromS2}]')

    #Get the neg and pos kaon kinematics from the scalars
    df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_pT', f'genParPt[{isnKfromS1}]').Define('nKS1_eta', f'genParEta[{isnKfromS1}]').Define('nKS1_phi', f'genParPhi[{isnKfromS1}]')\
            .Define('pKS1_pT', f'genParPt[{ispKfromS1}]').Define('pKS1_eta', f'genParEta[{ispKfromS1}]').Define('pKS1_phi', f'genParPhi[{ispKfromS1}]')\
            .Define('nKS2_pT', f'genParPt[{isnKfromS2}]').Define('nKS2_eta', f'genParEta[{isnKfromS2}]').Define('nKS2_phi', f'genParPhi[{isnKfromS2}]')\
            .Define('pKS2_pT', f'genParPt[{ispKfromS2}]').Define('pKS2_eta', f'genParEta[{ispKfromS2}]').Define('pKS2_phi', f'genParPhi[{ispKfromS2}]')

    #Make gen LVs
    df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_lv', 'ROOT::Math::PtEtaPhiMVector(nKS1_pT[0], nKS1_eta[0], nKS1_phi[0], 0.493677)')\
            .Define('pKS1_lv', 'ROOT::Math::PtEtaPhiMVector(pKS1_pT[0], pKS1_eta[0], pKS1_phi[0], 0.493677)')\
            .Define('nKS2_lv', 'ROOT::Math::PtEtaPhiMVector(nKS2_pT[0], nKS2_eta[0], nKS2_phi[0], 0.493677)')\
            .Define('pKS2_lv', 'ROOT::Math::PtEtaPhiMVector(pKS2_pT[0], pKS2_eta[0], pKS2_phi[0], 0.493677)')

    df_diCh_genmatch = df_diCh_genmatch.Define('chsPx', '''ROOT::VecOps::RVec<double> v = {ch1_px, ch2_px, ch3_px, ch4_px}; return v;''')\
            .Define('chsPy', '''ROOT::VecOps::RVec<double> v = {ch1_py, ch2_py, ch3_py, ch4_py}; return v;''')\
            .Define('chsPz', '''ROOT::VecOps::RVec<double> v = {ch1_pz, ch2_pz, ch3_pz, ch4_pz}; return v;''')\
            .Define('chsCharge', '''ROOT::VecOps::RVec<int> v = {ch1_charge, ch2_charge, ch3_charge, ch4_charge}; return v;''')\
            .Define('chsIdx', '''ROOT::VecOps::RVec<Int_t> v = {1, 2, 3, 4}; return v;''')\

    #Matching
    df_diCh_genmatch = df_diCh_genmatch.Define('nKS1_matchedIdx', 'MatchGenbis(nKS1_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, -1)')\
            .Define('pKS1_matchedIdx', 'MatchGenbis(pKS1_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, 1)')\
            .Define('nKS2_matchedIdx', 'MatchGenbis(nKS2_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, -1)')\
            .Define('pKS2_matchedIdx', 'MatchGenbis(pKS2_lv, chsPx, chsPy, chsPz, 0.493677, chsCharge, 1)')\

    #Filter out unmatched events
    df_diCh_genmatch_filtered = df_diCh_genmatch.Filter('nKS1_matchedIdx != -1 && pKS1_matchedIdx != -1 && nKS2_matchedIdx != -1 && pKS2_matchedIdx != -1 ')


    df_diCh_genmatch_filtered = df_diCh_genmatch_filtered.Define('genMatchPair1', '''ROOT::VecOps::RVec<Int_t> v = {chsIdx[nKS1_matchedIdx], chsIdx[pKS1_matchedIdx]}; return v;''')\
            .Define('genMatchPair2', '''ROOT::VecOps::RVec<Int_t> v = {chsIdx[nKS2_matchedIdx], chsIdx[pKS2_matchedIdx]}; return v;''')

    sys.stderr.write('\nEvents before matching conditions on pair selected events: ' + str(df_diCh_genmatch.Count().GetValue()))
    sys.stderr.write('\nAmount of events able to be genmatched: ' + str(df_diCh_genmatch_filtered.Count().GetValue()))



############ GEN MATCHING ON THE MUON SELECTED DATA #############

    sys.stderr.write('\n\nGen matching of the reco chs in the muon selected data...')
    df_genKMatch = df_recoMu_cut.Define('genParId_nKS1', f'genParId[{isnKfromS1}]').Define('genParId_pKS1', f'genParId[{ispKfromS1}]')\
            .Define('genParCharge_nKS1', f'genParCharge[{isnKfromS1}]').Define('genParCharge_pKS1', f'genParCharge[{ispKfromS1}]')\
            .Define('genParId_nKS2', f'genParId[{isnKfromS2}]').Define('genParId_pKS2', f'genParId[{ispKfromS2}]')\
            .Define('genParCharge_nKS2', f'genParCharge[{isnKfromS2}]').Define('genParCharge_pKS2', f'genParCharge[{ispKfromS2}]')

    #Get the neg and pos kaon kinematics from the scalars
    df_genKMatch = df_genKMatch.Define('nKS1_pT', f'genParPt[{isnKfromS1}]').Define('nKS1_eta', f'genParEta[{isnKfromS1}]').Define('nKS1_phi', f'genParPhi[{isnKfromS1}]')\
            .Define('pKS1_pT', f'genParPt[{ispKfromS1}]').Define('pKS1_eta', f'genParEta[{ispKfromS1}]').Define('pKS1_phi', f'genParPhi[{ispKfromS1}]')\
            .Define('nKS2_pT', f'genParPt[{isnKfromS2}]').Define('nKS2_eta', f'genParEta[{isnKfromS2}]').Define('nKS2_phi', f'genParPhi[{isnKfromS2}]')\
            .Define('pKS2_pT', f'genParPt[{ispKfromS2}]').Define('pKS2_eta', f'genParEta[{ispKfromS2}]').Define('pKS2_phi', f'genParPhi[{ispKfromS2}]')

    #Make gen LVs
    df_genKMatch = df_genKMatch.Define('nKS1_lv', 'ROOT::Math::PtEtaPhiMVector(nKS1_pT[0], nKS1_eta[0], nKS1_phi[0], 0.493677)')\
            .Define('pKS1_lv', 'ROOT::Math::PtEtaPhiMVector(pKS1_pT[0], pKS1_eta[0], pKS1_phi[0], 0.493677)')\
            .Define('nKS2_lv', 'ROOT::Math::PtEtaPhiMVector(nKS2_pT[0], nKS2_eta[0], nKS2_phi[0], 0.493677)')\
            .Define('pKS2_lv', 'ROOT::Math::PtEtaPhiMVector(pKS2_pT[0], pKS2_eta[0], pKS2_phi[0], 0.493677)')

    #Get the matched indices of the reco Kaons
    df_recoCh_matched = df_genKMatch.Define('nKS1_matchedIdx', 'MatchGen(nKS1_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, -1)')\
            .Define('pKS1_matchedIdx', 'MatchGen(pKS1_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, 1)')\
            .Define('nKS2_matchedIdx', 'MatchGen(nKS2_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, -1)')\
            .Define('pKS2_matchedIdx', 'MatchGen(pKS2_lv, packedCandsPdgId, packedCandsPx, packedCandsPy, packedCandsPz, 0.493677, packedCandsCharge, 1)')

    sys.stderr.write('\nEvents before matching conditions: ' + str(df_recoCh_matched.Count().GetValue()))

    #Filter out unmatched events
    df_recoCh_matched_filtered = df_recoCh_matched.Filter('nKS1_matchedIdx != -1 && pKS1_matchedIdx != -1 && nKS2_matchedIdx != -1 && pKS2_matchedIdx != -1 ')
    #Create a df with the unmatched events
    df_recoCh_nomatch = df_recoCh_matched.Filter('nKS1_matchedIdx == -1 || pKS1_matchedIdx == -1 || nKS2_matchedIdx == -1 || pKS2_matchedIdx == -1 ')

    df_recoCh_nomatch = df_recoCh_nomatch.Define('genKS1NoMatch_absEta', 'abs(nKS1_eta[0]) > abs(pKS1_eta[0]) ? abs(nKS1_eta[0]) : abs(pKS1_eta[0])')\
            .Define('genKS2NoMatch_absEta', 'abs(nKS2_eta[0]) > abs(pKS2_eta[0]) ? abs(nKS2_eta[0]) : abs(pKS2_eta[0])')\
            .Define('genKS1NoMatch_pT_leading', 'nKS1_pT[0] > pKS1_pT[0] ? nKS1_pT[0] : pKS1_pT[0]')\
            .Define('genKS1NoMatch_pT_subleading', 'nKS1_pT[0] < pKS1_pT[0] ? nKS1_pT[0] : pKS1_pT[0]')\
            .Define('genKS2NoMatch_pT_leading', 'nKS2_pT[0] > pKS2_pT[0]? nKS2_pT[0] : pKS2_pT[0]')\
            .Define('genKS2NoMatch_pT_subleading', 'nKS2_pT[0] < pKS2_pT[0] ? nKS2_pT[0] : pKS2_pT[0]')


    #Check if the gen pairs are track pairs, check crossover as well
    df_recoCh_matched_TkCheck = df_recoCh_matched_filtered.Define('TkPair1Check', f'checkTrkPair(nKS1_matchedIdx, pKS1_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
            .Define('TkPair2Check', f'checkTrkPair(nKS2_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
            .Define('TkPairCrossCheck_n1p2', f'checkTrkPair(nKS1_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
            .Define('TkPairCrossCheck_n2p1', f'checkTrkPair(nKS2_matchedIdx, pKS1_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
            .Define('TkPairCrossCheck_n1n2', f'checkTrkPair(nKS1_matchedIdx, nKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')\
            .Define('TkPairCrossCheck_p1p2', f'checkTrkPair(pKS1_matchedIdx, pKS2_matchedIdx, numChsTrackPairs, chsTkPairIndex1, chsTkPairIndex2)')
     
    sys.stderr.write('\nNumber of unmatched events: ' + str(df_recoCh_nomatch.Count().GetValue()))
    sys.stderr.write('\nNumber of matched events: ' + str(df_recoCh_matched_filtered.Count().GetValue()))
    sys.stderr.write('\nNumber of matched events and one correct Tkpair found: ' + str(df_recoCh_matched_TkCheck.Filter('TkPair1Check != -1 || TkPair2Check != -1').Count().GetValue()))
    sys.stderr.write('\nNumber of matched events and both correct Tkpair found: ' + str(df_recoCh_matched_TkCheck.Filter('TkPair1Check != -1 && TkPair2Check != -1').Count().GetValue()))
    sys.stderr.write('\n-----------CrossPairs-----------')
    sys.stderr.write('\nNumber of matched events where n1p2 or n2p1 is a track pair: ' + str(df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1p2 != -1 || TkPairCrossCheck_n2p1 != -1').Count().GetValue()))
    sys.stderr.write('\nNumber of matched events where n1p2 and n2p1 is a track pair: ' + str(df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1p2 != -1 && TkPairCrossCheck_n2p1 != -1').Count().GetValue()))
    sys.stderr.write('\nNumber of matched events where n1n2 or p1p2 is a track pair: ' + str(df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1n2 != -1 || TkPairCrossCheck_p1p2 != -1').Count().GetValue()))
    sys.stderr.write('\nNumber of matched events where n1n2 andp1p2 is a track pair: ' + str(df_recoCh_matched_TkCheck.Filter('TkPairCrossCheck_n1n2 != -1 && TkPairCrossCheck_p1p2 != -1').Count().GetValue()))
    sys.stderr.write('\n-------------------------------')
    sys.stderr.write('\nnum of matched gen events with 1 or more svVertexChi2: ' + str(df_recoCh_matched_filtered.Filter('svVertexChi2.size() >= 1').Count().GetValue()))
    sys.stderr.write('\nnum of matched gen events with 2 or more svVertexChi2: ' + str(df_recoCh_matched_filtered.Filter('svVertexChi2.size() >= 2').Count().GetValue()))
    sys.stderr.write('\nnum of total events with 1 or more svVertexChi2: ' + str(df.Filter('svVertexChi2.size() >= 1').Count().GetValue()))
    sys.stderr.write('\nnum of total events with 2 or more svVertexChi2: ' + str(df.Filter('svVertexChi2.size() >= 2').Count().GetValue()))

    #Select the matched reco px, py, pz
    df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_px', 'packedCandsPx[nKS1_matchedIdx]')\
            .Define('nKS1Matched_py', 'packedCandsPy[nKS1_matchedIdx]').Define('nKS1Matched_pz', 'packedCandsPz[nKS1_matchedIdx]')\
            .Define('pKS1Matched_px', 'packedCandsPx[pKS1_matchedIdx]').Define('pKS1Matched_py', 'packedCandsPy[pKS1_matchedIdx]').Define('pKS1Matched_pz', 'packedCandsPz[pKS1_matchedIdx]')\
            .Define('nKS2Matched_px', 'packedCandsPx[nKS2_matchedIdx]').Define('nKS2Matched_py', 'packedCandsPy[nKS2_matchedIdx]').Define('nKS2Matched_pz', 'packedCandsPz[nKS2_matchedIdx]')\
            .Define('pKS2Matched_px', 'packedCandsPx[pKS2_matchedIdx]').Define('pKS2Matched_py', 'packedCandsPy[pKS2_matchedIdx]').Define('pKS2Matched_pz', 'packedCandsPz[pKS2_matchedIdx]')

    #Make matched reco LVs
    df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_lv', 'ROOT::Math::PxPyPzMVector(nKS1Matched_px, nKS1Matched_py, nKS1Matched_pz, 0.493677)')\
            .Define('pKS1Matched_lv', 'ROOT::Math::PxPyPzMVector(pKS1Matched_px, pKS1Matched_py, pKS1Matched_pz, 0.493677)')\
            .Define('nKS2Matched_lv', 'ROOT::Math::PxPyPzMVector(nKS2Matched_px, nKS2Matched_py, nKS2Matched_pz, 0.493677)')\
            .Define('pKS2Matched_lv', 'ROOT::Math::PxPyPzMVector(pKS2Matched_px, pKS2Matched_py, pKS2Matched_pz, 0.493677)')

    #Get pT, Eta, Phi from the mathced reco LVs
    df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('nKS1Matched_pT', 'nKS1Matched_lv.Pt()')\
            .Define('nKS1Matched_eta', 'nKS1Matched_lv.Eta()').Define('nKS1Matched_phi', 'nKS1Matched_lv.Phi()')\
            .Define('pKS1Matched_pT', 'pKS1Matched_lv.Pt()').Define('pKS1Matched_eta', 'pKS1Matched_lv.Eta()').Define('pKS1Matched_phi', 'pKS1Matched_lv.Phi()')\
            .Define('nKS2Matched_pT', 'nKS2Matched_lv.Pt()').Define('nKS2Matched_eta', 'nKS2Matched_lv.Eta()').Define('nKS2Matched_phi', 'nKS2Matched_lv.Phi()')\
            .Define('pKS2Matched_pT', 'pKS2Matched_lv.Pt()').Define('pKS2Matched_eta', 'pKS2Matched_lv.Eta()').Define('pKS2Matched_phi', 'pKS2Matched_lv.Phi()')\
            .Define('S1_lv', 'nKS1Matched_lv + pKS1Matched_lv').Define('S2_lv', 'nKS2Matched_lv + pKS2Matched_lv')\
            .Define('s1Matched_invmass', 'S1_lv.M()').Define('s2Matched_invmass', 'S2_lv.M()')\
            .Define('s1Matched_pT', 'S1_lv.Pt()').Define('s2Matched_pT', 'S2_lv.Pt()')


    #Seperate leading and subleading contributions
    df_recoCh_matched_filtered = df_recoCh_matched_filtered.Define('ch1Matched_pT_leading', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_pT:pKS1Matched_pT')\
            .Define('ch2Matched_pT_subleading', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_pT:nKS1Matched_pT')\
            .Define('ch1Matched_eta', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_eta:pKS1Matched_eta')\
            .Define('ch2Matched_eta', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_eta:nKS1Matched_eta')\
            .Define('ch1Matched_phi', 'nKS1Matched_pT > pKS1Matched_pT ? nKS1Matched_phi:pKS1Matched_phi')\
            .Define('ch2Matched_phi', 'nKS1Matched_pT > pKS1Matched_pT ? pKS1Matched_phi:nKS1Matched_phi')\
            .Define('ch3Matched_pT_leading', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_pT:pKS2Matched_pT')\
            .Define('ch4Matched_pT_subleading', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_pT:nKS2Matched_pT')\
            .Define('ch3Matched_eta', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_eta:pKS2Matched_eta')\
            .Define('ch4Matched_eta', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_eta:nKS2Matched_eta')\
            .Define('ch3Matched_phi', 'nKS2Matched_pT > pKS2Matched_pT ? nKS2Matched_phi:pKS2Matched_phi')\
            .Define('ch4Matched_phi', 'nKS2Matched_pT > pKS2Matched_pT ? pKS2Matched_phi:nKS2Matched_phi')


    #Check quality cut on the matched reco ch
    sys.stderr.write('\nPerforming pair selection on the matched ch')
    df_recoCh_matched_pairsel = df_recoCh_matched_filtered.Define('LVs','''ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> v = {nKS1Matched_lv, pKS1Matched_lv, nKS2Matched_lv, pKS2Matched_lv}; return v;''' )\
            .Define('lv_trk_pt','getKinematics(LVs, "pt")')\
            .Define('lv_trk_eta','getKinematics(LVs, "eta")')\
            .Define('ch_trk_pt',f'lv_trk_pt[{ch_cuts}]')\
            .Define('ch_trk_eta',f'lv_trk_eta[{ch_cuts}]')\
            .Define('LVs_sel',f'LVs[{ch_cuts}]')\
            .Define('ch_trk_sel','LVs_sel.size()>=4')


    #Filter on the events which have 4 or more Charged Hadrons
    df_recoCh_matched_pairsel_qualcut = df_recoCh_matched_pairsel.Filter('ch_trk_sel')

    #sys.stderr.write('\nAmount of events with ch reco larger than 4 per events ' + str(df_recoCh_noMu.Filter('recoChNoMu_charge.size() >= 4').Count().GetValue()))
    sys.stderr.write('\nAmount of matched events with >= 4 ch after quality cuts (pT & Eta): ' + str(df_recoCh_matched_pairsel_qualcut.Count().GetValue()))

    sys.stderr.write('\nnum of matched events after quality cut  with 1 or more svVertexChi2: ' +  str(df_recoCh_matched_pairsel_qualcut.Filter('svVertexChi2.size() >= 1').Count().GetValue()))
    sys.stderr.write('\nnum of matched events after quality cut with 2 or more svVertexChi2: ' + str(df_recoCh_matched_pairsel_qualcut.Filter('svVertexChi2.size() >= 2').Count().GetValue()))



##### Reconstructed Scalar/Higss ######

sys.stderr.write('\n\nReconstructed Scalar and Higgs')

#Scalar kinematics
df_s1s2 = df_diCh_check.Define('s1_lv', 'ch1_lv + ch2_lv' ).Define('s2_lv', 'ch3_lv + ch4_lv')\
        .Define('s1_mass', 's1_lv.M()').Define('s2_mass', 's2_lv.M()')\
        .Define('s1_pT', 's1_lv.Pt()').Define('s2_pT', 's2_lv.Pt()')\
        .Define('s12_mass', '''double s = (s1_mass + s2_mass)/2.0; return s;''')\
        .Define('s12_dPhi', 'ROOT::VecOps::DeltaPhi(s1_lv.Phi(), s2_lv.Phi())')\
        .Define('s12_dR', 'ROOT::VecOps::DeltaR(s1_lv.Eta(), s2_lv.Eta(), s1_lv.Phi(), s2_lv.Phi())')

sys.stderr.write('\nNo scalar cuts currently')


############## LXY CALCULATION FOR SCALARS ########################
pv_sel = 'pvChi2!=0 && pvNdof!=0'
df_vertex = df_s1s2.Define('PVCov00',f'pvCov00[{pv_sel}][0]').Define('PVCov01',f'pvCov01[{pv_sel}][0]').Define('PVCov02',f'pvCov02[{pv_sel}][0]').Define('PVCov10',f'pvCov10[{pv_sel}][0]').Define('PVCov11',f'pvCov11[{pv_sel}][0]').Define('PVCov12',f'pvCov12[{pv_sel}][0]').Define('PVCov20',f'pvCov20[{pv_sel}][0]').Define('PVCov21',f'pvCov21[{pv_sel}][0]').Define('PVCov22',f'pvCov22[{pv_sel}][0]').Define('PVX',f'pvX[{pv_sel}][0]').Define('PVY',f'pvY[{pv_sel}][0]').Define('PVZ',f'pvZ[{pv_sel}][0]')\
                        .Define('s1_lxyInfo','''getLxy(PVCov00,PVCov01,PVCov02,PVCov10,PVCov11,PVCov12,PVCov20,PVCov21,PVCov22,PVX,PVY,
                                        chsTkPairTkVtxCov00[ch_pair_idx1[2]],chsTkPairTkVtxCov01[ch_pair_idx1[2]],chsTkPairTkVtxCov02[ch_pair_idx1[2]],
                                        chsTkPairTkVtxCov10[ch_pair_idx1[2]],chsTkPairTkVtxCov11[ch_pair_idx1[2]],chsTkPairTkVtxCov12[ch_pair_idx1[2]],
                                        chsTkPairTkVtxCov20[ch_pair_idx1[2]],chsTkPairTkVtxCov21[ch_pair_idx1[2]],chsTkPairTkVtxCov22[ch_pair_idx1[2]],
                                        chsTkPairTkVx[ch_pair_idx1[2]],chsTkPairTkVy[ch_pair_idx1[2]])''')\
                        .Define('s2_lxyInfo','''getLxy(PVCov00,PVCov01,PVCov02,PVCov10,PVCov11,PVCov12,PVCov20,PVCov21,PVCov22,PVX,PVY,
                                        chsTkPairTkVtxCov00[ch_pair_idx2[2]],chsTkPairTkVtxCov01[ch_pair_idx2[2]],chsTkPairTkVtxCov02[ch_pair_idx2[2]],
                                        chsTkPairTkVtxCov10[ch_pair_idx2[2]],chsTkPairTkVtxCov11[ch_pair_idx2[2]],chsTkPairTkVtxCov12[ch_pair_idx2[2]],
                                        chsTkPairTkVtxCov20[ch_pair_idx2[2]],chsTkPairTkVtxCov21[ch_pair_idx2[2]],chsTkPairTkVtxCov22[ch_pair_idx2[2]],
                                        chsTkPairTkVx[ch_pair_idx2[2]],chsTkPairTkVy[ch_pair_idx2[2]])''')\
                        .Define('s1_lxy_vector','ROOT::Math::XYVector(chsTkPairTkVx[ch_pair_idx1[2]] - PVX, chsTkPairTkVy[ch_pair_idx1[2]] - PVY)')\
                        .Define('s2_lxy_vector','ROOT::Math::XYVector(chsTkPairTkVx[ch_pair_idx2[2]] - PVX, chsTkPairTkVy[ch_pair_idx2[2]] - PVY)')\
                        .Define('pt_dphi','ROOT::Math::VectorUtil::DeltaPhi(s1_lv,s2_lv)')\
                        .Define('lxy_dphi','ROOT::Math::VectorUtil::DeltaPhi(s1_lxy_vector,s2_lxy_vector)')\
                        .Define('s1_vz_tmp','chsTkPairTkVz[ch_pair_idx1[2]]')\
                        .Define('s1_chi2_tmp','chsTkPairTkVtxChi2[ch_pair_idx1[2]]')\
                        .Define('s1_ndof_tmp','chsTkPairTkVtxNdof[ch_pair_idx1[2]]')\
                        .Define('s1_dz_tmp','chsTkPairTkVz[ch_pair_idx1[2]] - PVZ')\
                        .Define('s2_vz_tmp','chsTkPairTkVz[ch_pair_idx2[2]]')\
                        .Define('s2_chi2_tmp','chsTkPairTkVtxChi2[ch_pair_idx2[2]]')\
                        .Define('s2_ndof_tmp','chsTkPairTkVtxNdof[ch_pair_idx2[2]]')\
                        .Define('s2_dz_tmp','chsTkPairTkVz[ch_pair_idx2[2]] - PVZ')\
                        .Define('s1_chi2ndof_tmp','s1_chi2_tmp/s1_ndof_tmp')\
                        .Define('s2_chi2ndof_tmp','s2_chi2_tmp/s2_ndof_tmp')\
                        .Define('s1_nhits_tmp','abs(ch2_nhits_tmp-ch1_nhits_tmp)')\
                        .Define('s2_nhits_tmp','abs(ch4_nhits_tmp-ch3_nhits_tmp)')\
                        .Define('s1_lxysign_tmp','s1_lxyInfo[2]')\
                        .Define('s2_lxysign_tmp','s2_lxyInfo[2]')\
                        .Define('s1_sim_val','s1_lxyInfo[3]')\
                        .Define('s2_sim_val','s2_lxyInfo[3]')\
                        .Define('s1_lxy','s1_lxyInfo[0]')\
                        .Define('s2_lxy','s2_lxyInfo[0]')

#NO CUTS ON SCALAR MASS CURRENTLY


#Higgs kinematics
df_higgs_preblind = df_vertex.Define('higgs_lv', 's1_lv + s2_lv')\
        .Define('higgs_mass', 'higgs_lv.M()')\
        .Define('higgs_pT', 'higgs_lv.Pt()')\
        .Define('recohiggs_mass_peak', higgs_mass_peak)\
        .Define('recohiggs_mass_blinded', higgs_mass_peak_blinded)\
        .Define('recohiggs_mass_check_loose', higgs_mass_cuts_loose)\
        .Define('recohiggs_mass_check', higgs_mass_cuts)\
        .Define('ch1_iso','''PFIsolation("hadron", ch1_lv, ch_pair_idx1[0], ch_pair_idx1[1], packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, 0.4)''')\
        .Define('ch2_iso','''PFIsolation("hadron", ch2_lv, ch_pair_idx1[1], ch_pair_idx1[0], packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, 0.4)''')\
        .Define('ch3_iso','''PFIsolation("hadron", ch3_lv, ch_pair_idx2[0], ch_pair_idx2[1], packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, 0.4)''')\
        .Define('ch4_iso','''PFIsolation("hadron", ch4_lv, ch_pair_idx2[1], ch_pair_idx2[0], packedCandsPx, packedCandsPy, packedCandsPz, packedCandsE, packedCandsCharge, packedCandsPdgId, packedCandsFromPV, numPackedCands, 0.4)''')\
        .Define('ch1_reliso', 'ch1_iso[0]')\
        .Define('ch2_reliso', 'ch2_iso[0]')\
        .Define('ch3_reliso', 'ch3_iso[0]')\
        .Define('ch4_reliso', 'ch4_iso[0]')

df_higgs_preblind_loose = df_higgs_preblind.Filter('recohiggs_mass_check_loose', f'm(s1s2) in ]{hmass_low},{hmass_high}[')

sys.stderr.write('\nAmount of Events after scalar selection/cuts: ' + str(df_higgs_preblind.Count().GetValue()))
sys.stderr.write('\nAmount of Events with m(s1s2)/m(higgs) in ]110, 140[: ' + str(df_higgs_preblind_loose.Count().GetValue()))
sys.stderr.write('\nAmount of Events with m(s1s2)/m(higgs) in [-inf, 110] U [140, +inf[: ' + str(df_higgs_preblind.Filter(f'higgs_mass <= {hmass_low}'+' || '+f'higgs_mass >= {hmass_high}').Count().GetValue()))


############## Scale Factors ########################
#HERE SCALEFACTORS ARE ADDED TO DF

df_higgs_preblind_loose = df_higgs_preblind_loose.Define('prompt_check',R1)\
    .Define('displaceds1_check',R2h1h2)\
    .Define('displaceds2_check',R2h3h4)\
    .Define('displaced_check',R3)


############## (Un)Blinding ########################
#Data is blinded by default, MC is unblinded by default
#MC can be blinded, Data can be unblinded but would need to comment out in if statement below 

if (unblind):
    if isData=='true':
        sys.stderr.write("\nUNBLIND == TRUE BUT LOOKING AT DATA --------- MANUALLY UNBLIND THE PEAK IN DATA - COMMENT OUT ABOVE")
        quit()
    else:
        sys.stderr.write("\nUNBLIND == TRUE AND LOOKING AT MC --------- EVENTS REMAIN AS IS")
        df_higgs_precat = df_higgs_preblind_loose
        df_higgs_precat_test = df_higgs_preblind_loose
else:
    sys.stderr.write("\nUNBLIND == FALSE --------- EVENTS WILL BE BLINDED, EVENTS OMITTED IN mHiggs [122.5,127.5]")
    df_higgs_precat = df_higgs_preblind_loose.Filter('recohiggs_mass_blinded',f'm(mumuhh) blinded in [122.5,127.5]')
    df_higgs_precat_test = df_higgs_preblind_loose.Filter('recohiggs_mass_blinded',f'm(mumuhh) blinded in [122.5,127.5]')

#Look at signal region if unblinded [122.5, 127.5]
if (unblind):
    sys.stderr.write('\nUNBLINDED, LOOKING AT [122.5, 127.5] PEAK REGION')
    df_higgs_precat = df_higgs_precat.Filter('recohiggs_mass_peak','m(mumuhh) in [122.5,127.5] peak')
#Look at sidebands if blinded: Tight = [120, 122.5] U [127.5, 130] || Loose = [110, 122.5] U [127.5, 140]
else:
    #Tight CR condition
    sys.stderr.write('\nBLINDED, LOOKING AT TIGHT SIDEBANDS[120, 122.5] U [127.5, 130]')
    df_higgs_precat = df_higgs_precat.Filter('recohiggs_mass_check','m(mumuhh) in [120,130] sidebands')
    #Loose CR condition
    #sys.stderr.write('\n BLINDED, LOOKING AT LOOSE SIDEBANDS[110, 122.5] U [127.5, 140]')
    #df_higgs_precat = df_higgs_precat.Filter('recohiggs_mass_check_loose','m(mumuhh) in [110,140] sidebands')



############## Histograms Precat ########################

outFile = ROOT.TFile(args.output, "RECREATE")
# outFile.SetBit(ROOT.TFile.k630forwardCompatibility)
outFile.cd()

hists_1d_test={}
hists_2d_test={}
hists_1d_precat={}
hists_2d_precat={}

s12m_range = [0.8, 3.6] if isData else [0.8, 1.6]
s12m_bins = 280 if isData else 80
s12lxy_range = [0, 60]
s12lxy_bins = 300
s12lxysig_range = [0, 1000]
s12lxysig_bins = 50


precat_count = df_higgs_precat.Count().GetValue()
precattest_count = df_higgs_precat_test.Count().GetValue()

sys.stderr.write('\nPRECAT (PEAK SR or TIGHT CR): ' + str(precat_count))
sys.stderr.write('\nPRECAT TEST (LOOSE CR): ' + str(precattest_count))

hists_1d_test['h_'+'ChargedHadron1'+'relIso'] = df_higgs_precat.Histo1D(('h_'+'ChargedHadron1'+'relIso','', 100, 0, 10), 'ch1_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron2'+'relIso'] = df_higgs_precat.Histo1D(('h_'+'ChargedHadron2'+'relIso','', 100, 0, 10), 'ch2_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron3'+'relIso'] = df_higgs_precat.Histo1D(('h_'+'ChargedHadron3'+'relIso','', 100, 0, 10), 'ch3_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron4'+'relIso'] = df_higgs_precat.Histo1D(('h_'+'ChargedHadron4'+'relIso','', 100, 0, 10), 'ch4_reliso', 'weight')
sys.stderr.write('\nhere1')
hists_1d_test['h_'+'ZBoson'+'Mass']=df_recoMu_cut_invcut.Histo1D(('h_'+'ZBoson'+'Mass','', 400, 70, 110), 'Z_mass', 'weight')
hists_1d_test['h_'+'Scalar1'+'Mass'] = df_higgs_precat.Histo1D(('h_'+'Scalar1'+'Mass', '', s12m_bins, s12m_range[0], s12m_range[1]), 's1_mass','weight')
hists_1d_test['h_'+'Scalar2'+'Mass'] = df_higgs_precat.Histo1D(('h_'+'Scalar2'+'Mass', '', s12m_bins, s12m_range[0], s12m_range[1]), 's2_mass','weight')
hists_1d_test['h_'+'Scalar1'+'Lxy'] = df_higgs_precat.Histo1D(('h_'+'Scalar1'+'Lxy', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's1_lxy', 'weight')
hists_1d_test['h_'+'Scalar1'+'LxySig'] = df_higgs_precat.Histo1D(('h_'+'Scalar1'+'LxySig', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's1_lxysign_tmp', 'weight')
hists_1d_test['h_'+'Scalar2'+'Lxy'] = df_higgs_precat.Histo1D(('h_'+'Scalar2'+'Lxy', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's2_lxy', 'weight')
hists_1d_test['h_'+'Scalar2'+'LxySig'] = df_higgs_precat.Histo1D(('h_'+'Scalar2'+'LxySig', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's2_lxysign_tmp', 'weight')

sys.stderr.write('\nhere2')
hists_2d_test['h_'+'Scalar'+'Lxy'] = df_higgs_precat.Histo2D(('h_'+'Scalar'+'Lxy', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1], s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's1_lxy', 's2_lxy', 'weight')
hists_2d_test['h_'+'Scalar'+'LxySig'] = df_higgs_precat.Histo2D(('h_'+'Scalar'+'LxySig', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1], s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's1_lxysign_tmp', 's2_lxysign_tmp', 'weight')


sys.stderr.write('\nhere3')
hists_1d_test['h_'+'ChargedHadron1'+'relIso'+ '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'ChargedHadron1'+'relIso'+ '_Loose','', 100, 0, 10), 'ch1_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron2'+'relIso'+ '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'ChargedHadron2'+'relIso'+ '_Loose','', 100, 0, 10), 'ch2_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron3'+'relIso'+ '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'ChargedHadron3'+'relIso'+ '_Loose','', 100, 0, 10), 'ch3_reliso', 'weight')
hists_1d_test['h_'+'ChargedHadron4'+'relIso'+ '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'ChargedHadron4'+'relIso'+ '_Loose','', 100, 0, 10), 'ch4_reliso', 'weight')

sys.stderr.write('\nhere4')
hists_1d_test['h_'+'HiggsBoson'+'Mass' + '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'HiggsBoson'+'Mass' + '_Loose', '', 300, 110, 140), 'higgs_mass','weight')
hists_1d_test['h_'+'Scalar1'+'Mass' + '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar1'+'Mass' + '_Loose', '', s12m_bins, s12m_range[0], s12m_range[1]), 's1_mass','weight')
hists_1d_test['h_'+'Scalar2'+'Mass' + '_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar2'+'Mass' + '_Loose', '', s12m_bins, s12m_range[0], s12m_range[1]), 's2_mass','weight')
hists_1d_test['h_'+'Scalar1'+'Lxy'+'_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar1'+'Lxy'+'_Loose', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's1_lxy', 'weight')
hists_1d_test['h_'+'Scalar1'+'LxySig'+'_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar1'+'LxySig'+'_Loose', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's1_lxysign_tmp', 'weight')
hists_1d_test['h_'+'Scalar2'+'Lxy'+'_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar2'+'Lxy'+'_Loose', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's2_lxy', 'weight')
hists_1d_test['h_'+'Scalar2'+'LxySig'+'_Loose'] = df_higgs_precat_test.Histo1D(('h_'+'Scalar2'+'LxySig'+'_Loose', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's2_lxysign_tmp', 'weight')

sys.stderr.write('\nhere5')
hists_2d_test['h_'+'Scalar'+'Lxy'+'_Loose'] = df_higgs_precat_test.Histo2D(('h_'+'Scalar'+'Lxy'+'_Loose', '', s12lxy_bins, s12lxy_range[0], s12lxy_range[1], s12lxy_bins, s12lxy_range[0], s12lxy_range[1]), 's1_lxy', 's2_lxy', 'weight')
hists_2d_test['h_'+'Scalar'+'LxySig'+'_Loose'] = df_higgs_precat_test.Histo2D(('h_'+'Scalar'+'LxySig'+'_Loose', '', s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1], s12lxysig_bins, s12lxysig_range[0], s12lxysig_range[1]), 's1_lxysign_tmp', 's2_lxysign_tmp', 'weight')



#Tight
if precat_count != 0:
    sys.stderr.write('\nhere6')
    hists_1d_test['h_ChargedHadron1relIso'].Write()
    hists_1d_test['h_ChargedHadron2relIso'].Write()
    hists_1d_test['h_ChargedHadron3relIso'].Write()
    hists_1d_test['h_ChargedHadron4relIso'].Write()
    hists_1d_test['h_ZBosonMass'].Write()
    hists_1d_test['h_Scalar1Mass'].Write()
    hists_1d_test['h_Scalar2Mass'].Write()
    hists_1d_test['h_Scalar1Lxy'].Write()
    hists_1d_test['h_Scalar2Lxy'].Write()
    hists_1d_test['h_Scalar1LxySig'].Write()
    hists_1d_test['h_Scalar2LxySig'].Write()
    #2D
    hists_2d_test['h_ScalarLxy'].Write()
    hists_2d_test['h_ScalarLxySig'].Write()

#Loose
if precattest_count != 0:
    sys.stderr.write('\nhere7')
    hists_1d_test['h_ChargedHadron1relIso_Loose'].Write()
    hists_1d_test['h_ChargedHadron2relIso_Loose'].Write()
    hists_1d_test['h_ChargedHadron3relIso_Loose'].Write()
    hists_1d_test['h_ChargedHadron4relIso_Loose'].Write()
    hists_1d_test['h_HiggsBosonMass_Loose'].Write()
    hists_1d_test['h_Scalar1Mass_Loose'].Write()
    hists_1d_test['h_Scalar2Mass_Loose'].Write()
    hists_1d_test['h_Scalar1Lxy_Loose'].Write()
    hists_1d_test['h_Scalar2Lxy_Loose'].Write()
    hists_1d_test['h_Scalar1LxySig_Loose'].Write()
    hists_1d_test['h_Scalar2LxySig_Loose'].Write()
    #2D
    hists_2d_test['h_ScalarLxy_Loose'].Write()
    hists_2d_test['h_ScalarLxySig_Loose'].Write()
############## Categorization ########################

reg_={}
#reg_['prompt'] = higgs_definitions_blinded_loose_iso.Filter('prompt_check','only prompt region')
#reg_['displaceds1'] = higgs_definitions_blinded_loose_iso.Filter('displaceds1_check','only displaced s1 region')
#reg_['displaceds2'] = higgs_definitions_blinded_loose_iso.Filter('displaceds2_check','only displaced s2 region')
#reg_['displaced'] = higgs_definitions_blinded_loose_iso.Filter('displaced_check','both displaced region')


############## Histograms Postcat ########################

hists_1d_={}
hists_2d_={}

outFile.Close()

sys.stderr.write("\n\n\nTime taken: --- %s seconds ---" % (time.time() - start_time)+'\n')
sys.stderr.flush()








