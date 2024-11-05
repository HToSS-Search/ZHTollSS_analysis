import numpy as np
import os,argparse
import re
import yaml
import ROOT
import time
import sys

cpu_count = 16 # give the same for request_cpus on condor script
ROOT.ROOT.EnableImplicitMT(cpu_count)

def maxfilenumber(path):
	list_of_files = os.listdir(path)
	n = [int(re.findall("\d+", str(file))[0]) for file in list_of_files]
	return max(n)
getLV_code = '''
ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> makeLVs(ROOT::VecOps::RVec<Int_t>& mu_ch,ROOT::VecOps::RVec<Float_t>& mu_pt,ROOT::VecOps::RVec<Float_t>& mu_eta,ROOT::VecOps::RVec<Float_t>& mu_phi, double mu_mass) {
    ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> lvs;
    ROOT::Math::PtEtaPhiMVector lVec1 {mu_pt[0], mu_eta[0], mu_phi[0], mu_mass};
    ROOT::Math::PtEtaPhiMVector lVec2 {mu_pt[1], mu_eta[1], mu_phi[1], mu_mass};
    ROOT::Math::PtEtaPhiMVector dilep = lVec1+lVec2;
    lvs.push_back(lVec1);
    lvs.push_back(lVec2);
    lvs.push_back(dilep);
    return lvs;
}
'''
getKinematics_code ='''
using FourVector = ROOT::Math::PtEtaPhiMVector;
using namespace ROOT::VecOps;
RVec<Float_t> getKinematics(const RVec<FourVector> &tracks, TString var = "pt")
{
   auto pt = [](const FourVector &v) { return v.Pt(); };
   auto eta = [](const FourVector &v) { return v.Eta(); };
   auto phi = [](const FourVector &v) { return v.Phi(); };
   auto mass = [](const FourVector &v) { return v.M(); };
   if (var.Contains("eta")) return Map(tracks, eta);
   else if (var.Contains("phi")) return Map(tracks, phi);
   else if (var.Contains("mass")) return Map(tracks, mass);
   else return Map(tracks, pt);
}
'''
ROOT.gInterpreter.Declare(getKinematics_code)
ROOT.gInterpreter.Declare(getLV_code)
start_time = time.time()
parser = argparse.ArgumentParser(description='Plot stacked histogram')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("-y", "--year", dest="year",   help="Enter year/data/MC to process", type=str)
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)

args = parser.parse_args()

fin = open(args.config,'r')
conf_pars = yaml.safe_load(fin)

# hardcoding the locations - ideally would take from the yaml config file
if '2017B' in args.year:
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DoubleMuon/Run2017B_UL2017_MiniAODv2_v1_DataUL2017B_ZJetsResiduals_May2023/TOTAL.root"
elif "2017C" in args.year:
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DoubleMuon/Run2017C_UL2017_MiniAODv2_v1_DataUL2017C_ZJetsResiduals_May2023/TOTAL.root"
elif "2017D" in args.year:
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DoubleMuon/Run2017D_UL2017_MiniAODv2_v1_DataUL2017D_ZJetsResiduals_May2023/TOTAL.root"
elif "2017E" in args.year:
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DoubleMuon/Run2017E_UL2017_MiniAODv2_v2_DataUL2017E_ZJetsResiduals_May2023/TOTAL.root"
elif "2017F" in args.year:
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DoubleMuon/Run2017F_UL2017_MiniAODv2_v1_DataUL2017F_ZJetsResiduals_May2023/TOTAL.root"
else:
	# print("enters here")
	fname="/pnfs/iihe/cms/store/user/sdansana/JetMETStudies/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/TOTAL_UL2017.root"

treeName = "jmeanalyzer/tree"
dfile = ROOT.TFile(fname)
weightPlot = dfile.Get("jmeanalyzer/h_Counter").Clone()
sum_wts = 1 if 'Run' in args.config else weightPlot.GetBinContent(1) # Sum of weights taken from weight histogram; ideally taken from config file
cross_section = 1 if 'Run' in args.config else 5398 
lumi = 1 if 'Run' in args.config else 41474 #2017 for now - hardcoded again, should depend on year

# data_loc = conf_pars['locations'][0]

# print("sum of entries from weight plots:",total_entries)
rdf = ROOT.RDataFrame(treeName,fname)
print('RDF loaded')
entries_total = rdf.Count()
# auto myHist2 = myDf.Histo1D<float>({"histName", "histTitle", 64u, 0., 128.}, "myColumn");
sys.stderr.write("entries total:"+str(entries_total.GetValue())+'\n')
sys.stderr.flush()

# NOTE: cuts are hard-coded below (and their corresponding numbers) - ideally create "cuts" yaml file which reads the trigger cut / pT cuts depending on single muon or dimuon

# trigger_cuts = 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 > 0'
trigger_cuts = 'HLT_IsoMu27 > 0'
muon_cuts = 'abs(_lEta) <2.4 && _lPassTightID && _lPt >= 20'
leadingmu_cut = 'Max(mu_pt)>20' #ordered collection
dilepton_cut = 'dilep_mass > 70 && dilep_mass < 110'

# Requiring only muons with pT > 20 GeV might be too high for a subleading muon: Try reducing the pT threshold maybe? 

trigger_definitions = rdf.Define('mu_trig',trigger_cuts)
				# .Define('met_filters', 'Flag_goodVertices > 0 && Flag_globalTightHalo2016Filter > 0 && Flag_HBHENoiseFilter > 0 && Flag_HBHENoiseIsoFilter > 0 && Flag_EcalDeadCellTriggerPrimitiveFilter > 0 && Flag_BadPFMuonFilter > 0 && Flag_ecalBadCalibFilter > 0')
if not '2017B' in args.config:
	cut_trig = trigger_definitions.Filter('mu_trig')
else:
	cut_trig = trigger_definitions # For 2017B accept all events because of some issue with the trigger
cut_met = cut_trig

muon_definitions = cut_met.Define('mu_pt',f'_lPt[{muon_cuts}]')\
				.Define('mu_ch',f'_lpdgId[{muon_cuts}]/13')\
				.Define('mu_eta',f'_lEta[{muon_cuts}]')\
				.Define('mu_phi',f'_lPhi[{muon_cuts}]')\
				.Define('mu_sel','mu_pt.size()==2' + ' && ' + leadingmu_cut + ' && mu_ch[0]*mu_ch[1] < 0')
cut_mu_sel = muon_definitions.Filter('mu_sel')
dimuon_definitions = cut_mu_sel.Define('LVs','makeLVs(mu_ch,mu_pt,mu_eta,mu_phi,0.105)')\
				.Define('lv_pt','getKinematics(LVs, "pt")')\
				.Define('lv_eta','getKinematics(LVs, "eta")')\
				.Define('lv_phi','getKinematics(LVs, "phi")')\
				.Define('lv_mass','getKinematics(LVs, "mass")')\
				.Define('dilep_mass','lv_mass[2]')\
				.Define('dilep_sel',dilepton_cut)
cut_dimuon = dimuon_definitions.Filter('dilep_sel')
if not 'Run' in args.config: #'Run' in args.config is identifier for MC or otherwise
	cut_dimuon = cut_dimuon.Define('sample_wt',f'({cross_section}*{lumi}/{sum_wts})')\
				.Define('evt_wt',f'({cross_section}*{lumi}/{sum_wts})*_weight') 
else: 
	cut_dimuon = cut_dimuon.Define('evt_wt','1')
	# If MC, define event weight 
entries_trig = cut_trig.Count()
entries_met= cut_met.Count()
entries_mu_sel= cut_mu_sel.Count()
entries_dimuon= cut_dimuon.Count()
# entries_test= cut_test.Count()

sys.stderr.write ('total:'+str(entries_total.GetValue())+'\n')
sys.stderr.write('After trigger cut:'+str(entries_trig.GetValue())+'\n')
sys.stderr.write('After MET filters:'+str(entries_met.GetValue())+'\n')
sys.stderr.write('After muon selection:'+str(entries_mu_sel.GetValue())+'\n')
sys.stderr.write('After dimuon cut:'+str(entries_dimuon.GetValue())+'\n')
sys.stderr.flush()


# Define histograms below & write
# NOTE: Add histograms below and then write them to the ROOT file
hists_1d_ = {}

# if not 'Run' in args.config:
hists_1d_["h_dilep_mass"] = cut_dimuon.Histo1D(("h_dilep_mass","dilepton mass",160,70,110),"dilep_mass","evt_wt")
# histA.Draw()
# time.sleep(10)
outFile = ROOT.TFile(args.output, "RECREATE")
outFile.cd()
for hname in hists_1d_:
	hists_1d_[hname].Write()

sys.stderr.write("Time taken: --- %s seconds ---" % (time.time() - start_time)+'\n')
sys.stderr.flush()


