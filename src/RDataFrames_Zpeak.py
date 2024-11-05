import numpy as np
import os,argparse
import re
import yaml
import ROOT
import time
import sys

cpu_count = 4 # give the same for request_cpus on condor script
ROOT.ROOT.EnableImplicitMT(cpu_count)

def maxfilenumber(path):
	list_of_files = os.listdir(path)
	n = [int(re.findall("\d+", str(file))[0]) for file in list_of_files]
	return min(n),max(n)
getLV_code = '''
ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> makeLVs(ROOT::VecOps::RVec<Int_t>& mu_ch,ROOT::VecOps::RVec<Float_t>& mu_px,ROOT::VecOps::RVec<Float_t>& mu_py,ROOT::VecOps::RVec<Float_t>& mu_pz, double mu_mass) {
    ROOT::VecOps::RVec<ROOT::Math::PxPyPzMVector> lvs;
    ROOT::Math::PxPyPzMVector lVec1 {mu_px[0], mu_py[0], mu_pz[0], mu_mass};
    ROOT::Math::PxPyPzMVector lVec2 {mu_px[1], mu_py[1], mu_pz[1], mu_mass};
    ROOT::Math::PxPyPzMVector dilep = lVec1+lVec2;
    lvs.push_back(lVec1);
    lvs.push_back(lVec2);
    lvs.push_back(dilep);
    return lvs;
}
'''
getKinematics_code ='''
using FourVector = ROOT::Math::PxPyPzMVector;
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
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
parser.add_argument("-t","--totalfromcfg", dest="total", help="Take total from config or not", action='store_true')


args = parser.parse_args()

fin = open(args.config,'r')
conf_pars = yaml.safe_load(fin)

data_loc = conf_pars['locations'][0]
cross_section = 1 if 'Run' in args.config else conf_pars['cross_section']
sum_wts = 1 if 'Run' in args.config else conf_pars['sum_weights']
# lumi = 1 if 'Run' in args.config else 41474 #2017 for now
lumi = 1 if 'Run' in args.config else 4247.682053046 #2017D for now


if data_loc[-1] != '/':
	data_loc = data_loc+'/'
directories = [data_loc+d+"/" for d in os.listdir(data_loc) if os.path.isdir(os.path.join(data_loc, d))]

print(directories)
data_name = conf_pars['name']
treeName = "makeTopologyNtupleMiniAOD/tree"
list_of_files = []
for dir in directories:
	print(dir)
	flow, fhigh = maxfilenumber(dir)
	for i in range(flow, fhigh+1):
		fno = str(i)
		fistr = dir+"output_"+fno+".root"
		if not os.path.exists(fistr):
			continue
		list_of_files.append(fistr)
if not args.total: # calculate sum of weights
	print("enters sum weights calculation")
	if not 'Run' in args.config:
		file = ROOT.TFile(list_of_files[0])
		weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
		weightPlot.SetDirectory(0)
		file.Close()
		# print(weightPlot)
		# total_entries = 0
		for i,fistr in enumerate(list_of_files):
			# print(fistr)
			# if i == 10:
			# 	break
			if i==0:
				continue
			if not os.path.exists(fistr):
				continue
			# print("Processing: ",fistr)
			try:
				# Open the ROOT file
				file = ROOT.TFile.Open(fistr)

				# Check if the file was opened successfully
				if not file or file.IsZombie() or file.TestBit(ROOT.TFile.kRecovered):
					raise Exception(f"Error opening file: {fistr}")

				# Process the file
				tmpPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
				# total_entries=total_entries+tmpPlot.GetEntries()
				weightPlot.Add(tmpPlot)
				# Add your code here to work with the file

				# Close the file
				file.Close()
				# continue

			except Exception as e:
				print(f"Error processing file {fistr}: {e}")
				continue  # Continue to the next file in case of an error
		totalEvents_ = weightPlot.GetBinContent(1)
		sum_wts = totalEvents_
sys.stderr.write("sum of weights:"+str(sum_wts)+"\n")

rdf = ROOT.RDataFrame(treeName,list_of_files)
print('RDF loaded')
entries_total = rdf.Count()
sys.stderr.write("entries total:"+str(entries_total.GetValue())+'\n')
sys.stderr.flush()

trigger_cuts = 'HLT_IsoMu27_v > 0' # for the new ntuples
muon_cuts = 'abs(muonPF2PATEta) <2.4 && muonPF2PATTightCutId && muonPF2PATPt >= 20'
iso_cuts = 'muonPF2PATComRelIsodBeta < 0.15'
leadingmu_cut = 'Max(mu_pt)>20' #ordered collection
dilepton_cut = 'dilep_mass > 70 && dilep_mass < 110'

trigger_definitions = rdf.Define('mu_trig',trigger_cuts)
if not '2017B' in args.config:
	cut_trig = trigger_definitions.Filter('mu_trig') 
else:
	cut_trig = trigger_definitions # there was a filter problem here
cut_met = cut_trig
muon_iso_cuts = muon_cuts + ' && ' + iso_cuts
muon_definitions = cut_met.Define('mu_pt',f'muonPF2PATPt[{muon_iso_cuts}]')\
				.Define('mu_ch',f'muonPF2PATCharge[{muon_iso_cuts}]')\
				.Define('mu_px',f'muonPF2PATPX[{muon_iso_cuts}]')\
				.Define('mu_py',f'muonPF2PATPY[{muon_iso_cuts}]')\
				.Define('mu_pz',f'muonPF2PATPZ[{muon_iso_cuts}]')\
				.Define('mu_sel','mu_pt.size()==2' + ' && ' + leadingmu_cut + ' && mu_ch[0]*mu_ch[1] < 0')
cut_mu_sel = muon_definitions.Filter('mu_sel') # selecting events with at least 2 opp. charged muons passing the cuts
# below, using C++ code defined at the start to access ROOT lorentz vectors
dimuon_definitions = cut_mu_sel.Define('LVs','makeLVs(mu_ch,mu_px,mu_py,mu_pz,0.105)')\
				.Define('lv_pt','getKinematics(LVs, "pt")')\
				.Define('lv_eta','getKinematics(LVs, "eta")')\
				.Define('lv_phi','getKinematics(LVs, "phi")')\
				.Define('lv_mass','getKinematics(LVs, "mass")')\
				.Define('dilep_mass','lv_mass[2]')\
				.Define('dilep_sel',dilepton_cut)\
				.Define('evt_wt',f'({cross_section}*{lumi}/{sum_wts})*processMCWeight')
cut_dimuon = dimuon_definitions.Filter('dilep_sel') # selecting events with dilepton mass cuts

entries_trig = cut_trig.Count()
# entries_met= cut_met.Count()
entries_mu_sel= cut_mu_sel.Count()
entries_dimuon= cut_dimuon.Count()
# entries_test= cut_test.Count()

sys.stderr.write ('total:'+str(entries_total.GetValue())+'\n')
sys.stderr.write('After trigger cut:'+str(entries_trig.GetValue())+'\n')
sys.stderr.write('After muon selection:'+str(entries_mu_sel.GetValue())+'\n')
sys.stderr.write('After dimuon cut:'+str(entries_dimuon.GetValue())+'\n')
sys.stderr.flush()

if not 'Run' in args.config:
	histA = cut_dimuon.Histo1D(("dilep_mass","dilepton mass",160,70,110),"dilep_mass","evt_wt")
else:
	histA = cut_dimuon.Histo1D(("dilep_mass","dilepton mass",160,70,110),"dilep_mass")
# histA.Draw()
# time.sleep(10)
outFile = ROOT.TFile(args.output, "RECREATE")
outFile.cd()
histA.Write()
if not 'Run' in args.config and not args.total:
	weightPlot.Write()

sys.stderr.write("Time taken: --- %s seconds ---" % (time.time() - start_time)+'\n')
sys.stderr.flush()





