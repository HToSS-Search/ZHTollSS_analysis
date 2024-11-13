# In this .py file I will try to write my own version of the Zpeak analysis provided in RDataFrames_Zpeak_alt

import ROOT
import os, argparse
import numpy as np
import sys
import re
import time
import yaml

#########################################################################################
# C++ Helper Functions
#########################################################################################



cpp = R"""

using FourVector = ROOT::Math::PxPyPzMVector;
using RVecFloat = ROOT::RVec<Float_t>;

float getLeading(RVecFloat vec){
	auto idxmax = ROOT::VecOps::ArgMax(vec);
	return vec[idxmax];

}

float getTrailing(RVecFloat vec){
	auto idxmin = ROOT::VecOps::ArgMin(vec);
	return vec[idxmin];
}
auto myInvariantMass(RVecFloat& pt, RVecFloat& eta, RVecFloat& phi, float mass){
        Float_t px1 = pt[0]*cos(phi[0]);Float_t py1 = pt[0]*sin(phi[0]);Float_t pz1 = pt[0]*cos(eta[0]);
	Float_t px2 = pt[1]*cos(phi[1]);Float_t py2 = pt[1]*sin(phi[1]);Float_t pz2 = pt[1]*cos(eta[1]);
	FourVector P1 {px1, py1, pz1, mass};
	FourVector P2 {px2, py2, pz2, mass};
	FourVector dilep = P1 + P2;
	auto mass_ = dilep.M();
	return mass_;
}

RVecFloat getDeltaR(RVecFloat& eta, RVecFloat& phi){
/** Get all DeltaR from all possible pairs in the data
 * Input should be the relevant eta and phi columns respectively
 * Returns Vec with length equal to amount of combinations where each entry is the DeltaR for a pair
 */
	if (eta.size() != phi.size()) {
		std::cout << "Eta and phi are not the same size!" << std::endl;
		return RVecFloat();
	} else {
		auto idx = ROOT::VecOps::Combinations(eta, 2);
		RVecFloat deltaR(idx[0].size()); // Initialize deltaR with the correct size
					           
		auto phi1 = Take(phi, idx[0]); auto phi2 = Take(phi, idx[1]);
		auto eta1 = Take(eta, idx[0]); auto eta2 = Take(eta, idx[1]);

		deltaR = ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
		return deltaR;
	}

}
"""


#########################################################################################
# Setting up the environment for loading in the data from configs
#########################################################################################

ROOT.gInterpreter.Declare(cpp)

cpu_count = 16 # give the same for request_cpus on condor script
#ROOT.ROOT.EnableImplicitMT(cpu_count)

#Setting up the loading in of the appropriate dataset along with the right cuts, weights, etc...
parser = argparse.ArgumentParser(description='Analysis of Zpeak using Z to muons')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
# For now the cuts will be defined in this .py file 
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)

args = parser.parse_args()
#Open the config file and set a variable to be dictionary corresponding to content of the config file
fconfig = open(args.config, 'r')
conf_pars = yaml.safe_load(fconfig)

#Get data location from the config file (in my config file location is directly to .root file)
fname = conf_pars['locations'][0]
dataFile = ROOT.TFile(fname)

#This is the  luminosity for the total 2017UL run (see file name).
#These values can be found in the config file and are idealy taken from here in an automated way dependent on which sampleset is called in the commandline to analyze.
givenLuminosity = 41474

#Get the weights, cross section, and luminosity from MonteCarlo. Set to 1 if Data is not MC (Non_MC will always contain 'Run' in name?)
sumWeights = 1 if 'Run' in args.config else conf_pars['sum_weights']
crossSection = 1 if 'Run' in args.config else conf_pars['cross_section']
luminosity = 1 if 'Run' in args.config else givenLuminosity
#Lets perform a check
weightPlot = dataFile.Get("jmeanalyzer/h_Counter").Clone()
if sumWeights != weightPlot.GetBinContent(1):
	sumWeights = weightPlot.GetBinContent(1)
	print("\nSum of weights in config file did not match actual sum of weights, change the value in the config file accordingly")
	print("\ncorrect sum of weights is: {} ".format(sumWeights))
else:
	None
	print("\nSum of weights in config matches sum of weights calculated from simulation")
	print("\ncorrect sum of weights is: {} ".format(sumWeights))


#########################################################################################
# Helper Functions
#########################################################################################


def genTurnOn(df, triggerName, mmin, mmax, steps):
	#Filter out the muons from the dataset before applying trigger to get a better gauge on trigger efficiency
	#If this is not done, max efficiency plateau at 50% due to presence of electrons
	df_muons = df.Filter('_nEles == 0')
	df_muontrigger = df_muons.Filter(triggerName)
	
	pTBins = np.arange(mmin,mmax,steps)

	eff_hist = df_muontrigger.Histo1D(('triggered_muons', 'Triggered muons', len(pTBins)-1, pTBins), '_lPt')
	eff_hist_total = df_muons.Histo1D(('total_muons', 'Total  muons', len(pTBins)-1, pTBins), '_lPt')

	canvas = ROOT.TCanvas('canvas_turnon_' + triggerName, "Turn On Curve " + triggerName, 800, 600)
	canvas.SetGrid()

	turnon = ROOT.TEfficiency(eff_hist.GetValue(), eff_hist_total.GetValue())
	turnon.SetTitle("Turn on curve " + triggerName + ";Muon pT (GeV);Efficiency")
	turnon.Draw("AP")
	
	#Print the maximum reached efficiency
	max_efficiency = 0
	for i in range(1, turnon.GetTotalHistogram().GetNbinsX() + 1):
	    if turnon.GetEfficiency(i) > max_efficiency:
	            max_efficiency = turnon.GetEfficiency(i)
	print("HERE")
	print(triggerName + " Maximum Efficiency:" + str(max_efficiency))
	
	#Generate a vertical line at the trigger threshold
	xpos = 0
	if triggerName[-2:].isdigit():
		xpos = int(triggerName[-2:])
	vLine = ROOT.TLine(xpos, 0, xpos, canvas.GetUymax())
	vLine.SetLineStyle(2)
	vLine.SetLineWidth(2)
	vLine.SetLineColor(ROOT.kBlue)
	vLine.Draw()

	canvas.Write('canvas_turnon_' + triggerName)

def genPtPlot(df, mmin, mmax, bins):

	df_leading_subleading = df.Define('leadingPt', 'getLeading(mu_pt)')\
			.Define('subleadingPt', 'getTrailing(mu_pt)')
	
	# Create a canvas to plot the histograms
	canvas = ROOT.TCanvas("canvas_pt", "Muon PT Distributions", 800, 600)

	# Histograms
	mu_leadingPt_hist = df_leading_subleading.Histo1D(('mu_leadingpt_hist', 'Leading Muon PTs', bins, mmin, mmax), 'leadingPt')
	mu_subleadingPt_hist = df_leading_subleading.Histo1D(('mu_subleadingpt_hist', 'Subleading Muon PT', bins, mmin, mmax), 'subleadingPt')

	mu_leadingPt_hist.GetXaxis().SetTitle("Muon pT (GeV)")  # X-axis label
	mu_leadingPt_hist.GetYaxis().SetTitle("Events")         # Y-axis label

	mu_leadingPt_hist.SetLineColor(ROOT.kBlue)  # Set color for the leading pt histogram
	mu_subleadingPt_hist.SetLineColor(ROOT.kRed)  # Set color for the subleading pt histogram

	mu_leadingPt_hist.Draw()
	mu_subleadingPt_hist.Draw("same")  # "same" will overlay it on the existing plot

	canvas.Write('canvas_pt')

def genXPlot(df, varName, mmin, mmax, bins):
	canvas = ROOT.TCanvas('canvas_' + varName, varName + " Distribution", 800, 600)

	hist = df.Histo1D(('hist_' + varName, varName, bins, mmin, mmax), varName)
	hist.GetXaxis().SetTitle(varName)
	hist.GetYaxis().SetTitle("Events")

	hist.Draw()
	
	canvas.Write('canvas_' + varName)

def genInvMassPlot(df, mmin, mmax, bins):
	canvas = ROOT.TCanvas('canvas_invmass', 'Invariant Z Mass', 800, 600)

	hist = df.Histo1D(('hist_invmass', 'Invariant mass distribution of Z', bins, mmin, mmax), 'dimuon_mass', 'evt_weight')
	hist.GetXaxis().SetTitle('m#_{\mu\mu} (GeV)')
	hist.GetYaxis().SetTitle('Events')

	hist.Draw()

	canvas.Write('canvas_invmass')

#########################################################################################
# Loading in the Tree and performing selection on Dataframe
#########################################################################################

#Load in the tree from the data in the .root file
treeName = "jmeanalyzer/tree"

df = ROOT.RDataFrame(treeName, fname)
print('\nTree loaded in succesfully')

displayList = ['_lEta', '_lPhi', '_lPt', '_nEles', '_nMus', '_lPassTightID', 'HLT_IsoMu27']
df.Display(displayList, 10).Print()

totalEntries = df.Count().GetValue()
print(totalEntries)

#note that _lpt of muons passing the trigger is not necessarily > 27 GeV.
df_muons = df.Filter('_nEles == 0', 'Filter out electrons')
df_muontrigger = df_muons.Filter('HLT_IsoMu27', 'MuonTriggerCut')

#perform a cut
muon_cuts = 'abs(_lEta) < 2.4 && _lPt > 10 && _lPassTightID'
df_muons_aftercut = df_muontrigger.Define('mu_pt', f'_lPt[{muon_cuts}]')\
			.Define('mu_eta', f'_lEta[{muon_cuts}]')\
			.Define('mu_phi', f'_lPhi[{muon_cuts}]')\
			.Define('mu_ch',f'_lpdgId[{muon_cuts}]/13')

#perform a selection
muon_sel = "mu_pt.size() == 2 && mu_ch[0]*mu_ch[1] < 0 && Max(mu_pt) > 30"
df_muons_aftercutselection = df_muons_aftercut.Filter(muon_sel)

#display the new columns made from cut and selection
displayList2 = ['mu_eta', 'mu_phi', 'mu_ch', 'mu_pt', '_lPt', 'HLT_IsoMu27']
df_muons_aftercutselection.Describe().Print()
df_muons_aftercutselection.Display(displayList2).Print()

#get the dimuon mass and put into a df
llim_zmass = 70
ulim_zmass = 130
muon_massVal = 0.105 #GeV

dimuon_masscut = 'dimuon_mass > ' + str(llim_zmass) + ' &&  dimuon_mass < '+ str(ulim_zmass)

df_dimuon = df_muons_aftercutselection.Define('mu_mass', str(muon_massVal))\
		.Define('dimuon_mass', 'myInvariantMass(mu_pt, mu_eta, mu_phi, mu_mass)')\
		.Filter(dimuon_masscut)
df_dimuon = df_dimuon.Define('deltaR', 'getDeltaR(mu_eta, mu_phi)')

df_dimuon.Display(['mu_pt', 'dimuon_mass', 'deltaR']).Print()

#Setup the weights
df_dimuon = df_dimuon.Define("evt_weight", f'({crossSection}*{luminosity}/{sumWeights})*_weight')


if False:
	#Print out the efficiencies for each step:
	totalMuEntries = df_muons.Count().GetValue()
	eff_triggercut = df_muontrigger.Count().GetValue()/totalMuEntries
	eff_cuts = df_muons_aftercut.Count().GetValue()/totalMuEntries
	eff_sel = df_muons_aftercutselection.Count().GetValue()/totalMuEntries
	eff_dimuon = df_dimuon.Count().GetValue()/totalMuEntries
	print("\nTotal events with only muons in dataset: " + str(totalMuEntries))
	print("\nEfficiency after HLTIsoMu27: " + str(eff_triggercut))
	print("\nEfficiency after cuts: " + str(eff_cuts))
	print("\nEfficiency after cuts and selection: " + str(eff_sel))
	print("\nEfficiency after dimuon mass selection: " + str(eff_dimuon))

if (False):
	c1 = ROOT.TCanvas('c1', 'mu_pT before HLT/Cuts/Selections', 800, 600)
	h1 = df_muons.Histo1D(('hist_ptmuons', 'mu_pT before HLT/Cuts/Selections', 320, 0, 80), '_lPt')
	h1.Draw()

	c2 = ROOT.TCanvas('c2', 'mu_eta before HLT/Cuts/Selectiosn', 800, 600)
	h2 = df_muons.Histo1D(('hist_etamuons', 'mu_eta before HLT/Cuts/Selections', 320, -4,4), '_lEta')
	h2.Draw()

	c3 = ROOT.TCanvas('c3', 'mu_phi before HLT/Cuts/Selectiosn', 800, 600)
	h3 = df_muons.Histo1D(('hist_phimuons', 'mu_phi before HLT/Cuts/Selections', 320, -4, 4), '_lPhi')
	h3.Draw()

	c1.Draw()
	c1.Update()
	c2.Draw()
	c2.Update()
	c3.Draw()
	c3.Update()

#print(df.Filter('HLT_IsoTkMu24').Count().GetValue())
#print(df_muons.Filter('HLT_IsoTkMu24').Count().GetValue())


needGeneratePlots = False
if(needGeneratePlots):

	#generate the plots and save to output file
	outFile = ROOT.TFile(args.output, "RECREATE")

	genTurnOn(df, 'HLT_IsoMu27', 0, 80, 0.1)
	genTurnOn(df, 'HLT_IsoMu24', 0, 80, 0.1) 
	genTurnOn(df, 'HLT_IsoTkMu24', 0, 80, 0.1)
	genPtPlot(df_muons_aftercutselection, 0, 80, 320)
	genXPlot(df_muons_aftercutselection, 'mu_eta', -4, 4, 320)
	genXPlot(df_muons_aftercutselection, 'mu_phi', -4, 4, 320)
	genXPlot(df_dimuon, 'deltaR', 0, 5, 320)
	genInvMassPlot(df_dimuon, llim_zmass, ulim_zmass, 320)

	outFile.Close()

print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")













