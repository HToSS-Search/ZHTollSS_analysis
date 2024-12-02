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
start_time = time.time()
cpu_count = 16 # give the same for request_cpus on condor script
ROOT.ROOT.EnableImplicitMT(cpu_count)

#Setting up the loading in of the appropriate dataset along with the right cuts, weights, etc...
parser = argparse.ArgumentParser(description='Analysis of Zpeak using Z to muons')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
parser.add_argument("-y", "--year", help="Enter the year of the dataset you want to analyse", type=str)
# For now the cuts will be defined in this .py file 
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)

args = parser.parse_args()
#Open the config file and set a variable to be dictionary corresponding to content of the config file
fconfig = open(args.config, 'r')
conf_pars = yaml.safe_load(fconfig)

#Get the cuts
fcuts = open(conf_pars['cuts'])
cuts_pars = yaml.safe_load(fcuts)

#Get data location from the config file (in my config file location is directly to .root file)
fname = ''
locations_list = conf_pars['locations']
for location in locations_list:
    #assumes location names are ordered: 2017, 2017B, 2017B?, 2017B??.., 2017C, etc.
    if args.year in location:
        fname = location
        break
dataFile = ROOT.TFile(fname)

#This is the  luminosity for the total 2017UL run (see file name).
#These values can be found in the config file and are idealy taken from here in an automated way dependent on which sampleset is called in the commandline to analyze.
#givenLuminosity = conf_pars['luminosity'][args.year]

#Get the weights, cross section, and luminosity from MonteCarlo. Set to 1 if Data is not MC (Non_MC will always contain 'Run' in name?)
sumWeights = 1 if 'Run' in fname else conf_pars['sum_weights']
crossSection = 1 if 'Run' in fname  else conf_pars['cross_section']
luminosity = 1 if 'Run' in fname else 1

#Lets perform a check
if not 'Run' in fname:
    weightPlot = dataFile.Get("jmeanalyzer/h_Counter").Clone()
    if sumWeights != weightPlot.GetBinContent(1):
        sumWeights = weightPlot.GetBinContent(1)
        print("\nSum of weights in config file did not match actual sum of weights, change the value in the config file accordingly")
        print("\ncorrect sum of weights is: {} ".format(sumWeights))
    else:
        print("\nSum of weights in config matches sum of weights calculated from simulation")
        print("\ncorrect sum of weights is: {} ".format(sumWeights))


#########################################################################################
# Helper Functions
#########################################################################################

def genTurnOn(df, dataTag, triggerName, mmin, mmax, steps):
    dataStr = 'Data' + args.year if dataTag else 'MC'
    #Filter out the muons from the dataset before applying trigger to get a better gauge on trigger efficiency
    #If this is not done, max efficiency plateau at 50% due to presence of electrons
    df_muons = df.Filter('_nEles == 0')
    df_muontrigger = df_muons.Filter(triggerName)
    
    #Get the leading muon pT since this is most likely the muon the event trigger triggered on
    df_muons = df_muons.Define('leadingPt', 'getLeading(_lPt)')\
            .Define('subleadingPt', 'getTrailing(_lPt)')
    df_muontrigger = df_muontrigger.Define('leadingPt', 'getLeading(_lPt)')\
            .Define('subleadingPt', 'getTrailing(_lPt)')
    
    pTBins = np.arange(mmin,mmax,steps)

    eff_hist = df_muontrigger.Histo1D(('hist_triggerleading', 'Triggered muons', len(pTBins)-1, pTBins), 'leadingPt')
    eff_hist_total = df_muons.Histo1D(('hist_totalleading', 'Total  muons', len(pTBins)-1, pTBins), 'leadingPt')

    turnon = ROOT.TEfficiency(eff_hist.GetValue(), eff_hist_total.GetValue())
    turnon.SetName("turnon_" + triggerName + '_' + dataStr)
    turnon.SetTitle("Turn on curve " + triggerName + dataStr + ";Muon pT (GeV);Efficiency")
    
    #Print the maximum reached efficiency
    max_efficiency = 0
    for i in range(1, turnon.GetTotalHistogram().GetNbinsX() + 1):
        if turnon.GetEfficiency(i) > max_efficiency:
                max_efficiency = turnon.GetEfficiency(i)
    print(triggerName + " Maximum Efficiency:" + str(max_efficiency))
    
    turnon.Write()

def genPtPlot(df, datatag, mmin, mmax, bins):
    dataStr = 'Data' + args.year if dataTag else 'MC'

    df_leading_subleading = df.Define('leadingPt', 'getLeading(mu_pt)')\
        .Define('subleadingPt', 'getTrailing(mu_pt)')
    # Histograms
    if dataTag:

        mu_leadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_leading_' + dataStr, 'Leading Muon PTs ' + dataStr, bins, mmin, mmax), 'leadingPt')
        mu_subleadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_subleading_' + dataStr, 'Subleading Muon PT ' + dataStr, bins, mmin, mmax), 'subleadingPt')

        mu_leadingPt_hist.GetXaxis().SetTitle("Muon pT (GeV)")  # X-axis label
        mu_leadingPt_hist.GetYaxis().SetTitle("Events")         # Y-axis label
        mu_leadingPt_hist.SetLineColor(ROOT.kBlue)  # Set color for the leading pt histogram
        mu_subleadingPt_hist.SetLineColor(ROOT.kRed)  # Set color for the subleading pt histogram

        mu_leadingPt_hist.Write()
        mu_subleadingPt_hist.Write()
    else:
        for year in years:
            mu_leadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_leading_' + dataStr + '_scaled' + year, 'Leading Muon PTs ' + dataStr, bins, mmin, mmax), 'leadingPt', 'evt_weight_' + year)
            mu_subleadingPt_hist = df_leading_subleading.Histo1D(('hist_mu_pT_subleading_' + dataStr + '_scaled' + year, 'Subleading Muon PT ' + dataStr, bins, mmin, mmax), 'subleadingPt', 'evt_weight_' + year)

            mu_leadingPt_hist.GetXaxis().SetTitle("Muon pT (GeV)")  # X-axis label
            mu_leadingPt_hist.GetYaxis().SetTitle("Events")         # Y-axis label
            mu_leadingPt_hist.SetLineColor(ROOT.kBlue)  # Set color for the leading pt histogram
            mu_subleadingPt_hist.SetLineColor(ROOT.kRed)  # Set color for the subleading pt histogram

            mu_leadingPt_hist.Write()
            mu_subleadingPt_hist.Write()


def genXPlot(df, dataTag, varName, mmin, mmax, bins):
    dataStr = 'Data'+ args.year if dataTag else 'MC'
    if dataTag:
        hist = df.Histo1D(('hist_' + varName + '_' + dataStr, varName + dataStr , bins, mmin, mmax), varName)
        hist.GetXaxis().SetTitle(varName)
        hist.GetYaxis().SetTitle("Events")
        hist.Write()
    else:
        for year in years:
            hist = df.Histo1D(('hist_' + varName + '_' + dataStr + '_scaled' + year, varName + dataStr , bins, mmin, mmax), varName, 'evt_weight_' + year)
            hist.GetXaxis().SetTitle(varName)
            hist.GetYaxis().SetTitle("Events")
            hist.Write()


def genInvMassPlot(df, dataTag, mmin, mmax, bins):
    dataStr = 'Data' + args.year if dataTag else 'MC'
    if dataTag:
        hist = df.Histo1D(('hist_mu_invmass_' + dataStr, 'Invariant mass distribution of Z', bins, mmin, mmax), 'dimuon_mass')
        hist.GetXaxis().SetTitle('m#_{\mu\mu} (GeV)')
        hist.GetYaxis().SetTitle('Events')
        hist.Write()

    else:
        for year in years:
            hist = df.Histo1D(('hist_mu_invmass_' + dataStr + '_scaled' + year, 'Invariant mass distribution of Z', bins, mmin, mmax), 'dimuon_mass', 'evt_weight_' + year)
            hist.GetXaxis().SetTitle('m#_{\mu\mu} (GeV)')
            hist.GetYaxis().SetTitle('Events')
            hist.Write()


#########################################################################################
# Loading in the Tree and performing selection on Dataframe
#########################################################################################

#Load in the tree from the data in the .root file
treeName = "jmeanalyzer/tree"

df = ROOT.RDataFrame(treeName, fname)
print('\nTree loaded in succesfully')

totalEntries = df.Count().GetValue()
print(f"\nTotal Entries : {totalEntries}")

#note that _lpt of muons passing the trigger is not necessarily > 27 GeV.
df_muons = df.Filter('_nEles == 0', 'Filter out electrons')
df_muontrigger = df_muons.Filter('HLT_IsoMu27', 'MuonTriggerCut')
df_muontrigger2 = df_muons.Filter('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL')
df_muontrigger3 = df_muons.Filter('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')

#SKIM INCLUDES:
#pT > 20, _lPassTightId, events only 2 leptons, m_ll < 110 GeV

#Perform a cut
#eta < 2.4
#_lPassTightID
muon_cuts = cuts_pars['muon']
muCut1 = 'abs(_lEta) < ' + str(muon_cuts['eta']) + ' && _lPassTightID'
df_muons_aftercut = df_muontrigger.Define('mu_pt', f'_lPt[{muCut1}]')\
            .Define('mu_eta', f'_lEta[{muCut1}]')\
            .Define('mu_phi', f'_lPhi[{muCut1}]')\
            .Define('mu_ch',f'_lpdgId[{muCut1}]/13')

#perform a selection
#Leading pT > 30 GeV (turnon at 27), oppositely charged, 2 muons
muon_sel = "mu_pt.size() == 2 && mu_ch[0]*mu_ch[1] < 0 && Max(mu_pt) > " + str(muon_cuts['leadingPt'])
df_muons_aftercutselection = df_muons_aftercut.Filter(muon_sel)

#get the dimuon mass and put into a df, mass window need be symmetric around ~90GeV
llim_zmass = muon_cuts['massLower']
ulim_zmass = muon_cuts['massUpper']
muon_massVal = 0.105 #GeV

dimuon_masscut = 'dimuon_mass > ' + str(llim_zmass) + ' &&  dimuon_mass < '+ str(ulim_zmass)

df_dimuon = df_muons_aftercutselection.Define('mu_mass', str(muon_massVal))\
        .Define('dimuon_mass', 'myInvariantMass(mu_pt, mu_eta, mu_phi, mu_mass)')\
        .Filter(dimuon_masscut)
df_dimuon = df_dimuon.Define('mu_deltaR', 'getDeltaR(mu_eta, mu_phi)')


#Setup the weights
years = []
if not 'Run' in fname:
    luminosity_ = conf_pars.get("luminosity", {})
    for year, lumi in luminosity_.items():
        year = str(year)
        df_dimuon = df_dimuon.Define("evt_weight_" + year, f'({crossSection}*{lumi}/{sumWeights})*_weight')
        years.append(year)


#########################################################################################
# Extracting information and distributions
#########################################################################################

if True:
    #Print out the efficiencies for each step:
    totalMuEntries = df_muons.Count().GetValue()
    eff_triggercut = df_muontrigger.Count().GetValue()/totalMuEntries
    eff_triggercut2 = df_muontrigger2.Count().GetValue()/totalMuEntries
    eff_triggercut3 = df_muontrigger3.Count().GetValue()/totalMuEntries
    eff_cuts = df_muons_aftercut.Count().GetValue()/totalMuEntries
    eff_sel = df_muons_aftercutselection.Count().GetValue()/totalMuEntries
    eff_dimuon = df_dimuon.Count().GetValue()/totalMuEntries
    print("\nTotal events with only muons in dataset: " + str(totalMuEntries))
    print("\nEfficiency after HLTIsoMu27: " + str(eff_triggercut))
    print("\nEfficiency after HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL: " + str(eff_triggercut2))
    print("\nEfficiency after HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ: " + str(eff_triggercut3))
    print("\nEfficiency after cuts: " + str(eff_cuts))
    print("\nEfficiency after cuts and selection: " + str(eff_sel))
    print("\nEfficiency after dimuon mass selection: " + str(eff_dimuon))


#Generate Histograms and put them into a Dictonary
dataTag = True if 'Run' in fname else False

needGeneratePlots = True
if(needGeneratePlots):
    #generate the plots and save to output file
    outFile = ROOT.TFile(args.output, "UPDATE")
    print("\n Generating TurnOn curves...")
    genTurnOn(df, dataTag, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 0, 80, 0.1)
    genTurnOn(df, dataTag, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 0, 80, 0.1)
    genTurnOn(df, dataTag, 'HLT_IsoMu27', 0, 80, 0.1) 
    genTurnOn(df, dataTag, 'HLT_IsoMu24', 0, 80, 0.1)
    print("\n Generating Pt Plot...")
    genPtPlot(df_dimuon, dataTag, 0, 80, 160)
    print("\n Generating Eta Plot...")
    genXPlot(df_dimuon, dataTag, 'mu_eta', -4, 4, 160)
    print("\n Generating Phi Plot...")
    genXPlot(df_dimuon, dataTag, 'mu_phi', -4, 4, 160)
    print("\n Generating DeltaR Plot...")
    genXPlot(df_dimuon, dataTag, 'mu_deltaR', 0, 5, 160)
    print("\n Generating InvMass Plot...")
    genInvMassPlot(df_dimuon, dataTag, llim_zmass, ulim_zmass, 160)

    outFile.Close()

sys.stderr.write("\nTime taken: --- %s seconds ---" % (time.time() - start_time))

print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")













