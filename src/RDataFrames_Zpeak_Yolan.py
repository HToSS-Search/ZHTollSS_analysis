# In this .py file I will try to write my own version of the Zpeak analysis provided in RDataFrames_Zpeak_alt

import ROOT
import os, argparse
import numpy as np
import sys
import re
import time
import yaml

cpu_count = 16 # give the same for request_cpus on condor script
ROOT.ROOT.EnableImplicitMT(cpu_count)

#Setting up the loading in of the appropriate dataset along with the right cuts, weights, etc...
parser = argparse.ArgumentParser(description='Analysis of Zpeak using Z to muons')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file to process", type=str)
parser.add_argument("-o","--output", dest="output", help="Destination directory", type=str)
# For now the cuts will be defined in this .py file 
# parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)

args = parser.parse_args()
print("Here")
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
print("Here")
#Lets perform a check
weightPlot = dataFile.Get("jmeanalyzer/h_Counter").Clone()
if sumWeights != weightPlot.GetBinContent(1):
	sumWeights = weightPlot.GetBinContent(1)
	print("Sum of weights in config file did not match actual sum of weights, change the value in the config file accordingly")
	print("correct sum of weights is: {} ".format(sumWeights))
else:
	None
	print("Sum of weights in config matches sum of weights calculated from simulation")
	print("correct sum of weights is: {} ".format(sumWeights))


#Load in the tree from the data in the .root file
treeName = "jmeanalyzer/tree"

df = ROOT.RDataFrame(treeName, fname)
print('Tree loaded in succesfully')
df.Describe().Print()





