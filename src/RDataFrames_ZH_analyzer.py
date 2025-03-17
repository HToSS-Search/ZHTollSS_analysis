#This is the RDF Analyzer for the ZH analysis
import numpy as np
from array import array
import os,argparse
import re
import yaml
import ROOT
import time
import sys

from ZH_analysis_helper import gInterpreter_lv, gInterpreter_getIndices, gInterpreter_getKinematics, gInterpreter_pairselection, gInterpreter_matching
from sf_helper import get_pileup, gInterpreter_SF

ROOT.gROOT.SetBatch(True)

cpu_count = 4 # give the same for request_cpus on condor script
# ROOT.ROOT.EnableImplicitMT(cpu_count) #comment it when testing

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
lumi = 1 if 'Run' in args.config else lumi_factor
# lumi = 1 if 'Run' in args.config else 4247.682053046 #2017D for now
isData = 'true' if 'Run' in args.config else 'false'



###################### LOADING ALL FILES FOR PROCESSING ############################
isOldNtuple = False
if 'almorton' in data_loc:
	isOldNtuple = True
# data_loc = data_loc[0]
if data_loc[-1] != '/':
	data_loc = data_loc+'/'
print(data_loc)
if isOldNtuple:
	directories = [data_loc] #below only for old ntuples
else:
	directories = [data_loc+d+"/" for d in os.listdir(data_loc) if os.path.isdir(os.path.join(data_loc, d))]

print(directories)
print(directories2)
# print(directories+directories2)
directories=directories+directories2

treeName = "makeTopologyNtupleMiniAOD/tree"
list_of_files = []
for dirtmp in directories:
	print(dirtmp)
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
			# 	continue
		except Exception as e:
			print(f"Error processing file {fistr}: {e}")
			continue
		list_of_files.append(fistr)

print(list_of_files)



###################### CALCULATION OF SUM OF WEIGHTS ####################
if args.onlyweights:
	# sum_wts calculated here
	print("enters sum weights calculation")
	print(list_of_files)
	if not 'Run' in args.config:
		if not isinstance(list_of_files, list):
			file = ROOT.TFile(list_of_files)
			weightPlot = file.Get("makeTopologyNtupleMiniAOD/weightHisto").Clone()
			weightPlot.SetDirectory(0)
			file.Close()
		else:
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
					print(f"Error processing file {fistr}: {e}")
					continue  # Continue to the next file in case of an error
		totalEvents_ = weightPlot.GetBinContent(2) - weightPlot.GetBinContent(3) # bins filled from 1, but bins available from 0
		sum_wts = totalEvents_
	sys.stderr.write("sum of weights:"+str(sum_wts)+"\n")
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
			print(f"Error processing file {fistr}: {e}")
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

#ADD TRIGGER CONDITIONS HERE WHEN AVAILABLE
#----

#UPDATE WITH RECOSCALAR ONCE IMPLEMENTED
mu_cuts = cut_pars['muons']
Z_cuts = cut_pars['recoZ']
ch_cuts = cut_pars['hadrons']
dihadron_cuts = cut_pars['dihadron']
recohiggs_cuts = cut_pars['recohiggs']
#recoscalar_cuts =cut_pars['recoscalar']



###################### THE RDF ANALYZER ############################
gInterpreter_lv()
gInterpreter_getIndices()
gInterpreter_getKinematics() 
gInterpreter_pairselection()
gInterpreter_matching()
get_pileup()
gInterpreter_SF()

rdf = ROOT.RDataFrame(treeName, list_of_files)
print('\nTree loaded in succesfully')

totalEntries = df.Count().GetValue()
print(f"\nTotal Entries : {totalEntries}")










