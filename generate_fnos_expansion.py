import numpy as np
import os,argparse
import re
import yaml

def maxfilenumber(path):
	list_of_files = os.listdir(path)
	n = [int(re.findall("\d+", str(file))[0]) for file in list_of_files]
	return max(n)

parser = argparse.ArgumentParser(description='Plot stacked histogram')
parser.add_argument("-c", "--config", dest="config",   help="Enter config file (with all processes)", type=str)
parser.add_argument("-a", "--cuts", dest="cuts",   help="Enter cuts file to process", type=str)
parser.add_argument("-n","--njobs", dest="njobs", help="No. of files to process in one job", type=int)
parser.add_argument("-o","--output", dest="outdir", help="Destination directory", type=str)
parser.add_argument("-y","--year", dest="year", help="Year for processing", type=str) #UL2017
parser.add_argument("-f","--fname", dest="fname", help="Name of params file", type=str)
parser.add_argument("-d","--data",dest="isData", help="Data or not; false by default", action="store_true")
# parser.add_argument("-s","--skim",dest="justSkim", help="Skimming or not; false by default", action="store_true")
# parser.add_argument("-w","--weight",dest="justWeights", help="Generate only weight dist? true or false; false by default", action="store_true")
parser.add_argument("-j","--justWeights",dest="justWeights", help="Measure just weights; false by default", action="store_true")
# parser.add_argument("--withSkim",dest="withSkim", help="Use skimmed files instead of bare ntuples; false by default", action="store_true")

args = parser.parse_args()

cwd = os.getcwd()
fin = open(args.config,'r')
f = open(args.fname,'w')
pars = yaml.safe_load(fin)
for conf in pars['datasets']:
	# print(conf.split("/")[-2])
	#print(conf)
	# quit()
	conf_og=conf
	if conf.count('ctauS') == 2:
		# print(new_lt)
		if 'pion' in conf:
			new_lt = conf.split('/')[-1].split('_')[-2].replace('.yaml','')
			conf=conf.replace('_'+new_lt+'_pion.','_pion.')
		else:
			new_lt = conf.split('/')[-1].split('_')[-1].replace('.yaml','')
			conf=conf.replace('_'+new_lt+'.','.')
		print('ENTERS')
		# print(conf_og)
		# print(conf)
		# quit()
	if conf.count('ctauS') == 3:
		if 'pion' in conf:
			old_lt2, new_lt = conf.split('/')[-1].split('_')[-3],conf.split('/')[-1].split('_')[-2].replace('.yaml','')
			conf=conf.replace('_'+old_lt2+'_','').replace(new_lt+'_pion.','_pion.')
		else:
			old_lt2, new_lt = conf.split('/')[-1].split('_')[-2],conf.split('/')[-1].split('_')[-1].replace('.yaml','')
			# new_lt=conf.split('/')[-1].split('.')[0].split('_')[-1]
			conf=conf.replace('_'+old_lt2+'_','').replace(new_lt+'.','.')
		# print(conf)
		# quit()
	f_conf = open(conf,'r')
	conf_pars = yaml.safe_load(f_conf)
	data_loc = conf_pars['locations']
	if isinstance(data_loc,list):
		loc=data_loc[0]
		nmax = maxfilenumber(loc)
		out_name = sample_string[len(sample_string)-2] if sample_string[len(sample_string)-1]=='' else sample_string[len(sample_string)-1]
	else:
		loc = data_loc
		directories = [loc+"/"+d+"/" for d in os.listdir(data_loc) if os.path.isdir(os.path.join(data_loc, d))]
		# print(directories)
		nmax = maxfilenumber(directories[-1])
		out_name = conf_og.split("/")[-1].split(".yaml")[0].replace('_pion','')

	# print(out_name)
	if conf_og.count('ctauS') > 1:
		# data_name = conf_pars['name']+'_'+new_lt
		data_name = conf_og.split('/')[-1].replace('.yaml','').replace('_pion','')
		# data_name = conf_pars['name']+'_'+new_lt
	else:
		data_name = conf_pars['name']
	sample_string = loc.split("/")
	# sample_string_2 = conf.split("/")
	# sample_type = sample_string_2[len(sample_string_2)-2] if ".yaml" in sample_string_2[len(sample_string_2)-1] else sample_string_2[len(sample_string_2)-1]
	outdir_extra = conf_og.split("/datasets/")[1].split(".yaml")[0]

	mk_dir = cwd +'/'+ args.outdir +'/'+outdir_extra
	# print(mk_dir)
	if (args.justWeights): # generating only weight distribution
			if conf_og.count('ctauS') <= 1:
				os.makedirs(mk_dir, exist_ok=True)
	else:
		if 'HToSS' not in data_name:
			os.makedirs(mk_dir, exist_ok=True)
		else:
			os.makedirs(cwd+'/'+args.outdir+'/'+'ggH', exist_ok=True)

	start, end = str(0),str(0) #dummy values
	if 'DYJetsToLL_Pt-650ToInf' in data_name or 'DYJetsToLL_Pt-400To650' in data_name:
		njobs = 2
	else:
		njobs = args.njobs
	for i in range(1, nmax+1, njobs):
		start = str(i)
		end = str(min(nmax,i+njobs-1))
		if (args.justWeights): # generating only weight distribution
			if conf_og.count('ctauS') <= 1:
				#print (conf + "," + "dummy" +"," +args.outdir+"/"+outdir_extra+"/output_"+out_name +"_" + start +"-"+ end + ".root" + "," + start + "," + end + "," + data_name + ','+ args.year)
				line = conf + "," + "dummy" +"," +args.outdir+"/"+outdir_extra+"/output_"+out_name +"_" + start +"-"+ end +  ".root" + "," + start + "," + end + "," + data_name + ','+ args.year
			# else:
			# 	print('problem')
		else:
			# $(cfg) $(output) $(flow) $(fhigh) $(dataset) $(year)
			#print (conf + "," + args.cuts + "," +args.outdir+"/"+outdir_extra+"/output_"+out_name +"_" + start +"-"+ end + ".root" + "," + start + "," + end + "," + data_name + ','+ args.year)
			if 'HToSS' in out_name:
				line = conf + "," + args.cuts + "," +args.outdir+"/"+"ggH"+"/output_"+out_name + ".root" + "," + start + "," + end + "," + data_name + ','+ args.year
			else:
				line = conf + "," + args.cuts + "," +args.outdir+"/"+outdir_extra+"/output_"+out_name +"_" + start +"-"+ end +  ".root" + "," + start + "," + end + "," + data_name + ','+ args.year
			#line = conf + "," + args.cuts + "," +args.outdir+"/"+outdir_extra+"/output_"+out_name + ".root" + "," + start + "," + end + "," + data_name + ','+ args.year
		print(line)
		f.write("\n"+line)
f.close()

