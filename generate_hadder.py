import numpy as np
import os,argparse
import re
import yaml
import shutil

def maxfilenumber(path):
	list_of_files = os.listdir(path)
	n = [int(re.findall("\d+", str(file))[0]) for file in list_of_files]
	return max(n)

parser = argparse.ArgumentParser(description='Plot stacked histogram')
parser.add_argument("-d", "--dir", dest="dir",   help="Enter plots folder to hadd elements", type=str)
parser.add_argument("-l", "--lvl", dest="lvl",   help="1 for hadd-ing subprocess root files, 2 for hadding subprocess into process", type=str)
# parser.add_argument("-j","--justWeights",dest="justWeights", help="Measure just weights; false by default", action="store_true")
# parser.add_argument("--withSkim",dest="withSkim", help="Use skimmed files instead of bare ntuples; false by default", action="store_true")

args = parser.parse_args()

cwd = os.getcwd()
params_hadd = open('hadd_params_lv'+args.lvl+'.txt','w')

all_files = os.listdir(args.dir)
# print(os.listdir(args.dir))
all_dirs = []
for f in all_files:
	if '.' not in f:
		all_dirs.append(f)
# all_dirs =[f for f in all_files]
print(all_dirs)
if args.lvl == '2':
	if not os.path.exists(cwd+'/'+args.dir+'/'+'total'):
		os.makedirs(cwd+'/'+args.dir+'/'+'total')
for d in all_dirs:
	hadd_dir=[]
	cp_dir=[]
	print(d)
	if 'total' in d:
		continue
	proc_dir = cwd+'/'+args.dir+'/'+d
	ls_proc_dir = os.listdir(proc_dir)
	if args.lvl == '1':
		for tmp in ls_proc_dir:
			tmpdir = proc_dir+'/'+tmp
			if '.' not in tmpdir:
				if len(os.listdir(tmpdir)) == 1:
					cp_dir.append(tmpdir+'/'+os.listdir(tmpdir)[0])
				elif len(os.listdir(tmpdir)) > 1:
					hadd_dir.append(tmpdir)
				else:
					continue
		print('hadd_dir')
		print(hadd_dir)
		print('cp_dir')
		print(cp_dir)
		for hd in hadd_dir:
			line = proc_dir+'/output_'+hd.split('/')[-1]+'.root'+','+hd+','+hd.split('/')[-1]
			params_hadd.write("\n"+line)
		for cd in cp_dir:
			rt_f = cd.split('/')[-1]
			dir_f = cd.split('/')[-2]
			print('COPYING:')
			print(cd, proc_dir+'/output_'+dir_f+'.root')
			shutil.copy(cd, proc_dir+'/output_'+dir_f+'.root')
	elif args.lvl == '2':
		if 'ggH' not in proc_dir:
			print('look here')
			# print(ls_proc_dir)
			dirs, rootfs = [],[]
			for tmp in ls_proc_dir:
				if 'HToSS' in tmp:
					continue
				if '.root' in tmp:
					rootfs.append(tmp.replace('.root','').replace('output_',''))
				else:
					dirs.append(tmp)
			print(rootfs)
			print(dirs)
			print(set(rootfs)==set(dirs))
			if set(rootfs)!=set(dirs):
				continue
			total_dir = cwd+'/'+args.dir+'/'+'total'+'/output_'+d+'.root'
			line = total_dir+','+proc_dir+','+total_dir.split('/')[-1].replace('.root','').replace('output_','')
			print(line)
			params_hadd.write("\n"+line)
		else:
			for tmp in ls_proc_dir:
				if '.root' in tmp:
					print('COPYING:')
					print(proc_dir+'/'+tmp,cwd+'/'+args.dir+'/'+'total'+'/'+tmp)
					shutil.copy(proc_dir+'/'+tmp,cwd+'/'+args.dir+'/'+'total'+'/'+tmp)
params_hadd.close()
