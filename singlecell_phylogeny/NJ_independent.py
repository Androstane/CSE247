import numpy as np
import os 
import sys
from os import walk 
import argparse
import ete3
from ete3 import Tree
import matplotlib.pyplot as plt
import seaborn as sns 
import pandas as pd
import dendropy as dnd
sns.set_style("darkgrid")
import parse

class UnAcceptedCSVLine(Exception):
	def __init__(self, data):
		self.data = data
	def __str__(self):
		return repr(self.data)

def Nexus_state(the_string):
	#### The input must be an integer 
	#### The output is the character suitable for PAUP
	cut = 36
	alphabets = ["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
	state = None
	if int(the_string)<=cut-1:
		state = alphabets[int(the_string)]
	else:
		state = alphabets[-1]
	return state

def stateExtract(label):
	stateString = ""
	for i in label:
		if i.isdigit():
			stateString += i
	return int(stateString)

def NJ(bin_array, samples, output_address):
	output_f = open(output_address,"w")
	output_f.write("#NEXUS"+'\n')
	output_f.write("begin taxa;"+"\n"+"\t")
	#### Add one to the number of samples to account for the normal/diploid cell 
	output_f.write("dimensions ntax="+str(len(samples))+";"+"\n"+"\t")
	output_f.write("\t"+"taxlabels"+"\n")
	for cell in samples:
		output_f.write(cell+"\n")
	##### Add one more taxon as the diploid/normal cell 
	output_f.write(";"+"\n"+"end"+";"+"\n")
	output_f.write("begin characters;"+"\n"+"\t"+"dimensions nchar="+str(bin_array.shape[1])+";"+"\n"+"\t")
	output_f.write("format datatype=standard symbols=\"0~9 A~Z\";"+"\n"+"\t")
	output_f.write("matrix"+"\n")
	for i in range(bin_array.shape[0]):
		output_f.write(samples[i]+"\t")
		for j in range(bin_array.shape[1]):
			output_f.write(str(bin_array[i][j]))
		output_f.write('\n')
	##### Add the diploid cell to the matrix with all the bins' states equal to 2
	##### Specify the number of best trees 

	output_f.write(";"+"\n"+"end;"+"\n")
	output_f.write("begin paup;"+"\n")
	# output_f.write("\t"+"out group diploid;"+"\n")
	output_f.write("\t"+"NJ "+";"+"\n")
	output_f.write("\t"+"DerootTrees;"+"\n")
	# output_f.write("\t"+"RootTrees;"+"\n")
	output_f.write("\t"+"log start file=NJ.log replace;"+"\n")
	output_f.write("\t"+"DescribeTrees 1 /xout=Both plot=none brLens=yes;"+"\n")
	output_f.write("\t"+"log stop;"+"\n")
	output_f.write("\t"+"savetrees file=NJ.tre replace;"+"\n")
	output_f.write("\t"+"quit;"+"\n")
	output_f.write("end;")
	output_f.close()
	return 

def Parse_Ginkgo(address):
	sample_names = None
	bin_arr = []
	with open(address,"r") as f:
		for line in f:
			line = line.strip()
			if "CHR" in line:
				raw_names = line.split()[3:]
				sample_names = [name.strip() for name in raw_names]
			else:
				raw_states = line.split()[3:]
				# print len(set(raw_states))
				##### any state greater than or equal to 15 will be named as "F"
				##### states from 10 to 14 are named as "A" to "E"
				bin_arr.append([Nexus_state(state) for state in raw_states])
	bin_arr = np.array(bin_arr)
	return (bin_arr.T,sample_names)

def Parse_HMMcopy(address):
	return Parse_Ginkgo(address)

def call_PAUP(path_2_paup_dir, exec_name, input_address):
	# cmd = "cd "+path_2_paup_dir
	# os.system(cmd)
	##### now we are in the paup dir
	##### run paup on the nexus file(s) present in the directory
	cmd = path_2_paup_dir+"/"+exec_name+" "+input_address
	os.system(cmd)
	##### log file(s) ready

	return 
if __name__ == "__main__":
	#######################################################################
	####################### Parsing the arguments #########################
	#######################################################################
	infile = ""
	paup_path = ""
	gt_path = ""
	in_dir = "/Users/edrisi/Documents/CNV_medicc_project/data/march/"
	ap = argparse.ArgumentParser()
	ap.add_argument("-paup","--path to PAUP",required=True, help="Path to the command-line version of PAUP")

	args = vars(ap.parse_args())

	if args['path to PAUP']!=None:
		paup_path = args['path to PAUP']

	executive_file_name = ""
	working_dir = ""
	arr = paup_path.strip().split('/')
	executive_file_name = arr[-1]
	working_dir = '/'.join(map(str, arr[0:len(arr)-1]))
	working_dir += '/'

	df = pd.DataFrame()
	for ploidy in range(5):
		tmp_arr = []
		for repetition in range(5):
			infile = in_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/gt.cnp"
			(bin_array,samples) = Parse_Ginkgo(address=infile)

			NJ(bin_array=bin_array, samples=samples, output_address=working_dir+"NJ.nex")
			call_PAUP(path_2_paup_dir=working_dir, exec_name=executive_file_name, input_address=working_dir+"NJ.nex")
			gt_path = in_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/gt.newick"
			gt_tr = parse.read_newick(filename=gt_path)
			newick_str = ""
			translate = {}
			with open(working_dir+"NJ.tre", 'r') as result_f:
				for line in result_f:
					if "[&U]" in line:
						newick_str = line.strip().split("[&U] ")[-1]
					if "leaf" in line:
						val = line.strip().split(" ")[1]
						key = line.strip().split(" ")[0]
						if val.endswith(","):
							val = val.replace("," , "")
							translate[key]=val
						else:
							translate[key]=val
			paup_tr = Tree(newick_str, format=8)
			for leaf in paup_tr:
				leaf.name = translate[leaf.name]
			gt_tr.unroot()
			paup_tr.unroot()
			Distance = gt_tr.robinson_foulds(paup_tr, unrooted_trees=True)[0]
			print("RF distance: ", Distance)
			tmp_arr.append(Distance)
		df['ploidy'+str(ploidy+1)] = tmp_arr
	sns.boxplot(x="variable", y="value", data=pd.melt(df))
	plt.show()
				



