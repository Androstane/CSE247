import numpy as np
import os 
import sys
from os import walk 
import argparse


class UnAcceptedCSVLine(Exception):
	def __init__(self, data):
		self.data = data
	def __str__(self):
		return repr(self.data)

def checkEqualto2(arr):
	return list(set(arr))==["2"] 

def Nexus_state(the_string):
	#### The input must be an integer 
	#### The output is the character suitable for PAUP
	cut = 15
	alphabets = ["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"]
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

def to_Nexus(bin_array, samples, output_address, swapping, initial_tree, search_mode):
	output_f = open(output_address,"w")
	output_f.write("#NEXUS"+'\n')
	output_f.write("begin taxa;"+"\n"+"\t")
	#### Add one to the number of samples to account for the normal/diploid cell 
	output_f.write("dimensions ntax="+str(len(samples)+1)+";"+"\n"+"\t")
	output_f.write("\t"+"taxlabels"+"\n")
	for cell in samples:
		output_f.write(cell+"\n")
	##### Add one more taxon as the diploid/normal cell 
	output_f.write("diploid"+"\n")
	output_f.write(";"+"\n"+"end"+";"+"\n")
	output_f.write("begin characters;"+"\n"+"\t"+"dimensions nchar="+str(bin_array.shape[1])+";"+"\n"+"\t")
	output_f.write("format datatype=standard symbols=\"0~9 A~F\";"+"\n"+"\t")
	output_f.write("matrix"+"\n")
	for i in range(bin_array.shape[0]):
		output_f.write(samples[i]+"\t")
		for j in range(bin_array.shape[1]):
			output_f.write(str(bin_array[i][j]))
		output_f.write('\n')
	##### Add the diploid cell to the matrix with all the bins' states equal to 2
	##### Specify the number of best trees 
	nbest = 10

	treelist = [i+1 for i in range(nbest)]
	output_f.write("diploid"+"\t")
	for j in range(bin_array.shape[1]):
		output_f.write(str(2))
	output_f.write('\n')
	output_f.write(";"+"\n"+"end;"+"\n")
	output_f.write("begin paup;"+"\n")
	output_f.write("\t"+"outgroup diploid;"+"\n")
	output_f.write("\t"+initial_tree+";"+"\n")
	output_f.write("\t"+"DerootTrees;"+"\n")
	output_f.write("\t"+"set increase=auto;"+"\n")
	output_f.write("\t"+"set criterion=parsimony;"+"\n")
	if search_mode=='steepest':
		output_f.write("\t"+"Hsearch nbest="+str(0)+" swap="+swapping+" start=current steepest=yes;"+"\n")
	else:
		output_f.write("\t"+"Hsearch nbest="+str(nbest)+" swap="+swapping+" start=current steepest=no;"+"\n")
	output_f.write("\t"+"RootTrees;"+"\n")
	output_f.write("\t"+"log start file="+swapping+"_"+initial_tree+"_"+search_mode+".log replace;"+"\n")
	if search_mode=='steepest':
		output_f.write("\t"+"DescribeTrees /xout=Both plot=none brLens=yes root=outgroup;"+"\n")
	else:
		output_f.write("\t"+"DescribeTrees "+' '.join(map(str,treelist))+" /xout=Both plot=none brLens=yes root=outgroup;"+"\n")
	output_f.write("\t"+"log stop;"+"\n")
	output_f.write("\t"+"savetrees file=st_"+swapping+"_"+initial_tree+"_"+search_mode+".tre replace;"+"\n")
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
				if checkEqualto2(raw_states):
					pass
					# print(line)
				else:
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
	ap = argparse.ArgumentParser()
	ap.add_argument("-paup","--path to PAUP",required=True, help="Path to the command-line version of PAUP")
	ap.add_argument("-in","--path to input",required=True, help="Path to the input file containing the absolute copy number values for all the bins")
	args = vars(ap.parse_args())

	if args['path to PAUP']!=None:
		paup_path = args['path to PAUP']
	if args["path to input"]!=None:
		infile = args["path to input"]

	# working_dir = "/Users/edrisi/Documents/CNV_project/scripts"
	# Ginkgo_path = "/Users/edrisi/Documents/CNV_project/realdata/Ginkgo/sample102/sample102.SegCopy"
	# HMM_path = "/Users/edrisi/Documents/CNV_project/realdata/HMMcopy/sample615/SegCopy.fromsegscsv.csv"
	# nexus_output = working_dir+"/nexus_input_hmm_sample615.nex"
	# Parse_HMMcopy(HMM_path)

	# (bin_array,samples) = Parse_Ginkgo(address=Ginkgo_path)
	# (bin_array,samples) = Parse_HMMcopy(address=HMM_path)
	executive_file_name = ""
	working_dir = ""
	arr = paup_path.strip().split('/')
	executive_file_name = arr[-1]
	working_dir = '/'.join(map(str, arr[0:len(arr)-1]))
	working_dir += '/'
	(bin_array,samples) = Parse_Ginkgo(address=infile)


	initial_trees = ['NJ','UPGMA']
	swapping_methods = ['SPR','TBR','NNI']
	mode = ['steepest', 'nbest']
	# initial_trees = ['NJ']
	# swapping_methods = ['NNI']
	for m in mode:
		for tree in initial_trees:
			for method in swapping_methods:
				to_Nexus(bin_array=bin_array,samples=samples,output_address=working_dir+"experiment_"+tree+"_"+method+"_"+m+".nex", swapping=method, initial_tree=tree, search_mode=m)
				call_PAUP(path_2_paup_dir=working_dir, exec_name=executive_file_name, input_address=working_dir+"experiment_"+tree+"_"+method+"_"+m+".nex")
	# print(bin_array.shape)
