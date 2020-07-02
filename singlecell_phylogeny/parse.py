import numpy as np
import ete3
from ete3 import Tree

def check_state(the_string):
	#### The input must be an integer 
	#### The output is the character suitable for PAUP
	cut = 36
	state = None
	if int(the_string)<=cut-1:
		state = int(the_string)
	else:
		state = cut-1
	return state

def checkEqualto2(arr):
	return list(set(arr))==["2"] 

def chr_extract(chr):
	if chr[-1]=="X" or chr[-1]=="x":
		return 23
	elif chr[-1]=="Y" or chr[-1]=="y":
		return 24
	else:
		chrString = ""
		for i in chr:
			if i.isdigit():
				chrString+=i
		return int(chrString)

def Parse_WG_wt_sep(filename):
	sample_names = None
	bin_arr = []
	chr_number = -1
	with open(filename,"r") as f:
		for line in f:
			line = line.strip()
			if "CHR" in line:
				raw_names = line.strip().split()[3:]
				sample_names = [name.strip() for name in raw_names]
			else:
				arr = line.strip().split()
				raw_states = arr[3:]
				current_chr = chr_extract(arr[0])
				if chr_number==-1:
					chr_number = current_chr
					bin_arr.append([check_state(state) for state in raw_states])
				elif chr_number!=current_chr:
					bin_arr.append([2 for i in range(len(sample_names))])
					bin_arr.append([2 for i in range(len(sample_names))])
					bin_arr.append([check_state(state) for state in raw_states])
				else:
					bin_arr.append([check_state(state) for state in raw_states])
	bin_arr = np.array(bin_arr)
	return (bin_arr.T,sample_names)

def Parse_input_simulated(filename):
	''' reads the simulated input data from a file 
	returns an array containing the absolute copy numbers for each cell,
	each chromosome, and each bin '''
	''' We assume the chromosomes are in order '''
	chr_name = ""
	counter = 1
	chr_counter = -1
	seqs = None
	with open(filename, "r") as infile:
		for line in infile:
			if "CHR" in line:
				names = line.strip().split('\t')[3:]
				names = [name.strip() for name in names]
				seqs = [[] for i in range(len(names))]
			else:
				lst = line.strip().split('\t')
				tmp_lst = [str_.strip() for str_ in lst[3:]]
				if lst[0]!=chr_name:
					chr_counter+=1
					chr_name = lst[0]
					for smpl in range(len(names)):
						seqs[smpl].append([])
						seqs[smpl][chr_counter].append(check_state(tmp_lst[smpl]))
				else:
					for smpl in range(len(names)):
						seqs[smpl][chr_counter].append(check_state(tmp_lst[smpl]))
	seqs = np.array(seqs)
	return (seqs,names)

def Parse_input_real(filename):
	return Parse_input_simulated(filename)

def read_newick(filename):
	''' reads a newick string from a file '''
	in_f = open(filename, 'r')
	str_ = in_f.readline()
	str_ = str_.strip()
	if str_[-1]!=";":
		str_+=";"
	#### format=8 is applied when there are only node names on the tree
	t = Tree(str_, format=8)
	# print t.get_tree_root().name
	for leaf in t:
		leaf.name = "leaf"+leaf.name
	return t

def read_dist_mat(filename):
	dist_mat = []
	with open(filename, "r") as file:
		for line in file:
			### we assume the distances are integers and tab separated
			dist_mat.append([int(i) for i in line.strip().split('\t')])
	return dist_mat

def read_names(filename):
	names = []
	with open(filename,'r') as file:
		for line in file:
			if line.strip():
				names.append("leaf"+line.strip())
	return names

if __name__ == "__main__":
	read_newick(filename="/Users/edrisi/Documents/CNV_medicc_project/data/simulated/from_first_step.tree.newick")
	sequences, leafnames =Parse_input_simulated(filename="/Users/edrisi/Documents/CNV_medicc_project/data/simulated/gt.all.csv.segcopy")
	# print sequences.shape
	# print sequences[:,0]




