import numpy as np 
import ete3 
from ete3 import Tree, TreeStyle
import argparse
import shlex
import random
import sys
import copy
import SPR
import NNI
import TBR
import tree 
import distance 
import uuid 
import parse
import collections
import matplotlib.pyplot as plt
from matrix_dist import rel
from matrix_dist import fb

def Sort_Dict(d):
	return sorted(d.items(), key=lambda x: x[1])

def Merge(dict1, dict2):
	res = dict1.copy()
	res.update(dict2)
	return res

def change_leafnames(newick_, dictionary_):
	new_newick = newick_
	for key in dictionary_:
		new_newick = new_newick.replace(key, str(dictionary_[key]))
	return new_newick

def ThereAreDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def gen_unique_ids(num):
	ids = []
	while True:
		ids = []
		for i in range(num):
			ids.append(str(uuid.uuid4()))
		if not ThereAreDuplicates(ids):
			break
	# ids = random.sample(range(100), num)
	return ids

def Check_Tree(tre_):
	return_tr = None
	''' read a Newick string or an instance of Tree class 
	and return a Tree object with all the nodes having names '''
	''' check if the input is a Newick string'''
	if isinstance(tre_, str):
		return_tr = Tree(tre_, format=3)
	else:
		return_tr = tre_
	nameless_count = 0
	if not return_tr.get_tree_root().name:
		return_tr.get_tree_root().name = 'diploid'
	for node in return_tr.traverse():
		if not node.name:
			nameless_count+=1
	new_ids = gen_unique_ids(nameless_count)
	id_idx = 0
	for node in return_tr.traverse():
		node.dist = 0
		if not node.name:
			node.name = new_ids[id_idx]
			id_idx+=1
	return return_tr

def GET_ALL_SCORES(tr_lst, dict_names):
	''' the input is a list of ete trees '''
	''' it returns a dictionary in which the key is the tree and 
	the score is the value '''
	tr_dict = {}
	for elm in tr_lst:
		elm = Check_Tree(elm)
		elm_cpy = copy.deepcopy(elm)
		new_root = (elm_cpy.get_tree_root()).get_children()[0]
		new_root.detach()
		new_root.name = "diploid"
		new_newick_ = change_leafnames(new_root.write(format=5), dict_names)
		# print new_newick_
		score = tree.tree_score(new_newick_)
		tr_dict[elm.write(format=3)] = score
	return tr_dict

def Jump_TBR(trr, num_jumps):
	### make a copy of the original tree
	original_tr = copy.deepcopy(Check_Tree(trr))
	tmp_tr = copy.deepcopy(original_tr)
	count = 0
	while True:
		tmp_tr_lst = TBR.Main(tmp_tr, 1)
		if len(tmp_tr_lst)==0:
			pass
		else:
			count+=1
			tmp_tr = tmp_tr_lst[0]
			# print tmp_tr
		if count== num_jumps:
			if Check_Tree(original_tr).get_topology_id() == Check_Tree(tmp_tr).get_topology_id():
				count = 0
				print("reset--------")
			else:
				print("finished---------")
				break
	return Check_Tree(tmp_tr)

if __name__ == "__main__":
	#######################################################################
	####################### Parsing the arguments #########################
	#######################################################################
	in_path = ""
	iter_num = 10
	ap = argparse.ArgumentParser()
	# ap.add_argument("-ms","--minimum coverage",required=False, help="The minimum number of reads required to support the alternative")
	# ap.add_argument("-nmc","--minimum cells",required=False,
	# 	help="The number of cells having a minimum number of alternative reads needed to consider a genomic location a mutation candidate locus")
	# ap.add_argument("-thr","--missing data threshold",required=False, help="The minimum coverage below which the query is a missing data")
	# ap.add_argument("-fn","--false negative rate",required=False, help="False negative error rate")
	ap.add_argument("-iter","--number of iterations",required=False, help="Number of iterations in search for optimal tree")
	ap.add_argument("-in","--path to input",required=True, help="Path to the input file containing the absolute copy number values for all the bins")
	ap.add_argument("-gt","--path to ground truth",required= True, help="Path to the newick string containing the ground truth tree")
	ap.add_argument("-gtm", "--path to ground truth matrix", required = False, help = "path to the ground truth distance matrix")
	# ap.add_argument("-n","--number of cells",required=True, help="The number of the cells")
	args = vars(ap.parse_args())
	# if args['minimum coverage']!=None:
	# 	ms = args['minimum coverage']
	# if args['minimum cells']!=None:
	# 	nmc = args['minimum cells']
	# if args['missing data threshold']!=None:
	# 	missing_data_threshold = args['missing data threshold']
	# if args['false negative rate']!=None:
	# 	fn = args['false negative rate']
	if args['number of iterations']!=None:
		iter_num = int(args['number of iterations'])
	if args["path to input"]!=None:
		in_path = args["path to input"]
	if args["path to ground truth"]!=None:
		gt_path = args["path to ground truth"]
	if args["path to ground truth matrix"] != None:
		gtm_path = args["path to ground truth matrix"]

	##########################################################################
	############################ initialization ##############################
	##########################################################################
	seqs_parsed, names = parse.Parse_input_simulated(filename=in_path)
	# seqs_parsed = seqs_parsed[0:20,:]
	# names = names[0:20]
	print(len(names))
	print(seqs_parsed.shape)
	# seqs_parsed, names = parse.Parse_input_simulated(filename="/Users/edrisi/Documents/CNV_medicc_project/data/simulated/gt.all.csv.segcopy")
	seqs_ = []
	for indx in range(seqs_parsed.shape[0]):
		seqs_.append(seqs_parsed[:,0,:][indx].tolist())
		#seqs_.append(seqs_parsed[:,0][indx])
		#print(shlex.split(str(seqs_parsed[:,0][indx]).strip('[]')))
		#print(seqs_parsed[:,0][indx])
		#print(seqs_parsed[:,0][indx])
	if 'diploid' not in names:
		names.append('diploid')
	seqs_.append([2 for Bin in range(len(seqs_[0]))])


	# seqs_ = [[2,1,1,1,2],[4,5,4,4,3],[3,0,2,5,3],[1,0,0,4,3],[2,3,3,4,2],[0,1,1,2,1],[2,2,2,2,2]]
	# names = ['F','E','D','C','A','B', 'diploid']
	my_dict = {}
	for k in range(len(names)):
		my_dict[names[k]] = seqs_[k]

	string_frst, hier_tr,dist_matrix = tree.compute_tree(seqs_)
	dist_d = dist_matrix
	dist_matrix = dist_matrix[0:-1, 0:-1]
	gt_matrix = np.genfromtxt(gtm_path, delimiter = "\t")
	gt_matrix = gt_matrix[:,0:-1]
	#print(names)
	print(gt_matrix.shape)
	print(dist_matrix.shape)
	print("distance to gt matrix:", rel(dist_matrix, gt_matrix), fb(dist_matrix, gt_matrix))
	string_scnd = tree.getNewick(hier_tr, "", hier_tr.dist, names)
	upgma_tr = Tree(string_scnd)
	newroot = (upgma_tr&"diploid").up
	upgma_tr.set_outgroup(upgma_tr&"diploid")
	(upgma_tr&"diploid").detach()
	tr2 = Tree(name="diploid")

	tr2.add_child(upgma_tr)
	newroot.delete()
	tr2 = Check_Tree(tr2)
	#print(tr2)
	#groundtruth_tr = parse.read_newick(gt_path)
	#groundtruth_tr.unroot()
	tr2.unroot()
	#print("distance is:", groundtruth_tr.robinson_foulds(tr2, unrooted_trees=True)[0])

	##########################################################################
	################################# Search #################################
	##########################################################################
	#### max size of the stack
	max_stack_len = 10
	scores = []
	#### The number of TBR rearrangements to perform when jumping to new hill
	jump_number = 5
	#### The number of new topologies returned from each rearrangement 
	N_tplgy = 2
	#### Total number of iterations 
	T = iter_num
	#### The array to save the N_best number of the best results 
	best_stack = GET_ALL_SCORES([tr2], my_dict)
	best_score = best_stack[next(iter(best_stack))]
	#### The number of consecutive iterations in which no improvement is observed
	Waiting_time = 1000
	#### counter for local optima alert 
	local_alert = 0
	#### probability distribution of the tree rearrangements
	p_NNI = 0.8
	p_SPR = 0.1
	p_TBR = 0.1
	#### the current trees used for the next iteration
	stack_ = copy.deepcopy(best_stack)
	rearrangements = ['NNI', 'SPR', 'TBR']
	distribution = [p_NNI,p_SPR,p_TBR]

	# setup toolbar
	# sys.stdout.write("[%s]" % (" " * T))
	# sys.stdout.flush()
	# sys.stdout.write("\b" * (T+1)) 

	for iteration in range(T):
		print("iteration number %s" %iteration)
		tmp_dict = {}
		tr_list = None
		if len(stack_)>max_stack_len:
			stack_ = dict(list(stack_.items())[:max_stack_len]) 
		if local_alert<Waiting_time:
			for tr_str in stack_:
				tr_ = Check_Tree(tr_str)
				rearrangement_type = np.random.choice(rearrangements, p=distribution)
				print(rearrangement_type)
				if rearrangement_type=='NNI':
					''' propose new topologies '''
					tr_list = NNI.Main(tr_, N=N_tplgy)
				elif rearrangement_type=='TBR': 
					tr_list = TBR.Main(tr_, N=N_tplgy)
				else:
					tr_list = SPR.Main(tr_, N=N_tplgy)
				''' compute the scores of the new topologies '''
				''' merge the results with the trees in tmp_dict '''
				tmp_dict = Merge(GET_ALL_SCORES(tr_list, my_dict), tmp_dict)
			# tmp_dict = dict(Sort_Dict(tmp_dict))
		#### here we need to perform multiple TBR rearrangements on each of the trees in the stack
		else:
			#### reset the counter 
			local_alert = 0
			#### Save all the topology ids in a list 
			new_list = []
			topology_ids = []
			for tree_key in stack_:
				tplgy_ids.append(Check_Tree(tree_key).get_topology_id())
			tplgy_id_set = set(tplgy_ids)
			#### iterate over all the current topologies in stack_
			for topolg in stack_:
				jumped_tr = Jump_TBR(trr=topolg, num_jumps=jump_number)
				while jumped_tr.get_topology_id() in tplgy_id_set:
					jumped_tr = Jump_TBR(trr=topolg, num_jumps=jump_number)
				#### a new topology has been generated which was not previously in the set of topologies
				tplgy_id_set.add(jumped_tr.get_topology_id())
				new_list.append(jumped_tr)
			### replace the previous stack_ with the dictionary from the new list of trees
			stack_ = {}
			stack_ = GET_ALL_SCORES(new_list, my_dict)
			tmp_dict = Merge(GET_ALL_SCORES(new_list, my_dict), tmp_dict)
		''' prepare for next iteration '''
		''' check if there is any improvement in terms of the score '''
		sorted_tmp = Sort_Dict(tmp_dict)
		''' best score from the tmp dictionary '''
		best_of_tmp = sorted_tmp[0][1]
		print(best_of_tmp)
		print(best_score)
		scores.append(best_score)
		''' do the following steps when facing a new topology 
		with better score '''
		if best_of_tmp < best_score:
			print('Found Better Score........................')
			best_score = best_of_tmp
			stack_ = {}
			best_stack = {}
			for element in sorted_tmp:
				if element[1]==best_score:
					stack_[element[0]] = element[1]
					best_stack[element[0]] = element[1]
					#### ALERT ####
					#### keeping only the first one, remove the break command to save all of them 
					# break
			#### reset the counter because a better score weas found 
			local_alert = 0

		elif best_of_tmp==best_score:
			tplgy_ids = []
			for tr_key in stack_:
				tplgy_ids.append(Check_Tree(tr_key).get_topology_id())
			tplgy_id_set = set(tplgy_ids)
			print("Found Equal Score.........................")
			for element in sorted_tmp:
				if element[1]==best_score and (Check_Tree(element[0]).get_topology_id() not in tplgy_id_set):
					# print("new topology was added to the list of best trees")
					##### ALERT #####
					##### keeping only one of them, remove the break command to save all the options
					##### ALERT ######
					##### reseting the stack_ to keep only one of the new ones
					##### remove the following line to merge all the previous ones with the new ones
					# stack_ = {}
					# best_stack = {}
					stack_[element[0]] = element[1]
					best_stack[element[0]] = element[1]
					# break
			#### increase the counter since the score was not improved
			local_alert+=1
		else:
			print("No Better Score Was Found.................")
			#### increase the counter since the score was not improved
			local_alert+=1
			pass
		# update the bar
		# sys.stdout.write("-")
		# sys.stdout.flush()
	# sys.stdout.write("]\n")
	groundtruth_tr = parse.read_newick(gt_path)
	groundtruth_tr.unroot()
	for tr_key in best_stack:
		tree_format = Check_Tree(tr_key)
		tree_format.unroot()
		print("distance is:", groundtruth_tr.robinson_foulds(tree_format, unrooted_trees=True)[0])

	print(best_stack)



