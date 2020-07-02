import numpy as np 
import ete3 
from ete3 import Tree, TreeStyle
import argparse
import random
import sys
import copy

def destination_branches(tree, brnch_arr):
	branch_arr = []
	# print brnch_arr
	for node in tree.traverse("postorder"):
		if len(node.get_children())!=0:
			for child in node.get_children():
				if (node.name, child.name) not in brnch_arr:
					branch_arr.append((node, child))
		else:
			continue
	if len(branch_arr)==0:
		return None
	else:
		return branch_arr

def list_branches(tree):
	branch_arr = []
	for node in tree.traverse("postorder"):
		if not node.is_root() and not node.is_leaf():
			for child in node.get_children():
				branch_arr.append((node, child))
		else:
			continue
	if len(branch_arr)==0:
		return None
	else:
		return branch_arr

def perform_SPR(in_tr, selection):
	trees = []
	(parent, selected_node) = selection 
	# print selection
	selected_node.detach()
	tuple_ = [(parent.up.name, (parent.get_children()[0]).name)]
	parent.delete()
	options_arr = destination_branches(in_tr, tuple_)
	# print options_arr
	if options_arr == None:
		return None
	else:
		''' we pick the other branch randomly '''
		for k in random.sample(range(len(options_arr)), 1):
			in_tr_cpy = copy.deepcopy(in_tr)
			selected_node_cpy = copy.deepcopy(selected_node)
			destination_arr = destination_branches(in_tr_cpy, tuple_)	
			(p,c) = destination_arr[k]
			# print (p,c)
			c.detach()
			in_node = p.add_child()
			in_node.add_child(c)
			in_node.add_child(selected_node_cpy)
			# print in_tr_cpy

			trees.append(in_tr_cpy)
		return trees

def Main(in_tree, N):
	''' Since it is not possible to deal with all the SPR rearrangements 
	given a particular topology, we randomly sample N new topologies in each call'''
	tree_list = []
	tr = in_tree
	# print in_tree
	lst_branches = list_branches(tr)
	# print len(lst_branches)
	if len(lst_branches)==0:
		print("No valid SPR rearrangement for the given tree")
		return tree_list
	else:
		for i in random.sample(range(len(lst_branches)), k=N):
			tmp_cpy = copy.deepcopy(tr)
			tmp_lst = list_branches(tmp_cpy)
			# print tmp_lst[i]
			# print "-------------------------------------------------------"
			result_trs = perform_SPR(in_tr=tmp_cpy, selection=tmp_lst[i])
			if result_trs == None:
				print("This rearrangement resulted in the same tree")
			else:
				tree_list.extend(result_trs)
	# for Tr in tree_list:
		# print Tr
	# print len(tree_list)
	return tree_list

if __name__=="__main__":
	#### NOTE: each of the internal nodes must have a name for this rearrangement
	### Toy example
	### SPR must return 64 different topologies for this toy example
	t = Tree(name="root")
	Z = t.add_child(name="Z")
	Y = Z.add_child(name="Y")
	F = Z.add_child(name="F")

	X = Y.add_child(name="X")
	V = Y.add_child(name="V")

	A = V.add_child(name="A")
	B = V.add_child(name="B")

	E = X.add_child(name="E")

	W = X.add_child(name="W")
	D = W.add_child(name="D")
	C = W.add_child(name="C")
	# print(t)
	Main(in_tree=t)
