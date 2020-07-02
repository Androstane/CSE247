import numpy as np 
import ete3 
from ete3 import Tree, TreeStyle
import argparse
import random
import sys
import copy


def list_branches(tree_):
	branch_arr = []
	for node in tree_.traverse("postorder"):
		if not node.is_leaf():
			if len(node.get_children())==2:
				a,b = node.get_children()
				a.detach()
				# if not a.is_leaf() and len(a)>2 and (len(tree_)+1)>2:
				if len(a)>2:
					branch_arr.append((node, a))
				node.add_child(a)
				b.detach()
				# if not b.is_leaf() and len(b)>2 and (len(tree_)+1)>2:
				if len(b)>2:
					branch_arr.append((node, b))
				node.add_child(b)
			else:
				a = node.get_children()[0]
				a.detach()
				# if not a.is_leaf() and len(a)>2 and (len(tree_)+1)>2:
				if len(a)>2:
					branch_arr.append((node, a))
				node.add_child(a)
	if len(branch_arr)==0:
		return None
	else:
		return branch_arr

def select_branch_exclusively(tree__, brnch_arr):
	branch_arr = []
	# print brnch_arr
	for node in tree__.traverse("postorder"):
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

def perform_TBR(tree, selection):
	trees = []
	tr = tree
	# print tr
	(parent, selected_node) = selection
	# print (parent, selected_node)
	selected_node.detach()
	##### we are sure that the tree under the selected node has len>=2
	##### so, we can assume that it has two children 
	forbidden_arr_ = [(selected_node.name, selected_node.get_children()[0].name), (selected_node.name, selected_node.get_children()[1].name)]
	if len(tr)==1:
		if parent.is_root():
			##### keep the root and use this node as the only option left on T1
			print "root is the T1"
			option_arr = select_branch_exclusively(selected_node,forbidden_arr_)
			for k in random.sample(range(len(option_arr)), 1):
			# for k in range(1):
				tr1_cpy = copy.deepcopy(tr)
				tr2_cpy = copy.deepcopy(selected_node)
				(p1,c_1) = select_branch_exclusively(tr2_cpy,forbidden_arr_)[k]
				# print (p1, c_1)
				tr2_cpy.set_outgroup(c_1)
				# print c_1.up
				tr1_cpy.add_child(c_1.up)
				# print tr1_cpy
				trees.append(tr1_cpy)
				# print parent.name
				# for node in parent.traverse():
					# print node.name
					# print node.get_children()
		else:
			##### delete the parent node, then it remains only one branch to select from T1
			print "there is only one branch left on T1"
			# parent.delete()
			option_arr = select_branch_exclusively(selected_node,forbidden_arr_)
			for l in random.sample(range(len(option_arr)), 1):
				tr1_cpy = copy.deepcopy(tr)
				tr2_cpy = copy.deepcopy(selected_node)
				(p1,c_1) = select_branch_exclusively(tr2_cpy,forbidden_arr_)[l]
				# print (p1, c_1)
				tr2_cpy.set_outgroup(c_1)
				(tr1_cpy&str(parent.name)).add_child(c_1.up)
				# print tr1_cpy
				trees.append(tr1_cpy)
				# for node in tr1_cpy.traverse():
					# print node.name
					# print node.get_children()

	else:
		# print selected_node
		# print tr
		option_arr_t2 = select_branch_exclusively(selected_node,forbidden_arr_)
		parent.delete()
		option_arr_t1 = select_branch_exclusively(tr, [])
		for br_indx_t2 in random.sample(range(len(option_arr_t2)), 1):
			for br_indx_t1 in random.sample(range(len(option_arr_t1)), 1):
				tr1_cpy = copy.deepcopy(tr)
				tr2_cpy = copy.deepcopy(selected_node)
				(p1,c_1) = select_branch_exclusively(tr2_cpy,forbidden_arr_)[br_indx_t2]
				# print (p1, c_1)
				tr2_cpy.set_outgroup(c_1)
				(p2, c_2) = select_branch_exclusively(tr1_cpy, [])[br_indx_t1]
				# print (p2, c_2)
				c_2.detach()
				in_node = p2.add_child()
				in_node.add_child(c_2)
				in_node.add_child(c_1.up)
				# print tr1_cpy
				# for node in tr1_cpy.traverse():
				# 	print node.name
				# 	print node.get_children()
				# print "------------------------------------"
				trees.append(tr1_cpy)
	return trees

def Main(in_tree, N):
	''' Since it is not possible to deal with all the TBR rearrangements 
	given a particular topology, we randomly select N proposed topologies from 
	all the possible ones '''

	tree_list = []
	tr = in_tree
	# print tr
	if len(tr)<=3:
		print("ERROR: Tree bisection reconnection cannot be applied on this tree")
		return tree_list

	lst_branches = list_branches(tr)
	# print len(lst_branches)
	if len(lst_branches)==0:
		print("No valid SPR rearrangement for the given tree")
		return tree_list
	else:
		for i in random.sample(range(len(lst_branches)), k=N):
		# for i in range(1):
			tmp_cpy = copy.deepcopy(tr)
			tmp_lst = list_branches(tmp_cpy)
			# print "**************************************"
			result_trs = perform_TBR(tree=tmp_cpy, selection=tmp_lst[i])
			if result_trs == None:
				print("This rearrangement resulted in the same tree")
			else:
				tree_list.extend(result_trs)
	# for Tr in tree_list:
		# print Tr
	return tree_list

if __name__=="__main__":
	### Toy example

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
