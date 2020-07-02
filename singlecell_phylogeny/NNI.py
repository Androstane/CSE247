import numpy as np 
import ete3 
from ete3 import Tree, TreeStyle
import argparse
import random
import sys
import copy

def has_internal_branch(node):
	######## This function checks if the given node has 
	######## any internal branches and returns a boolean value 
	if node.is_root():
		return False 
	elif node.is_leaf():
		return False
	elif len(node.get_children())==1:
		return False
	else:
		a,b = node.get_children()
		if a.is_leaf() and b.is_leaf():
			return False
		else:
			return True

def list_branches(tree):
	branch_arr = []
	for node in tree.traverse("postorder"):
		if has_internal_branch(node):
			a,b = node.get_children()
			if not a.is_leaf():
				branch_arr.append((node, a, b))
			if not b.is_leaf():
				branch_arr.append((node, b, a))
	if len(branch_arr)==0:
		return None
	else:
		return branch_arr

def perform_NNI(tree, selection, choice_):
	######## The root of the tree is treated as a leaf during 
	######## the move but it still remains as the root of the tree after the move 
	Z = tree
	selected_node = None
	D_component = None
	left_child = None
	A_component = None

	(selected_node, child, A_component) = selection
	D_component = selected_node.up
	##### select which exchange to perform 
	if choice_==1:
		#### the first type of exchange
		B_component = child.get_children()[0]
		C_component = child.get_children()[1]
	else:
		#### the second type of exchange
		B_component = child.get_children()[1]
		C_component = child.get_children()[0]
	B_component.detach()
	C_component.detach()
	A_component.detach()
	selected_node.add_child(B_component)
	child.add_child(A_component)
	child.add_child(C_component)
	return Z
def Main(in_tree, N):
	''' Since it is not possible to deal with all the 
	NNI rearrangements given a topology, we randomly sample 
	2N tree from all the possible rearrangements in each NNI.py call'''
	tree_list = []
	tr = in_tree
	lst_branches = list_branches(tree=tr)
	if len(lst_branches)==0:
		print("No valid NNI rearrangement for the given tree")
		return tree_list
	else:
		for i in random.sample(range(len(lst_branches)), k=N):
			tmp_cpy_1 = copy.deepcopy(tr)
			tmp_cpy_2 = copy.deepcopy(tr)
			tmp_lst_1 = list_branches(tmp_cpy_1)
			tmp_lst_2 = list_branches(tmp_cpy_2)
			tree_list.append(perform_NNI(tree=tmp_cpy_1, selection=tmp_lst_1[i], choice_=1))
			tree_list.append(perform_NNI(tree=tmp_cpy_2, selection=tmp_lst_2[i], choice_=2))
	# for Tr in tree_list:
		# print Tr
	return tree_list

if __name__=="__main__":
	### Toy example
	### When using this module, please call the Main function
	### NNI must return 8 different topologies for this toy example
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
