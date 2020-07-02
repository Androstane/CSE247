import numpy as np
import ete3
from ete3 import Tree, TreeStyle
import argparse


def update_distance(dist_matrix, index_a, index_b):
	n = dist_matrix.shape[0]
	updated_dist = np.zeros((n-1,n-1))
	### The first row of the updated matrix belongs to the new node
	### For the new node:
	counter = 1
	for j in range(n):
		if j!=index_a and j!=index_b:
			updated_dist[0][counter] = 0.5*(dist_matrix[index_a][j]+dist_matrix[index_b][j]-dist_matrix[index_a][index_b])
			updated_dist[counter][0] = updated_dist[0][counter]
			counter+=1
	### For the rest of the nodes:
	row = 1
	for i in range(n):
		if i!=index_a and i!=index_b:
			col = row+1
			for j in range(i+1,n):
				if j!=index_a and j!=index_b:
					updated_dist[row][col]=dist_matrix[i][j]
					updated_dist[col][row]=updated_dist[row][col]
					col+=1
			row+=1
	return updated_dist


def attach_new_node(distance_matrix, node_a, node_b, index_a, index_b):
	new_node = Tree(name=node_a.name+","+node_b.name)
	#### branch length estimation 
	n = distance_matrix.shape[0]
	b_a_u = 0.5*(distance_matrix[index_a][index_b]+(1.0/(n-2))*(sum(distance_matrix[index_a])-sum(distance_matrix[index_b])))
	b_b_u = distance_matrix[index_a][index_b] - b_a_u
	#### Concatenate the subtrees to the new node as the new parent
	new_node.add_child(node_a)
	new_node.add_child(node_b)
	#### Assign the estiamted branch lengths to the nodes
	node_a.dist = b_a_u
	node_b.dist = b_b_u
	return new_node


def calculate_Q(dist_matrix):
	Q_mat = np.zeros(dist_matrix.shape)
	n = Q_mat.shape[0]
	row=0
	for row in range(n):
	# while row<n:
		# col=row+1
		# while col<n:
		for col in range(row+1,n):
			Q_mat[row][col]=(n-2)*dist_matrix[row][col]-sum(dist_matrix[row])-sum(dist_matrix[col])
			# Q_mat[col][row]=Q_mat[row][col]
			# col+=1
		# row+=1
	# print Q_mat
	return Q_mat

def Main(dist_mat, names):
	#### Initialization 
	dist_mat = np.array(dist_mat)
	current_taxa = []
	for i in range(len(dist_mat)):
		current_taxa.append(Tree(name=names[i]))
	# print dist_mat
	#### Iteration
	while len(current_taxa)>2:
		# print current_taxa
		Q = calculate_Q(dist_mat)
		a_indx = 0
		b_indx = 1
		# a_indx, b_indx = np.unravel_index(Q.argmin(), Q.shape)
		min_ = Q[0][1]
		for i in range(Q.shape[0]):
			for j in range(i+1,Q.shape[0]):
				if Q[i][j]<min_:
					min_=Q[i][j]
					a_indx=i
					b_indx=j
		# print a_indx
		# print b_indx
		new = attach_new_node(dist_mat, current_taxa[a_indx], current_taxa[b_indx], a_indx, b_indx)
		a = current_taxa[a_indx]
		b = current_taxa[b_indx]
		# print a.name
		# print b.name
		current_taxa.remove(a)
		current_taxa.remove(b)
		current_taxa.insert(0,new)
		# print current_taxa
		dist_mat = update_distance(dist_mat,a_indx,b_indx)
		# print dist_mat
		
	### Final step
	root = Tree()
	root.add_child(current_taxa[0])
	root.add_child(current_taxa[1])
	# print root
	# ts = TreeStyle()
	# ts.show_leaf_name = True
	# ts.show_branch_length = True
	# ts.show_branch_support = True
	# root.show(tree_style=ts)
	# root.render("mytree.png", w=183, units="mm", tree_style=circular_style)



	tr2 = root
	# root.set_outgroup(root&"diploid")
	# # print root
	# to_attach = (root&"diploid").up
	# (root&"diploid").detach()
	# # print root
	# tr2 = Tree(name="diploid")
	# tr2.add_child(to_attach)
	# # print tr2
	# to_remove = tr2.get_children()[0]
	# to_remove.delete()
	# # print tr2
	# # counter = 0
	# # for node in tr2.traverse():
	# 	# print counter
	# 	# counter+=1
	# 	# print node.name
	return tr2

if __name__=="__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-ms","--minimum coverage",required=False, help="The minimum number of reads required to support the alternative")
	# ap.add_argument("-nmc","--minimum cells",required=False,
	# 	help="The number of cells having a minimum number of alternative reads needed to consider a genomic location a mutation candidate locus")
	# ap.add_argument("-thr","--missing data threshold",required=False, help="The minimum coverage below which the query is a missing data")
	# ap.add_argument("-fn","--false negative rate",required=False, help="False negative error rate")
	# ap.add_argument("-fp","--false positive rate",required=False, help="False positive error rate")
	# ap.add_argument("-in","--input mpileup file",required=True, help="The input file (in mpileup format)")
	# ap.add_argument("-n","--number of cells",required=True, help="The number of the cells")
	# args = vars(ap.parse_args())
	# if args['minimum coverage']!=None:
	# 	ms = args['minimum coverage']
	# if args['minimum cells']!=None:
	# 	nmc = args['minimum cells']
	# if args['missing data threshold']!=None:
	# 	missing_data_threshold = args['missing data threshold']
	# if args['false negative rate']!=None:
	# 	fn = args['false negative rate']
	# if args['false positive rate']!=None:
		# fp = args['false positive rate']
	toy_example = [[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]]
	name_list = ["diploid", "1","2","3","4"]
	print Main(dist_mat=toy_example, names=name_list)