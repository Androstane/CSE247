import ete3 
from ete3 import TreeStyle
from ete3 import Tree
import parse
import NJ
import matplotlib.pyplot as plt
import seaborn as sns 
import pandas as pd
import dendropy as dnd
sns.set_style("darkgrid")


def write_csv(dist_matrix, names_, dir_):
	outf = open(dir_+"dendropy.csv","w")
	st = ""
	for name in names_:
		st+=(","+name)
	outf.write(st+"\n")	
	for i in range(len(dist_matrix)):
		st=names_[i]
		for val in dist_matrix[i]:
			st+=(","+str(val))
		outf.write(st+"\n")
	outf.close()


if __name__ == "__main__":
	working_dir = "/Users/edrisi/Documents/CNV_medicc_project/data/march/"
	df = pd.DataFrame()
	for ploidy in range(5):
		tmp_arr = []
		for repetition in range(5):
			gt_path = working_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/gt.newick"
			dist_mat_path = working_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/gt.pairwisedist.txt"
			names_path = working_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/leaves.txt"
			gt_tr = parse.read_newick(filename=gt_path)
			# gt_tr.get_tree_root().name = "diploid"
			dist_mat_gt = parse.read_dist_mat(filename=dist_mat_path)
			names = parse.read_names(filename=names_path)
			write_csv(dist_matrix=dist_mat_gt, names_=names, dir_=working_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/")
			pdm = dnd.PhylogeneticDistanceMatrix.from_csv(src=open(working_dir+"p"+str(ploidy+1)+"/rep"+str(repetition+1)+"/dendropy.csv"),delimiter=",")
			nj_tree = pdm.nj_tree()

			dnd_newick = nj_tree.as_string("newick", suppress_edge_lengths=True)
			dnd_newick = dnd_newick.replace("[&U] ", "")
			# print dnd_newick
			dendropy_tr = Tree(dnd_newick, format=8)
			dendropy_tr.unroot()
			gt_tr.unroot()
			# print len(dendropy_tr)
			# print len(gt_tr)
			NJ_tr = NJ.Main(dist_mat=dist_mat_gt, names=names)
			NJ_tr.unroot()
			Distance = gt_tr.robinson_foulds(NJ_tr, unrooted_trees=True)[0]
			Alternate_distance = gt_tr.robinson_foulds(dendropy_tr, unrooted_trees=True)[0]
			tmp_arr.append(Alternate_distance)
			print("RF distance :", Distance)	
			print("RF distance from Dendropy :" ,Alternate_distance)
		df['ploidy'+str(ploidy+1)] = tmp_arr
	sns.boxplot(x="variable", y="value", data=pd.melt(df))
	plt.show()

