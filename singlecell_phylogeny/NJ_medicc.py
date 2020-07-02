import NJ
import ete3 
from ete3 import Tree
import tree
import parse 
import argparse
import dendropy as dnd

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
	in_path = ""
	gt_path = ""
	chr_ = 1
	mod = "wg"
	seqs_parsed = None
	names = None
	ap = argparse.ArgumentParser()
	ap.add_argument("-in","--path to input",required=True, help="Path to the input file containing the copy number profiles ")
	ap.add_argument("-gt","--path to the ground truth",required=True, help="Path to the ground truth tree in Newick format")
	ap.add_argument("-chr","--chromosome index",required=False, help="Selects the chromosome from the whole genome profiles, it must a positive integer. Default: 1")
	ap.add_argument("-mode","--mode of the analysis",required=False, help="Mode of the analysis could be whole genome or at chromosome level, the options are either wg or chr. Default: wg")

	args = vars(ap.parse_args())
	if args['path to input']!=None:
		in_path = args['path to input']
	if args['path to the ground truth']!=None:
		gt_path = args['path to the ground truth']
	if args['chromosome index']!=None:
		chr_ = int(args['chromosome index'])-1
	if args['mode of the analysis']!=None:
		mod = args['mode of the analysis']

	if mod=="wg" or mod=="Wg" or mod=="WG" or mod=="wG":
		mod_ = 0
		###### parsing the whole genome with separators ######
		seqs_parsed, names = parse.Parse_WG_wt_sep(filename=in_path)
	elif mod=="chr" or mod=="CHR" or mod=="Chr":
		mod_ = 1
		##### parsing the copy number profiles with separating the chromosomes #######
		seqs_parsed, names = parse.Parse_input_simulated(filename=in_path)

	seqs_ = []
	for indx in range(seqs_parsed.shape[0]):
		if mod_==0:
			##### WG with seaprators
			seqs_.append(seqs_parsed[indx].tolist())
		else:
			##### only one chromosome
			seqs_.append(seqs_parsed[indx,chr_])
	string_, hier_tr,dist_matrix = tree.compute_tree(seqs_)
	dist_matrix = dist_matrix.astype(int)
	write_csv(dist_matrix=dist_matrix, names_=names, dir_="./")
	pdm = dnd.PhylogeneticDistanceMatrix.from_csv(src=open("./dendropy.csv"),delimiter=",")
	nj_tree = pdm.nj_tree()
	dnd_newick = nj_tree.as_string("newick", suppress_edge_lengths=True)
	dnd_newick = dnd_newick.replace("[&U] ", "")
	tr_medicc = Tree(dnd_newick, format=8)
	gt_tr = parse.read_newick(filename=gt_path)
	gt_tr.unroot()
	tr_medicc.unroot()
	Distance = gt_tr.robinson_foulds(tr_medicc, unrooted_trees=True)[0]
	print("RF distance: ", Distance)