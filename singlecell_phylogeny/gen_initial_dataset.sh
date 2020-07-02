ploidy=$1
rep=$2
codedir=~/github/SingleCellCNABenchmark
dir=p$1/rep$2
treenpy=$dir/from_first_step.tree.npy
gtallf=$dir/gt.all.csv
samplesegcopyf=$dir/segcopy
leafid=$dir/leaves.txt
# output the ground truth pairwise distance
python $codedir/read_tree.py -d -f $treenpy > $dir/gt.pairwisedist.txt 
# output the ground truth tree in newick format
python $codedir/gen_newick.py $treenpy > $dir/gt.newick
# output the ground truth copy number profile 
python $codedir/comparison/bin_groundtruth_woSegCopy.py -a $samplesegcopyf -b $gtallf --leafonly -l $leafid > $dir/gt.cnp

for i in `seq 1 5`; do for j in `seq 1 5`; do pushd p$i/rep${j}; python ~/github/SingleCellCNABenchmark/comparison/bin_groundtruth_woSegCopy.py -a ../../segcopy -b gt.all.csv > gt.all.cnp; python ~/github/SingleCellCNABenchmark/read_tree.py -O -F gt.all.cnp -f from_first_step.tree.npy > newoverlappingCNA.csv; popd; done; done
