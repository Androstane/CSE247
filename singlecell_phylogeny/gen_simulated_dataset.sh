total_rep=5
codedir=~/github/SingleCellCNABenchmark
ref=/projects/nakhleh/xf2/reference/hg19.fa
tree=from_first_step.tree.npy
samplesegcopyf=segcopy
#/storage/hpc/work/nakhleh/xf2/benchmark/sim_ploidy/p1/rep_2/ginkgo/SegCopy
for i in `seq 1 5`; do
    mkdir p$i
done

# ploidy 1.55 
ploidy=1
for rep in `seq 1 $total_rep`; do
    dir=p$ploidy/rep$rep
    # generate the tree
    python $codedir/main.par.py -S $codedir/wgsim-master -r $dir -n 100 -p 1 -X 25 -t $ref -W 0 -l 36 -m 2000000 -e 5000000 -d 1 -c 3
done

# ploidy 2.1
ploidy=2
for rep in `seq 1 $total_rep`; do
    dir=p$ploidy/rep$rep
    # generate the tree
    python $codedir/main.par.py -S $codedir/wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.05 -l 36 -m 2000000 -e 5000000
done

# ploidy 3
ploidy=3
for rep in `seq 1 $total_rep`; do
    dir=p$ploidy/rep$rep
    # generate the tree
    python $codedir/main.par.py -S $codedir/wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.5 -l 36 -m 2000000 -e 5000000 -E 1 
done

# ploidy 3.8
ploidy=4
for rep in `seq 1 $total_rep`; do
    dir=p$ploidy/rep$rep
    # generate the tree
    python $codedir/main.par.py -S $codedir/wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.9 -l 36 -e 5000000 -E 1 -m 10000000 
done

# ploidy 5.26
ploidy=5
for rep in `seq 1 $total_rep`; do
    dir=p$ploidy/rep$rep
    # generate the tree
    python $codedir/main.par.py -S $codedir/wgsim-master -r $dir -n 100 -p 1 -X 8 -t $ref -W 1 -C 0.9 -l 36 -e 5000000 -E 1 -m 10000000 -J 0.55 
done

for ploidy in `seq 1 5`; do
    for rep in `seq 1 $total_rep`; do
        dir=p$ploidy/rep$rep
        pushd $dir
        # generate gt.all.csv
        python $codedir/read_tree.py -s -f $tree > gt.all.csv
        # generate the leaves.txt
        python $codedir/read_tree.py -l $tree > leaves.txt 
        # copy sample segcopy file here
        cp $samplesegcopyf ./ 
        popd
    done
done 
