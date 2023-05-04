#!/bin/bash


log_dir="/home/nima/scratch/dual_rna/count"

for exp in ~/scratch/dual_rna/raw_reads/experiment*/*
do
    exp_type=$(basename $exp)
    exp_id=$(basename $(dirname $exp))

    mkdir -p $log_dir/$exp_id/$exp_type
    for sra_dir in ~/scratch/dual_rna/raw_reads/$exp_id/$exp_type/*
    do
        srid=$(basename $sra_dir)
        if [ -d $sra_dir/star ]
        then
            # echo $exp_id
            cp $sra_dir/star/star_ReadsPerGene.out.tab $log_dir/$exp_id/$exp_type/$srid.tab
        fi
    done
done
