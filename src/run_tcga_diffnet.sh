#!/bin/bash

src_dir=$(dirname "$0")
proj_dir=$(cd "$src_dir"/../ || exit; pwd)

data_dir="$proj_dir"/examples/
cancer_types=(BRCA PAAD)

output_dir="$proj_dir"/output

parallel=40

mkdir -p "$output_dir"
for canc in "${cancer_types[@]}"
do
    echo $canc
    canc_in_dir="$data_dir"/"$canc"/mean
    canc_out_dir="$output_dir"/"$canc"/mean
    mkdir -p "$canc_out_dir"
    for i in `ls "$canc_in_dir"`
    do
        ipref=${i/-psi\.txt/}
        (
            echo "$proj_dir"/src/splitpea.py -i "$canc_in_dir"/$i --skip 1 -o "$canc_out_dir"/$ipref -v
            python "$proj_dir"/src/splitpea.py -i "$canc_in_dir"/$i --skip 1 -o "$canc_out_dir"/$ipref -v
        ) &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge $parallel ]
        then
            wait
            NPROC=0
        fi
    done
done