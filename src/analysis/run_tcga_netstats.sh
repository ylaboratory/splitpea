#!/bin/bash

src_dir=$(dirname "$0")
proj_dir=$(cd "$src_dir"/../ || exit; pwd)

data_dir="$proj_dir"/examples/
cancer_types=(BRCA PAAD)

net_dir="$proj_dir"/output

parallel=40

for canc in "${cancer_types[@]}"
do
    echo $canc
    canc_dir="$data_dir"/"$canc"/mean
    net_in_dir="$net_dir"/"$canc"/mean
    for i in `ls "$canc_dir"`
    do
        ipref=${i/-psi\.txt/}
        (
            echo "$proj_dir"/src/calc_network_stats.py -i "$net_in_dir"/"$ipref".edges.pickle -o "$net_in_dir"/"$ipref".centralities.txt 
            python "$proj_dir"/src/calc_network_stats.py -i "$net_in_dir"/"$ipref".edges.pickle -o "$net_in_dir"/"$ipref".centralities.txt
        ) &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge $parallel ]
        then
            wait
            NPROC=0
        fi
    done
done

echo "$proj_dir"/src/calc_network_stats.py -i "$proj_dir"/reference/human_ppi_ddi_bg.pickle -o "$proj_dir"/reference/human_ppi_ddi_bg.centralities.txt --bg
python "$proj_dir"/src/calc_network_stats.py -i "$proj_dir"/reference/human_ppi_ddi_bg.pickle -o "$proj_dir"/reference/human_ppi_ddi_bg.centralities.txt --bg