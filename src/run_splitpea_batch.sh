#!/bin/bash
src_dir=$(dirname "$0")
proj_dir=$(cd "$src_dir"/../ || exit; pwd)
output_dir=$proj_dir"/nets-out"
parallel=10

usage() {
  echo "$(basename $0)
  Runs the Splitpea algorithm on a folder of delta_psi values.

  -h, --help              prints instructions

  Options: 
  -i <folder path>             path to folder containing delta psi values
  -o <folder path>             path to folder for output networks
  -p <int>                     number of cores for the run (default: 10)
  "
}

while getopts :i:o:p: option; do
  case "${option}" in
    i) data_dir=${OPTARG};;
    o) output_dir=${OPTARG};;
    p) para=${OPTARG};;
    :) echo "Error: -${OPTARG} requires an argument."; usage; exit -1;;
    *) usage; exit -1;;
  esac
done

# defaults
if [ -n "$para" ]; then parallel=$para; fi
mkdir -p "$output_dir"

for i in `ls "$data_dir"`
do
    ipref=${i/-psi\.txt/}
    (
        echo "$proj_dir"/src/splitpea.py -i "$data_dir"/$i --skip 1 -o "$output_dir"/$ipref -v
        python "$proj_dir"/src/splitpea.py -i "$data_dir"/$i --skip 1 -o "$output_dir"/$ipref -v
    ) &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge $parallel ]
    then
        wait
        NPROC=0
    fi
done
