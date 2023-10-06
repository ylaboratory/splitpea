#!/bin/bash
src_dir=$(dirname "$0")
proj_dir=$(cd "$src_dir"/../../ || exit; pwd)

echo "tabixing..."
tabix -b 7 -e 8 -s 6 $proj_dir/reference/human_pfam_genome_coords_sorted.txt.gz