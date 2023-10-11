#!/bin/bash
src_dir=$(dirname "$0")
proj_dir=$(cd "$src_dir"/../ || exit; pwd)
patients_b="https://zenodo.org/record/8401618/files/BRCA-patient-rewired-networks.zip?download=1"
patients_p="https://zenodo.org/record/8401618/files/PAAD-patient-rewired-networks.zip?download=1"
consensus_b="https://zenodo.org/record/8401618/files/BRCA_consensus_networks.zip?download=1"
consensus_p="https://zenodo.org/record/8401618/files/PAAD_consensus_networks.zip?download=1"
psi_b="https://zenodo.org/record/8401618/files/BRCA-psi.zip?download=1"
psi_p="https://zenodo.org/record/8401618/files/PAAD-psi.zip?download=1"

mkdir -p $proj_dir/IRIS
wget -O $proj_dir/IRIS/BRCA-patient-rewired-networks.zip $patients_b
unzip $proj_dir/IRIS/BRCA-patient-rewired-networks.zip -d $proj_dir/IRIS

wget -O $proj_dir/IRIS/PAAD-patient-rewired-networks.zip $patients_p
unzip $proj_dir/IRIS/PAAD-patient-rewired-networks.zip -d $proj_dir/IRIS

wget -O $proj_dir/IRIS/BRCA-consensus.zip $consensus_b
unzip $proj_dir/IRIS/BRCA-consensus.zip -d $proj_dir/IRIS

wget -O $proj_dir/IRIS/PAAD-consensus.zip $consensus_p
unzip $proj_dir/IRIS/PAAD-consensus.zip -d $proj_dir/IRIS

wget -O $proj_dir/IRIS/BRCA-psi.zip $psi_b
unzip -j $proj_dir/IRIS/BRCA-psi.zip -d $proj_dir/IRIS/BRCA-psi

wget -O $proj_dir/IRIS/PAAD-psi.zip $psi_p
unzip -j $proj_dir/IRIS/PAAD-psi.zip -d $proj_dir/IRIS/PAAD-psi

# clean up weird artifacts + no longer needed zip files
rm -rf $proj_dir/IRIS/__MACOSX
rm -rf $proj_dir/IRIS/PAAD-consensus.zip
rm -rf $proj_dir/IRIS/BRCA-consensus.zip
rm -rf $proj_dir/IRIS/BRCA-patient-rewired-networks.zip
rm -rf $proj_dir/IRIS/PAAD-patient-rewired-networks.zip
rm -rf $proj_dir/IRIS/BRCA-psi.zip
rm -rf $proj_dir/IRIS/PAAD-psi.zip