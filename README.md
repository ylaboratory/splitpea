# splitpea: SPLicing InTeractions PErsonAlized

Here we present Splitpea, a method for calculating network rewiring changes due to splicing.

Splitpea takes differential exons in the form of PSI values and combines them
with domain-domain interactions (DDI) and protein-protein interactions (PPI) to
produce a PPI network with edges that change due to rewiring between two conditions. 
As a proof of principle, we show how Splitpea can be used to identify novel cancer
subtypes and driver genes affected by splicing changes using splicing changes in 
normal tissue as the background condition. 

## Citation

> Splitpea: quantifying protein interaction network rewiring changes due to alternative splicing in cancer.
Dannenfelser R and Yao V. Preprint of an article published in Pacific Symposium on Biocomputing 2024. 
[https://doi.org/10.1101/2023.09.04.556262](https://doi.org/10.1101/2023.09.04.556262)

## Code Organization

All of the code needed to run the base functionality of Splitpea is found
in the `src` directory. Helper functions for installation and exporting
for different software tools are found in the `src/utils` directory. 
Scripts need to reproduce the analysis and tables in Splitpea manuscript
can be found in `src/analysis`. 

## Setup

Conda is required for installing all relevant dependencies needed to run Splitpea.
Once conda is installed, create a new environment with the dependencies
specified in `env.yml`.

```sh
conda env create -f env.yml
```

Activate the environment:

```sh
conda activate splitpea
```

### Tabix

This project makes heavy use of tabix (v1.17), available through htslib. Though htslib is available via bioconda,
we have found that it is easier to install it separately to avoid package conflicts. We include a helper script
for installing it. This file will download htslib-1.17 to your working directory and install everything
to your home directory. The htslib-1.17 in your working directory is unnecessary after the install as binaries can be called
directly. To keep things self-contained, you could also install the binaries into this project directory
and call / modify paths as needed.

```sh
bash src/utils/install_tabix.sh
bash src/utils/tabix_feats.sh
```

### R

This project also utilizes R for preprocessing PSI background and input
PSI values. 

### Reference Files

Splitpea requires protein-protein and domain-domain interaction reference files
as well as genomic positions for exon and protein families as annotated in Pfam. While
these reference files can be manually assembled, we have provided these in the
`reference` directory. 

### Spliced Exon Data

Splitpea currently uses spliced exon data from the IRIS project. This dataset is too
large to be stored on GitHub, so only a small sample set is provided. The full set of data used is
available with the final output files on [Zenodo](https://zenodo.org/record/8401618). If you wish to use
this data with Splitpea please first download from Zenodo and use the downloaded spliced exon files
in place of the example data below.

## Running Splitpea

To illustrate how to use Splitpea, we show how it can be run to generate both
patient-specific rewired networks for individual pancreatic cancer
samples, as well a consensus network of PPI rewiring events across
all pancreatic cancer samples.

As the first step, we will to create a background level summary of splicing changes in
the normal pancreas. The following script will take alternatively spliced exon data and 
create mean summaries over all exon level coordinates across normal and tumor samples.

```sh
python src/combine_spliced_exons.py test-data
```

Next, calculate changes in PSI values (delta PSI) between each cancer sample
relative to the summarized normal pancreatic data. Here, we also empirically
calculate a p-value for each delta PSI.

```sh
Rscript src/delta_psi.R -o psis -s test-data/spliced_exons_gtex_pancreas_test_combined_mean.txt -b test-data/spliced_exons_gtex_pancreas_test.txt -t test-data/spliced_exons_tcga_paad_test.txt
```

Construct the background PPI network needed for later downstream
analysis. This network represents what happens when there are no
changes in the PPI network due to alternative splicing.

```sh
python src/get_background_ppi.py
```

Splitpea takes the preprocessed delta PSI values to generate a
network with rewired edges (as both a `.pickle` and `.dat` file). In our
case, we have a delta PSI file for each pancreatic cancer sample. Thus, we
provide Splitpea with the directory containing these delta PSI files as
well as an output directory for the final networks. 

Here we use a bash script to parallelize and run Splitpea for each
sample's psis. The `-p` flag is used to choose the number of cores
for parallelization.

```sh
# run splitpea to take a directory of sample PSI values and construct
# sample level networks of rewiring changes
bash src/run_splitpea_batch.sh -i psis -o output -p 4
```

So far, we have generated sample level networks. To get one summary network for
pancreatic cancer (the consensus network of interactions), we need to run
a script to summarize and combine across the patient samples.

```sh
# combine
python src/get_consensus_network.py output
```

## Analysis

The scripts needed to reproduce the analysis seen in the paper are in the
`src/analysis` folder. These files assume that you have run Splitpea
or have downloaded the patient specific networks for both breast cancer
(BRCA) and pancreatic cancer (PAAD), as well as the consensus networks
for each tumor type. The easiest setup is to download 
`BRCA-patient-rewired-networks.zip`, `BRCA_consensus_networks.zip`,
`PAAD-patient-rewired-networks.zip`, and  `PAAD_consensus_networks.zip`
from [Zenodo](https://zenodo.org/record/8401618)
and place them in a directory called `IRIS` in the main
Splitpea folder. Once these files are downloaded or created you can generate the
the main figures using the following scripts:

- `clustering_psi.R`: contains the code for building the heatmaps in Figure 2
- `something`: build the gains / and losses in Figure 3
- `calc_network_stats.py`: stats for the individual networks
- `get_consensus_network_stats.py`: stats for the consensus networks
- `get_largestcc_sizes.py`:
- `clustering_networks.ipynb`: generates network embeddings using FEATHER
- `clustering_networks.R`: further analysis of the network embeddings (Figure 6)
