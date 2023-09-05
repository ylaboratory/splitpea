# splitpea: SPLicing InTeractions PErsonAlized

Here we present, Splitpea, a method for calculating network rewiring changes due to splicing.

Splitpea takes differential exons in the form of PSI values and combines them
with domain-domain interactions (DDI) and protein-protein interactions (PPI) to
produce a PPI network with edges that change due to rewiring between two conditions. 
As a proof of principle, we show how Splitpea can be used to identify novel cancer
subtypes and driver genes affected by splicing changes using splicing changes in 
normal tissue as the background condition. 

## Citation

> Splitpea: uncovering cancer patient-specific protein interaction network rewiring.
Dannenfelser R and Yao V. Preprint forthcoming via BioRxiv. 2023.

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
python src/utils/install_tabix.sh
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
large to be stored on GitHub so only a small sample set is provided. The full set of data used is
available in several processed forms via [Zenodo](https://zenodo.org). If you wish to use
this data with Splitpea please first download from Zenodo and upload as a replacement to the sample
data in the `examples` splitpea directory.

## Running Splitpea

To illustrate how to use Splitpea we show how it can be run to generate both
patient-specific rewired networks for individual pancreatic cancer
samples, as well a consensus network of PPI rewiring events across
all pancreatic cancer samples.


As the first step, we will to create a background level summary of splicing changes in
the normal pancreas. The following script will
take alternatively spliced exon data and create mean and median summaries
over all exon level coordinates across normal and tumor samples.

```sh
python src/combine_spliced_exons.py
```

Next, calculate changes in PSI values (delta PSI) between each cancer sample
relative to the summarized normal pancreas data. Here we also empirically
calculate a p-value for each delta PSI.

```sh
RScript delta_psi.R
```

Construct the background PPI network needed for later downstream
analysis. This network represents what happens when there are no
changes in the PPI network due to alternative splicing.

```sh
python src/get_background_ppi.py
```

Splitpea takes the preprocessed delta PSI values to generate a
network with rewired edges (as both a `.pickle` and `.dat` file). In our
case, we have a delta PSI file for each pancreatic cancer sample. Thus we
provide splitpea with the directory containing these delta PSI files as
well as an output directory for the final networks. 

```sh
# run splitpea to take a directory of sample PSI values and construct
# sample level networks of rewiring changes
python src/splitpea.py -i dir_of_sample_psi_files --skip 1 -o out_file_prefix -v
```

So far we have generated sample level networks, to get one summary network for
pancreatic cancer (the consensus network of interactions) we need to run
a script to summarize and combine across the patient samples.

```sh
# combine
python src/get_consensus_network.py
```

## Code Organization

All of the code needed to run the base functionality of Splitpea is found
in the `src` directory. Helper functions for installation and exporting
for different software tools are found in the `src/utils` directory. 
Scripts need to reproduce the analysis and tables in Splitpea manuscript
can be found in `src/analysis`.

