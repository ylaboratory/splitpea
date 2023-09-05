# splitpea: SPLicing InTeractions PErsonAlized

Here we present, Splitpea, a method for calculating network rewiring changes due to splicing.

Splitpea takes differential exons in the form of PSI values and combines them
with domain-domain interactions (DDI) and protein-protein interactions (PPI) to
produce a PPI network with edges that change due to rewiring between two conditions. 
As a proof of principle, we show how Splitpea can be used to identify novel cancer
subtypes and driver genes affected by splicing changes. 

Differential exon input in the form of PSI values is currently the only supported input.

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

This project makes heavy use of tabix (v1.17), available through htslib. Though htslib  is available via bioconda,
we have found that it is easier to install it separately to avoid package conflicts. We include a helper script
for installing it. This file will download htslib-1.17 to your working directory and install everything
to your home directory. The htslib-1.17 in your working directory is unnecessary after the install as binaries can be called
directly. To keep things self-contained, you could also install the binaries into this project directory
and call / modify paths as needed.

```sh
python src/utils/install_tabix.sh
```

## Running Splitpea

**Update** The main hook to run Splint is currently `diff_exon2.net.py`, which takes as
input a list of exons and several preprocessed files / directories. 
General information about how the references files are structured
can be found in the `Reference Information` section below.

**Coming Soon** Supporting files can be downloaded from [Zenodo](https://zenodo.org) 


```sh
# DexSeq input - human
python src/diff_exon2net.py -i examples/bodymap_tissue.dexseq.txt -o out/bodymap_tissue.dexseq.diff_isoform -s 0.05 --org human -d reference/Homo_sapiens.GRCh38.91.gff -p 2 --skip 1 -f dexseq -v

```

## Code Organization

**Update** Most of the code is located in `src/` (though sometimes, code relevant to
a specific part can be found in the corresponding subdirectory...).

The main hook to run splint is currently `diff_exon2.net.py`, which takes as
input a list of exons and several preprocessed files / directories.
At the time I originally wrote this up, I left a note for myself that I
found a bug in tabix (pytabix and tabix command line tool may have had
different indexing - 0-based vs 1-based), not sure if it's fixed in more
recent implementations.


## Reference Information

**Update**

`interactions`:
Split by type of interaction re any type of pairwise interaction, the big one
being ppi, but also includes several other types. Within each category, there
is usually some kind of `gold_standard` directory that has processed edge
information from various databases (often w/ additional organization by organism).

The `orthology` directories have pooled and transferred edges from all other
organisms for which there was gold standard.

`mappings`:
Various gene identifier mappings across different organisms.

`reference`:
Mostly information re gene + domain coordinates processed for each organism.
