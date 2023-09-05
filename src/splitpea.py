#!/usr/bin/python

'''
takes exon info + delta_psi + pval and 
uses pre-calculated exon - protein domain info, ppi, ddi
to output network w/ scores
'''


import argparse
import pickle

from collections import defaultdict
import networkx as nx
import os
from numpy import average
from functools import partial
from os.path import basename
from subprocess import Popen, PIPE

from exons import Exons
from pathlib import Path


import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                    datefmt="%m-%d-%Y %I:%M:%S %p",
                    level=logging.ERROR)
logger = logging.getLogger('diff_exon')

# file path constants
base_dir = str(Path(__file__).resolve().parent.parent)
ppif = base_dir + "/reference/human_ppi_0.5.dat"
ddif = base_dir + "/reference/ddi_0.5.dat"
entrezpfamf = base_dir + "/reference/human_entrez_pfam.txt"
pfamcoordsf = base_dir + "/reference/human_pfam_genome_coords_sorted.txt.gz"
tbf = base_dir + "/reference/human_pfam_genome_coords_sorted.txt.gz"


def tb_query(tb_file, chrom, start, end):
    '''
    Call tabix and generate an array of strings for each line it returns.
    adapted from https://github.com/slowkow/pytabix
    '''
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', tb_file, query], stdout=PIPE, text=True)
    for line in process.stdout:
        yield line.strip().split()

def calculate_edges(g, dex, tb, all_ddis, all_ppis, entrez_pfams):
    '''
    given a graph g, containing nodes corresponding to
    differentially expressed genes in dex, build a network
    with edges from the reference ppi and ddi

    importantly, within each domain for each gene, minimum dpsi is tabulated
    so we know as a whole, conservatively, that domain is net loss / net gain

    nothing is returned but g is updated
    '''
    for ex in dex.dexon2scores:
        ex_dpsi = dex.dexon2scores[ex][0]

        for r in tb_query(tb, ex[0], ex[1], ex[2]):
            if r[8] != ex[3]:  # check for same strand
                continue

            gi = r[0]
            pfi = r[2]
            pfi_start = r[3]
            pfi_end = r[4]

            if gi not in all_ppis or gi not in entrez_pfams: # need to have ppi and pfam info for gi
                continue

            if pfi not in all_ddis or pfi not in entrez_pfams: # need to have ddi and genes with pfi
                continue

            for gj in all_ppis[gi]:
                if gj not in entrez_pfams:
                    continue
                
                ddi_overlap = set(entrez_pfams[gj]).intersection(set(all_ddis[pfi]))
                if len(ddi_overlap) > 0:
                    if gi not in g:
                        g.add_node(gi)

                    if pfi not in g.nodes[gi]:
                        g.nodes[gi][pfi] = defaultdict(dict)
                    
                    if (pfi_start, pfi_end) not in g.nodes[gi][pfi]:
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['dexons'] = set()
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi'] = 1

                    g.nodes[gi][pfi][(pfi_start, pfi_end)]['dexons'].add((ex, ex_dpsi))
                    if ex_dpsi < g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi']:
                        g.nodes[gi][pfi][(pfi_start, pfi_end)]['min_dpsi'] = ex_dpsi

                    if gj not in g[gi]:
                        g.add_edge(gi, gj,
                                    ppi_weight = all_ppis[gi][gj]['weight'],
                                    ppi_ddis = defaultdict(partial(defaultdict, set))) # cannot use lambda if want to pickle later
                    
                    g.edges[gi,gj]['ppi_ddis'][gi][pfi] |= ddi_overlap


def deduce_final_edge_weights(g):
    '''
    final edge weight / directionality is "collapsed" based on calculate_edges:
    in the event that there are multiple domains, if the directionality is consistent,
    the overall direction is clear

    if there are conflicts, then directionality is "chaotic" - final weight is going to be
    difference between the avg of positive dpsis and avg of negative dpsis
    '''
    num_chaos = 0
    for gi, gj in g.edges():
        chaos = False
        pos_dpsi = list()
        neg_dpsi = list()

        for g_ddi in g[gi][gj]['ppi_ddis']: # maximum length 2 if both genes had diff exon
            for pf_g in g[gi][gj]['ppi_ddis'][g_ddi]: # all these loops are very shallow
                for pf_loc in g.nodes[g_ddi][pf_g]:
                    pf_dpsi = g.nodes[g_ddi][pf_g][pf_loc]['min_dpsi']
                    if pf_dpsi < 0: # with dpsi cutoff in place, will not have 0s
                        neg_dpsi.append(pf_dpsi)
                    else:
                        pos_dpsi.append(pf_dpsi)

        if len(pos_dpsi) > 0 and len(neg_dpsi) > 0:
            chaos = True
            num_chaos += 1

        g[gi][gj]['weight'] = average(pos_dpsi + neg_dpsi)
        g[gi][gj]['chaos'] = chaos
        g[gi][gj]['num_dpsi_pfs'] = len(pos_dpsi) + len(neg_dpsi)

    g.graph['num_chaos'] = num_chaos

def fileExists(f):
    '''
    check if file exists
    '''
    if not os.path.isfile(f):
        print("ERROR: the file " + f + " was not found.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="given a list of differentially \
        expressed exons, outputs splicing-specific network")
    parser.add_argument('-i', '--in_file', required=True, type=str,
                        help="list of differentially expressed exons, along w/ optional sig_scores \
                        + log-fold change")
    parser.add_argument('--skip', type=int, default=0, help="number of lines to skip in input file (default=0)")
    parser.add_argument('--dpsi_cut', type=float, default=0.05, \
                        help="absolute delta psi cutoff for differential exons to use in network (default=0.05)")
    parser.add_argument('--sigscore_cut', type=float, default=0.05,
                        help="pvalue / qvalue / fdr cutoff for significance \
                        of differential exons to use in network (default=0.05)")
    parser.add_argument('--include_nas', type=bool, default=True,
                        help="include NAs in pvalue estimates for differential exon selection (default=True)")
    parser.add_argument('-o', '--out_file_prefix', required=True, type=str,
                        help="prefix for output files")
    parser.add_argument('-v', '--verbose', action="store_true", help="verbosity")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    # check all the files
    fileExists(ppif)
    fileExists(entrezpfamf)
    fileExists(pfamcoordsf)
    fileExists(pfamcoordsf + ".tbi")

    # load the reference files
    logger.info("Loading DDIs....")
    all_ddis = nx.read_weighted_edgelist(ddif)
    logger.info("# pfams: %i, # ddis: %i", len(all_ddis), all_ddis.number_of_edges())

    logger.info("Loading PPIs....")
    all_ppis = nx.read_weighted_edgelist(ppif) 
    logger.info("# proteins: %i, # interactions: %i", len(all_ppis), all_ppis.number_of_edges())

    logger.info("Loading gene-protein domain info....")
    # actually a bipartite graph, but instead of labeling the nodes as classes
    # we will just read it in as a general undirected graph for convenience
    entrez_pfams = nx.bipartite.read_edgelist(entrezpfamf)

    # reading in differential exon results
    logger.info("Reading differential exon results....")
    dexons = Exons()
    dexons.read_dexons(args.in_file, args.skip, args.dpsi_cut, args.sigscore_cut, args.include_nas)

    logger.info("# of differentially expressed exons: %i", len(dexons))

    logger.info("Calculating network....")
    # constructing graph w/ ppi nodes incident on genes w/ diff-exp'ed exons + ddis
    diff_splice_g = nx.Graph(name=basename(args.in_file))
    calculate_edges(diff_splice_g, dexons, tbf, all_ddis, all_ppis, entrez_pfams)
    deduce_final_edge_weights(diff_splice_g)

    logger.info("# nodes: %i, # edges: %i", len(diff_splice_g), diff_splice_g.number_of_edges())

    logger.info("Outputting as pickle...")
    with open(args.out_file_prefix + '.edges.pickle', 'wb') as out:
        pickle.dump(diff_splice_g, out)
   
    with open(args.out_file_prefix + '.edges.dat', 'w') as out:
        header = ["node1", "node2", "weight", "chaos"]
        out.write('\t'.join(header) + '\n')
        for gi, gj in nx.edges(diff_splice_g):
            ppi_ddi_out = [gi, gj, str(diff_splice_g[gi][gj]['weight']),
                           str(diff_splice_g[gi][gj]['chaos'])]
            out.write('\t'.join(ppi_ddi_out) + '\n')

    logger.info("Done")
