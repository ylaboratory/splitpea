#!/usr/bin/python

# takes a network + outputs a file with the largest connected component
# that can be used easily with gephi

import networkx as nx

import pickle
import time

import argparse

parser = argparse.ArgumentParser(description="given pickled network output (e.g., from diff_exon2net.py), \
                                 outputs some node properties")

parser.add_argument('-i', '--in_file', required=True,
                    type=argparse.FileType('rb'),
                    help="input file w/ list of differentially expressed genes, along w/ optional sig_scores + log-fold change")
parser.add_argument('--bg', action='store_true',
                    help="input file is a background network, which means no weights")
parser.add_argument('-o', '--out_file', required=True,
                    type=argparse.FileType('w'),
                    help="output file name")

args = parser.parse_args()

print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading network....")
g = pickle.load(args.in_file)
print(g)

g_largest_cc = g.subgraph(max(nx.connected_components(g), key=len)).copy()
print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Largest connected component....")
print(g_largest_cc)

header = ["source", "target", "weight", "direction"]
args.out_file.write('\t'.join(header) + '\n')
for gi, gj in nx.edges(g_largest_cc):
    weight = g[gi][gj]['weight']
    direction = "positive"
    if weight < 0:
        direction = "negative"
    
    if g[gi][gj]['chaos']:
        direction = "chaos-" + direction

    ppi_ddi_out = [gi, gj, str(abs(weight)), direction]
    args.out_file.write('\t'.join(ppi_ddi_out) + '\n')