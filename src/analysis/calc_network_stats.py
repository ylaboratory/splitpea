#!/usr/bin/python

# takes a network + does some simple node property calculations 
# for the largest connected component

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

print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Calculating node stats....")
g_degree = nx.degree(g_largest_cc)
if not args.bg:
    g_degree_weighted = nx.degree(g_largest_cc, weight="weight")

    g_largest_cc_neg = g_largest_cc.copy()
    g_largest_cc_pos = g_largest_cc.copy()
    
    for gi, gj in g_largest_cc.edges():
        if g_largest_cc[gi][gj]['chaos']:
            g_largest_cc_neg.remove_edge(gi, gj)
            g_largest_cc_pos.remove_edge(gi, gj)
            continue

        if g_largest_cc_pos[gi][gj]['weight'] <= 0:
            g_largest_cc_pos.remove_edge(gi, gj)
        else:
            g_largest_cc_pos[gi][gj]['distance'] = 1 - g_largest_cc_pos[gi][gj]['weight']
            g_largest_cc_neg.remove_edge(gi, gj)

    # pos and neg may not be fully connected anymore
    g_largest_cc_neg = g_largest_cc_neg.subgraph(max(nx.connected_components(g_largest_cc_neg), key=len))
    g_largest_cc_pos = g_largest_cc_pos.subgraph(max(nx.connected_components(g_largest_cc_pos), key=len))

    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Largest connected component (neg)....")
    print(g_largest_cc_neg)
    g_neg_degree = nx.degree(g_largest_cc_neg)
    g_neg_degree_weighted = nx.degree(g_largest_cc_neg, weight="weight")

    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Largest connected component (pos)....")
    print(g_largest_cc_pos)
    g_pos_degree = nx.degree(g_largest_cc_pos)
    g_pos_degree_weighted = nx.degree(g_largest_cc_pos, weight="weight")
    g_between_cent = nx.betweenness_centrality(g_largest_cc_pos, weight="distance")

    args.out_file.write('\t'.join(["entrez",
                                   "degree",
                                   "weighted_degree",
                                   "neg_degree",
                                   "neg_weighted_degree",
                                   "pos_degree",
                                   "pos_weighted_degree",
                                   "pos_betweenness"]) + '\n')
else:
    g_between_cent = nx.betweenness_centrality(g_largest_cc)
    args.out_file.write('\t'.join(["entrez",
                                   "degree",
                                   "betweenness"]) + '\n')

for entrez in g_largest_cc.nodes:
    if not args.bg:
        g_neg = 0
        g_negw = 0

        g_pos = 0
        g_posw = 0
        g_bcent = 0

        if entrez in g_largest_cc_neg:
            g_neg = g_neg_degree[entrez]
            g_negw = g_neg_degree_weighted[entrez]

        if entrez in g_largest_cc_pos:
            g_pos = g_pos_degree[entrez]
            g_posw = g_pos_degree_weighted[entrez]
            g_bcent = g_between_cent[entrez]
        
        args.out_file.write('\t'.join([entrez, 
                                    str(g_degree[entrez]),
                                    str(g_degree_weighted[entrez]),
                                    str(g_neg),
                                    str(g_negw),
                                    str(g_pos),
                                    str(g_posw),
                                    str(g_bcent)]) + '\n')
    else:
        args.out_file.write('\t'.join([entrez,
                                       str(g_degree[entrez]),
                                       str(g_between_cent[entrez])]) + '\n')
