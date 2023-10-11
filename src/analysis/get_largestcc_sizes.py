#!/usr/bin/python

# takes a network + does some simple node property calculations 
# for the largest connected component

import networkx as nx

import pickle
import time

from pathlib import Path
import glob

base_dir = str(Path(__file__).resolve().parent.parent)
net_dir = base_dir + "/output/"
cancers = ['BRCA', 'PAAD']

out = open(base_dir + "/results/tcga_splitpea_network.largest_ccs.txt", 'w')
out.write('\t'.join(["cancer", "sample", "direction", "orig_nodes", "orig_edges", "lcc_nodes", "lcc_edges"]) + '\n')

g = pickle.load(open(base_dir + "/reference/human_ppi_ddi_bg.pickle", 'rb'))
g_largest_cc = g.subgraph(max(nx.connected_components(g), key=len)).copy()
out.write('\t'.join(["background",
                    "ppi_ddi_bg",
                    "all",
                    str(len(g)),
                    str(g.number_of_edges()),
                    str(len(g_largest_cc)),
                    str(g_largest_cc.number_of_edges())]) + '\n')

for c in cancers:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading " + c + " networks...")

    for f in glob.glob(net_dir + c + "/mean/*.pickle"):
        file_pref = f.split('/')[-1].split('.')[0]
        g = pickle.load(open(f, 'rb'))
        
        g_largest_cc = g.subgraph(max(nx.connected_components(g), key=len)).copy()

        out.write('\t'.join([c,
                             file_pref,
                             "all",
                             str(len(g)),
                             str(g.number_of_edges()),
                             str(len(g_largest_cc)),
                             str(g_largest_cc.number_of_edges())]) + '\n')
    
        g_largest_cc_neg = g_largest_cc.copy()
        g_largest_cc_pos = g_largest_cc.copy()
        g_largest_cc_chaos = g_largest_cc.copy()
        
        for gi, gj in g_largest_cc.edges():
            if g_largest_cc[gi][gj]['chaos']:
                g_largest_cc_neg.remove_edge(gi, gj)
                g_largest_cc_pos.remove_edge(gi, gj)
                continue
            else:
                g_largest_cc_chaos.remove_edge(gi, gj)

            if g_largest_cc_pos[gi][gj]['weight'] <= 0:
                g_largest_cc_pos.remove_edge(gi, gj)
            else:
                g_largest_cc_pos[gi][gj]['distance'] = 1 - g_largest_cc_pos[gi][gj]['weight']
                g_largest_cc_neg.remove_edge(gi, gj)
        
        # pos and neg may not be fully connected anymore
        g_largest_cc_neg.remove_nodes_from([n for (n, deg) in g_largest_cc_neg.degree() if deg == 0]) 
        g_largest_cc_pos.remove_nodes_from([n for (n, deg) in g_largest_cc_pos.degree() if deg == 0])
        g_largest_cc_chaos.remove_nodes_from([n for (n, deg) in g_largest_cc_chaos.degree() if deg == 0]) 

        out.write('\t'.join([c,
                             file_pref,
                             "negative",
                             str(len(g_largest_cc)),
                             str(g_largest_cc.number_of_edges()),
                             str(len(g_largest_cc_neg)),
                             str(g_largest_cc_neg.number_of_edges())]) + '\n')
    
        out.write('\t'.join([c,
                             file_pref,
                             "positive",
                             str(len(g_largest_cc)),
                             str(g_largest_cc.number_of_edges()),
                             str(len(g_largest_cc_pos)),
                             str(g_largest_cc_pos.number_of_edges())]) + '\n')

        out.write('\t'.join([c,
                             file_pref,
                             "chaos",
                             str(len(g_largest_cc)),
                             str(g_largest_cc.number_of_edges()),
                             str(len(g_largest_cc_chaos)),
                             str(g_largest_cc_chaos.number_of_edges())]) + '\n')
        
out.close()