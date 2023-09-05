#!/usr/bin/python

# takes consensus networks calculated from `get_consensus_network.py`
# calculates size of largest connected component depending on % support needed

import networkx as nx

import pickle
import time

from pathlib import Path

base_dir = str(Path(__file__).resolve().parent.parent)
cancers = ['BRCA', 'PAAD']
# thresholds = [0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1]
thresholds = [0.8]

# out = open(base_dir + "/results/consensus_threshold.lcc_sizes.txt", 'w')
# out.write('\t'.join(['cancer', 'direction', 'threshold', 'num_nodes', 'num_edges']) + '\n')
for c in cancers:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading " + c + " consensus networks...")

    c_consensus_neg = pickle.load(open(base_dir + "/results/" + c + "_consensus_neg.pickle", 'rb'))
    c_consensus_pos = pickle.load(open(base_dir + "/results/" + c + "_consensus_pos.pickle", 'rb'))

#     # cneg_largest_cc = c_consensus_neg.copy()
#     # cpos_largest_cc = c_consensus_pos.copy() 
    cneg_largest_cc = c_consensus_neg.subgraph(max(nx.connected_components(c_consensus_neg), key=len)).copy()
    cpos_largest_cc = c_consensus_pos.subgraph(max(nx.connected_components(c_consensus_pos), key=len)).copy()

#     out.write('\t'.join([c, 'negative', '0', str(len(cneg_largest_cc)), str(cneg_largest_cc.number_of_edges())]) + '\n')
#     out.write('\t'.join([c, 'positive', '0', str(len(cpos_largest_cc)), str(cpos_largest_cc.number_of_edges())]) + '\n')

    for t in thresholds:
        thres_neg = cneg_largest_cc.copy()
        for e in cneg_largest_cc.edges:
            if cneg_largest_cc.edges[e]['num_neg'] < t * cneg_largest_cc.graph['num_graphs']:
                thres_neg.remove_edge(e[0], e[1])

        thres_neg.remove_nodes_from([n for (n, deg) in thres_neg.degree() if deg == 0]) 
        if not nx.is_empty(thres_neg):
            thres_neg = thres_neg.subgraph(max(nx.connected_components(thres_neg), key=len)).copy()

        thres_pos = cpos_largest_cc.copy()
        for e in cpos_largest_cc.edges:
            if cpos_largest_cc.edges[e]['num_pos'] < t * cpos_largest_cc.graph['num_graphs']:
                thres_pos.remove_edge(e[0], e[1])

        thres_pos.remove_nodes_from([n for (n, deg) in thres_pos.degree() if deg == 0]) 
        if not nx.is_empty(thres_pos):
            thres_pos = thres_pos.subgraph(max(nx.connected_components(thres_pos), key=len)).copy()

        with open(base_dir + "/results/" + c + "_consensus_neg.thres_" + str(t) + ".pickle", 'wb') as out:
            print(thres_neg)
            pickle.dump(thres_neg, out)
        
        with open(base_dir + "/results/" + c + "_consensus_pos.thres_" + str(t) + ".pickle", 'wb') as out:
            print(thres_pos)
            pickle.dump(thres_pos, out)

        with open(base_dir + "/results/" + c + "_consensus_joint.thres_" + str(t) + ".pickle", 'wb') as out:
            print(nx.compose(thres_neg, thres_pos))
            pickle.dump(nx.compose(thres_neg, thres_pos), out) 

#         out.write('\t'.join([c, 'negative', str(t), str(len(thres_neg)), str(thres_neg.number_of_edges())]) + '\n')
#         out.write('\t'.join([c, 'positive', str(t), str(len(thres_pos)), str(thres_pos.number_of_edges())]) + '\n')