#!/usr/bin/python
'''
 takes all the outputted cancer networks
 builds a consensus network (separate for positive and negative edges)
 keeps track of how many networks support each edge and sums all weights

 takes as input the directory containing the sample level networks
'''
import networkx as nx
import pickle
import time
import glob
import sys, os
from pathlib import Path

base_dir = str(Path(__file__).resolve().parent.parent)

if len(sys.argv) > 1:
    net_dir = sys.argv[1]
    if not os.path.exists(net_dir):
        sys.exit('Error: Input directory containing sample-level networks does not exist.')
else:
    sys.exit('Error: Please include an input directory containing sample-level networks as a command line argument. See the README for more information.')


print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading networks...")

c_consensus_neg = nx.Graph(name="consensus_neg")
c_consensus_neg.graph['num_graphs'] = 0

c_consensus_pos = nx.Graph(name="consensus_pos")
c_consensus_pos.graph['num_graphs'] = 0

for f in glob.glob(net_dir + "/*.pickle"):
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Adding " + f + "....")
    g = pickle.load(open(f, 'rb'))
    
    g_neg = g.copy()
    g_pos = g.copy()
    
    for gi, gj in g.edges():
        if g[gi][gj]['chaos']:
            g_neg.remove_edge(gi, gj)
            g_pos.remove_edge(gi, gj)
            continue

        if g[gi][gj]['weight'] <= 0:
            g_pos.remove_edge(gi, gj)
            g_neg[gi][gj]['num_neg'] = 1
        else:
            g_neg.remove_edge(gi, gj)
            g_pos[gi][gj]['num_pos'] = 1
    
    g_neg.remove_nodes_from([n for (n, deg) in g_neg.degree() if deg == 0]) 
    g_pos.remove_nodes_from([n for (n, deg) in g_pos.degree() if deg == 0])

    comb_neg = nx.compose(g_neg, c_consensus_neg)
    neg_edge_weight = {e: c_consensus_neg.edges[e]['weight'] + g_neg.edges[e]['weight'] 
                        for e in c_consensus_neg.edges & g_neg.edges}
    neg_edge_num = {e: c_consensus_neg.edges[e]['num_neg'] + g_neg.edges[e]['num_neg'] 
                        for e in c_consensus_neg.edges & g_neg.edges}
    nx.set_edge_attributes(comb_neg, neg_edge_weight, 'weight')
    nx.set_edge_attributes(comb_neg, neg_edge_num, 'num_neg')
    c_consensus_neg = comb_neg.copy()
    c_consensus_neg.graph['num_graphs'] += 1

    comb_pos = nx.compose(g_pos, c_consensus_pos)
    pos_edge_weight = {e: c_consensus_pos.edges[e]['weight'] + g_pos.edges[e]['weight'] 
                        for e in c_consensus_pos.edges & g_pos.edges}
    pos_edge_num = {e: c_consensus_pos.edges[e]['num_pos'] + g_pos.edges[e]['num_pos'] 
                        for e in c_consensus_pos.edges & g_pos.edges}
    nx.set_edge_attributes(comb_pos, pos_edge_weight, 'weight')
    nx.set_edge_attributes(comb_pos, pos_edge_num, 'num_pos')
    c_consensus_pos = comb_pos.copy()
    c_consensus_pos.graph['num_graphs'] += 1

with open(base_dir + "/consensus_network_neg.pickle", 'wb') as out:
    pickle.dump(c_consensus_neg, out)    

with open(base_dir + "/consensus_network_pos.pickle", 'wb') as out:
    pickle.dump(c_consensus_pos, out)    
    

