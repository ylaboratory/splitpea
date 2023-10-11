# exploring and visualizing results of network analyses on tcga data

library(data.table)
library(ggplot2)

largest_ccs = fread("examples/tcga_splitpea_network.largest_ccs.txt")
largest_ccs_all = largest_ccs[direction=="all", .(sample, orig_nodes, orig_edges, lcc_nodes, lcc_edges)]
names(largest_ccs_all) = c("sample", "all_orig_nodes", "all_orig_edges", "all_lcc_nodes", "all_lcc_edges")

largest_ccs = merge(largest_ccs_all, largest_ccs, by="sample")
largest_ccs[,':='(prop_orig_nodes = orig_nodes/all_orig_nodes,
                  prop_orig_edges = orig_edges/all_orig_edges,
                  prop_nodes = lcc_nodes/all_lcc_nodes,
                  prop_edges = lcc_edges/all_lcc_edges,
                  prop_lcc_nodes = lcc_nodes/orig_nodes,
                  prop_lcc_edges = lcc_edges/orig_edges)]

largest_ccs = merge(tcga_metadata, largest_ccs, by.x="entity_submitter_id", by.y="sample", all.y=T)

# > summary(largest_ccs[cancer!="background" & direction=="all"]$prop_lcc_nodes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7362  0.8025  0.8162  0.8151  0.8287  0.8594 
# 
# > summary(largest_ccs[cancer!="background" & direction=="all"]$prop_lcc_edges)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8194  0.8872  0.8966  0.8958  0.9061  0.9340 


largest_ccs[,direction:=factor(direction, levels=c("positive", "negative", "chaos", "all"))]
largest_ccs[cancer=="BRCA", cancer_name:="breast cancer"]
largest_ccs[cancer=="PAAD", cancer_name:="pancreatic cancer"]
pdf("examples/prop_pos_edges.pdf", width=6, height=4)
ggplot(largest_ccs[cancer != "background" & direction!="all"], aes(x=direction, y=prop_edges)) + 
  geom_violin(aes(fill=direction),width=2) + 
  geom_boxplot(width=0.1, color="darkgrey", alpha=0.2) +
  scale_fill_manual(values = c("positive" = "#146CAD", "negative" = "#E63846", "chaos" = "#FED283")) + 
  ylim(0,1) + 
  xlab("") + ylab("proportion of edges") +
  facet_grid(~cancer_name) + 
  guides(fill="none") +
  theme_minimal(base_size=14)
dev.off()


# consensus network stats
consensus_stats_all = fread("examples/consensus_threshold.sizes.txt")

ggplot(consensus_stats_all, aes(x=threshold, y=num_nodes)) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

ggplot(consensus_stats_all, aes(x=threshold, y=num_edges)) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

# 80% for brca so that we have things in both clusters
pdf("examples/consensus_threshold.sizes.pdf", width=6, height=4)
ggplot(consensus_stats_all[cancer=="BRCA"],
       aes(x=threshold, y=num_nodes, fill=direction, color=direction)) + 
  geom_point() + geom_line() + geom_vline(xintercept=0.8, lty='dashed') +
  scale_color_manual(values = c("positive" = "#146CAD", "negative" = "#E63846")) + 
  xlab("consensus network threshold") + 
  ylab("# nodes") +
  theme_bw(base_size=16)
dev.off()

consensus_stats_all_full = consensus_stats_all[threshold==0][,.(cancer,direction,num_nodes,num_edges)]
names(consensus_stats_all_full) = c("cancer", "direction", "all_nodes", "all_edges")
consensus_stats_all = merge(consensus_stats_all_full, consensus_stats_all, by=c("cancer", "direction"))
consensus_stats_all[,':='(prop_nodes = num_nodes/all_nodes, prop_edges = num_edges/all_edges)]
pdf("examples/consensus_threshold.prop.pdf", width=6, height=4)
ggplot(consensus_stats_all[cancer=="BRCA"],
       aes(x=threshold, y=prop_nodes, fill=direction, color=direction)) + 
  geom_point() + geom_line() + geom_vline(xintercept=0.8, lty='dashed') +
  scale_color_manual(values = c("positive" = "#146CAD", "negative" = "#E63846")) + 
  xlab("consensus network threshold") + 
  ylab("proportion of nodes") +
  theme_bw(base_size=16)
dev.off()