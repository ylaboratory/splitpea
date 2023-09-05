# exploring and visualizing results of network analyses on tcga data

library(data.table)
library(ggplot2)
library(clusterProfiler)

tcga_metadata = fread("/grain/rad4/splitpea/examples/TCGA/TCGA_metadata.csv")

tcga_metadata = tcga_metadata[, c('entity_submitter_id', 'sample_types', 'center', 'ajcc_pathologic_stage',
                          'tissue_or_organ_of_origin', 'days_to_last_follow_up', 'days_to_death', 'primary_diagnosis',
                          'year_of_diagnosis', 'vital_status', 'race', 'gender', 'ethnicity', 'age_at_index',
                          'site_of_resection_or_biopsy')]

cosmic_tier1 = read.gmt("/grain/rad4/gmt_files/COSMIC-driver-tumor-type-tier1only.gmt")


largest_ccs = fread("results/tcga_splitpea_network.largest_ccs.txt")
largest_ccs[,':='(prop_nodes = lcc_nodes/orig_nodes, prop_edges = lcc_edges/orig_edges)]

largest_ccs = merge(tcga_metadata, largest_ccs, by.x="entity_submitter_id", by.y="sample", all.y=T)

# > summary(largest_ccs[cancer!="background" & direction=="all"]$prop_nodes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7362  0.8025  0.8163  0.8152  0.8290  0.8594 
# 
# > summary(largest_ccs[cancer!="background" & direction=="all"]$prop_edges)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8194  0.8872  0.8969  0.8958  0.9061  0.9340 



ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=cancer, y=prop_nodes)) + facet_grid(~direction)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=cancer, y=prop_edges)) + facet_grid(~direction)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=cancer, y=lcc_nodes)) + facet_grid(~direction)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=cancer, y=lcc_edges)) + facet_grid(~direction)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=direction, y=prop_nodes)) + facet_grid(cancer~ajcc_pathologic_stage)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=direction, y=lcc_edges)) + facet_grid(cancer~ajcc_pathologic_stage)

ggplot(largest_ccs[cancer!="background"]) + 
  geom_boxplot(aes(x=ajcc_pathologic_stage, y=prop_edges)) + facet_grid(cancer~direction)

ggplot(largest_ccs[cancer=="PAAD"]) + 
  geom_boxplot(aes(x=direction, y=prop_edges)) + facet_wrap(~primary_diagnosis)

ggplot(largest_ccs[cancer=="PAAD" & direction=="neg"]) + 
  geom_point(aes(x=lcc_edges, y=days_to_death))


largest_ccs[cancer=="BRCA", cancer_name:="breast cancer"]
largest_ccs[cancer=="PAAD", cancer_name:="pancreatic cancer"]
pdf("results/prop_pos_edges.pdf", width=4, height=3)
ggplot(largest_ccs[cancer != "background" & direction!="all"], aes(x=direction, y=prop_edges)) + 
  geom_violin(aes(fill=direction)) + 
  geom_boxplot(width=0.1, color="darkgrey", alpha=0.2) +
  scale_fill_manual(values = c("positive" = "#146CAD", "negative" = "#E63846", "chaos" = "#FED283")) + 
  ylim(0,1) + 
  xlab("") + ylab("proportion of edges") +
  facet_grid(~cancer_name) + 
  guides(fill="none") +
  theme_minimal(base_size=14)
dev.off()


hsa_gene_mappings = fread("/grain/resources/gene_mappings/output/current/hsa_mapping_all.txt")
hsa_symbol_entrez = unique(hsa_gene_mappings[,.(symbol, entrez)])

bg_centralities = fread("reference/human_ppi_ddi_bg.centralities.txt")
setnames(bg_centralities, c("entrez", "bg_degree", "bg_betweenness"))

canc_dir = "output/PAAD/mean/"
paad_centralities = rbindlist(lapply(list.files(canc_dir, pattern="*centralities.txt"),
                                     function(x) data.table(tcga_samp_id = gsub(".centralities.txt", "", x),
                                                            fread(paste0(canc_dir, x)))))

paad_centralities = merge(tcga_metadata, paad_centralities, 
                          by.x="entity_submitter_id", by.y="tcga_samp_id", all.y=T)

paad_centralities = merge(hsa_symbol_entrez, paad_centralities, 
                          by="entrez", all.y=T)

setnames(paad_centralities, "entity_submitter_id", "tcga_samp_id")

# curious what the top PPI gain/loss is for everything by diff metrics
paad_centralities[order(-weighted_degree)][,head(.SD,1),by=tcga_samp_id][,.N,by=.(symbol,entrez)][order(-N)]

paad_centralities[order(-weighted_degree)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]

paad_centralities[order(weighted_degree)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]

paad_centralities[order(-pos_betweenness)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]


paad_centralities = merge(paad_centralities, bg_centralities, by="entrez")
paad_centralities[,':='(deg_bg=bg_degree-degree, pos_deg_bg=bg_degree-pos_degree, neg_deg_bg=bg_degree-neg_degree)]

paad_clusters = fread("/grain/rad4/splitpea/out/PAAD-spectral-clust.txt")
paad_centralities = merge(paad_clusters, paad_centralities, by.x="sam", by.y="tcga_samp_id")

paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol,entrez)][order(-mean_wdegree)][,head(.SD,5),by=cluster]
paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol,entrez)][order(mean_wdegree)][,head(.SD,5),by=cluster]

paad_centralities[,.(mean_betweenness=mean(pos_betweenness)),by=.(cluster,symbol,entrez)][order(-mean_betweenness)][,head(.SD,5),by=cluster]

paad_centralities_simp = unique(paad_centralities[sam=="TCGA-2J-AAB1-01A-11R-A41B-07",
                                                  .(entrez,weighted_degree)])
paad_centralities_enrich_list = paad_centralities_simp$weighted_degree
names(paad_centralities_enrich_list) = as.character(paad_centralities_simp$entrez)
paad_centralities_enrich_list = sort(paad_centralities_enrich_list, decreasing=T)
paad_centralities_wdeg_enrich = GSEA(paad_centralities_enrich_list,
                                     TERM2GENE = cosmic_tier1)
library(org.Hs.eg.db)

paad_centralities_wdeg_go_enrich =  gseGO(geneList     = paad_centralities_enrich_list,
                                          OrgDb        = org.Hs.eg.db,
                                          ont          = "BP",
                                          minGSSize    = 100,
                                          maxGSSize    = 500,
                                          pvalueCutoff = 0.05,
                                          verbose      = FALSE)

top5_per_clust = paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]
top5_per_clust = rbind(top5_per_clust,
                       paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster])

plot_top5 = paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
pdf("results/paad_clusters.top_nodes.pdf", width=5, height=3)
ggplot(plot_top5) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster)),stat="identity", width=0.5, position="dodge") +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12)
dev.off()


canc_dir = "output/BRCA/mean/"
brca_centralities = rbindlist(lapply(list.files(canc_dir, pattern="*centralities.txt"),
                                     function(x) data.table(tcga_samp_id = gsub(".centralities.txt", "", x),
                                                            fread(paste0(canc_dir, x)))))

brca_centralities = merge(hsa_symbol_entrez, brca_centralities, 
                          by="entrez", all.y=T)

brca_centralities = merge(tcga_metadata, brca_centralities, 
                          by.x="entity_submitter_id", by.y="tcga_samp_id", all.y=T)
setnames(brca_centralities, "entity_submitter_id", "tcga_samp_id")

brca_centralities[order(-weighted_degree)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]

brca_centralities[order(weighted_degree)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]

brca_centralities[order(-pos_betweenness)][,head(.SD,1),by=.(tcga_samp_id,primary_diagnosis)][,.N,by=.(symbol,entrez,primary_diagnosis)][order(-N)]

brca_clusters = fread("/grain/rad4/splitpea/out/BRCA-spectral-clust.txt")
brca_centralities = merge(brca_clusters, brca_centralities, by.x="sam", by.y="tcga_samp_id")

brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol,entrez)][order(-mean_wdegree)][,head(.SD,5),by=cluster]
brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol,entrez)][order(mean_wdegree)][,head(.SD,5),by=cluster]

brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster]
brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]

brca_centralities[,.(mean_betweenness=mean(pos_betweenness)),by=.(cluster,symbol,entrez)][order(-mean_betweenness)][,head(.SD,5),by=cluster]


top5_per_clust = brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]
top5_per_clust = rbind(top5_per_clust,
                       brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster])

plot_top5 = brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
pdf("results/brca_clusters.top_nodes.pdf", width=5, height=3)
ggplot(plot_top5) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster)),stat="identity", width=0.5, position="dodge") +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12)
dev.off()


# consensus network stats
consensus_stats_all = fread("results/consensus_threshold.sizes.txt")
consensus_stats_lcc = fread("results/consensus_threshold.lcc_sizes.txt")

ggplot(consensus_stats_all, aes(x=threshold, y=num_nodes)) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

ggplot(consensus_stats_all, aes(x=threshold, y=num_edges)) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

ggplot(consensus_stats_lcc, aes(x=threshold, y=num_nodes)) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

ggplot(consensus_stats_lcc, aes(x=threshold, y=log2(num_edges))) + geom_point() + geom_line() +
  facet_grid(cancer~direction)

# 80% for brca so that we have things in both clusters
pdf("results/consensus_threshold.lcc_sizes.pdf", width=6, height=4)
ggplot(consensus_stats_lcc[cancer=="BRCA"],
       aes(x=threshold, y=num_nodes, fill=direction, color=direction)) + 
  geom_point() + geom_line() + geom_vline(xintercept=0.8, lty='dashed') +
  scale_color_manual(values = c("positive" = "#146CAD", "negative" = "#E63846")) + 
  xlab("consensus network threshold") + 
  ylab("# nodes") +
  theme_bw(base_size=16)
dev.off()

