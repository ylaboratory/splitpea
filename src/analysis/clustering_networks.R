# code for the analysis depicted in figure 6
library('data.table')
library('ggplot2')
library('dbscan')

combo_embed = fread('examples/BRCA-PAAD-FEATHER-matrix.txt', header=T)

tcga_meta = fread('examples/TCGA/TCGA_metadata.csv')
brca_meta = fread('examples/TCGA/BRCA_additional_annotations.txt')

tcga_meta = tcga_meta[, c('entity_submitter_id', 'sample_types', 'center', 'ajcc_pathologic_stage',
                                  'tissue_or_organ_of_origin', 'days_to_last_follow_up', 'days_to_death', 'primary_diagnosis',
                                  'year_of_diagnosis', 'vital_status', 'race', 'gender', 'ethnicity', 'age_at_index',
                                  'site_of_resection_or_biopsy', 'disease')]


# get centralities
hsa_gene_mappings = fread("examples/hsa_mapping_all.txt")
hsa_symbol_entrez = unique(hsa_gene_mappings[,.(symbol, entrez)])

bg_centralities = fread("reference/human_ppi_ddi_bg.centralities.txt")
setnames(bg_centralities, c("entrez", "bg_degree", "bg_betweenness"))

canc_dir = "/grain/vy3/splitpea/output/PAAD/mean/"
paad_centralities = rbindlist(lapply(list.files(canc_dir, pattern="*centralities.txt"),
                                     function(x) data.table(tcga_samp_id = gsub(".centralities.txt", "", x),
                                                            fread(paste0(canc_dir, x)))))

paad_centralities = merge(hsa_symbol_entrez, paad_centralities, 
                          by="entrez", all.y=T)

paad_centralities[order(-weighted_degree)][,head(.SD,1),by=tcga_samp_id][,.N,by=.(symbol,entrez)][order(-N)]
paad_centralities = merge(paad_centralities, bg_centralities, by="entrez")
paad_centralities[,':='(deg_bg=bg_degree-degree, pos_deg_bg=bg_degree-pos_degree, neg_deg_bg=bg_degree-neg_degree)]



canc_dir = "/grain/vy3/splitpea/output/BRCA/mean/"
brca_centralities = rbindlist(lapply(list.files(canc_dir, pattern="*centralities.txt"),
                                     function(x) data.table(tcga_samp_id = gsub(".centralities.txt", "", x),
                                                            fread(paste0(canc_dir, x)))))

brca_centralities = merge(hsa_symbol_entrez, brca_centralities, 
                          by="entrez", all.y=T)

brca_centralities[order(-weighted_degree)][,head(.SD,1),by=tcga_samp_id][,.N,by=.(symbol,entrez)][order(-N)]
brca_centralities = merge(brca_centralities, bg_centralities, by="entrez")
brca_centralities[,':='(deg_bg=bg_degree-degree, pos_deg_bg=bg_degree-pos_degree, neg_deg_bg=bg_degree-neg_degree)]



sams = PAAD_embeddings[,sam]

tmp = PAAD_embeddings
tmp = tmp[, sam := NULL]
pca = prcomp(tmp, scale = F)

paad.dim.red = data.table(sam=sams, pc1=pca$x[,1],pc2=pca$x[,2])
paad.dim.red = merge(paad.dim.red, tcga_meta, by.x='sam', by.y='entity_submitter_id')

ggplot(paad.dim.red, aes(x=pc1, y=pc2, color=sample_types)) + geom_point() + theme_bw()
ggplot(paad.dim.red, aes(x=pc1, y=pc2, color=primary_diagnosis)) + geom_point() + theme_bw()
ggplot(paad.dim.red, aes(x=pc1, y=pc2, color=age_at_index)) + geom_point() + theme_bw()


db = hdbscan(paad.dim.red[,2:3], minPts=10)
db

paad.dim.red[, cluster := db$cluster]
ggplot(paad.dim.red, aes(x=pc1, y=pc2, color=as.factor(cluster))) + geom_point() + theme_bw()

paad_centralities = merge(paad_centralities, paad.dim.red, by.x="tcga_samp_id", by.y='sam')

top5_per_clust = paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]
top5_per_clust = rbind(top5_per_clust,
                       paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster])

plot_top5 = paad_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
ggplot(plot_top5) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster)),stat="identity", width=0.5, position="dodge") +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12)



# also analyze brca data
sams = BRCA_embeddings[,sam]

tmp = BRCA_embeddings
tmp = tmp[, sam := NULL]
pca = prcomp(tmp, scale = F)

brca.dim.red = data.table(sam=sams, pc1=pca$x[,1],pc2=pca$x[,2])
brca.dim.red = merge(brca.dim.red, tcga_meta, by.x='sam', by.y='entity_submitter_id')

ggplot(brca.dim.red, aes(x=pc1, y=pc2, color=sample_types)) + geom_point() + theme_bw()
ggplot(brca.dim.red, aes(x=pc1, y=pc2, color=primary_diagnosis)) + geom_point() + theme_bw()
ggplot(brca.dim.red, aes(x=pc1, y=pc2, color=ajcc_pathologic_stage)) + geom_point() + theme_bw()

db = hdbscan(brca.dim.red[,2:3], minPts=10)
db

brca.dim.red[, cluster := db$cluster]
ggplot(brca.dim.red, aes(x=pc1, y=pc2, color=as.factor(cluster))) + geom_point() + theme_bw()

brca_centralities = merge(brca_centralities, brca.dim.red, by.x="tcga_samp_id", by.y='sam')

top5_per_clust = brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]
top5_per_clust = rbind(top5_per_clust,
                       all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster])

plot_top5 = brca_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
ggplot(plot_top5) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster)),stat="identity", width=0.5, position="dodge") +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12)


# all together
sams = combo_embed[,sam]

tmp = combo_embed
tmp = tmp[, sam := NULL]
pca = prcomp(tmp, scale = F)

combo.dim.red = data.table(sam=sams, pc1=pca$x[,1],pc2=pca$x[,2])
combo.dim.red = merge(combo.dim.red, tcga_meta, by.x='sam', by.y='entity_submitter_id')

ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=disease)) + geom_point() + theme_bw()
ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=primary_diagnosis)) + geom_point() + theme_bw()
ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=age_at_index)) + geom_point() + theme_bw()


db = hdbscan(combo.dim.red[,2:3], minPts=10)
db

combo.dim.red[, cluster := db$cluster]
combo.dim.red[, cluster2 := cluster]
combo.dim.red[, cluster2 := ifelse(cluster2 == 0, paste0(cluster2, disease), cluster2)]

ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=as.factor(cluster))) + geom_point() + theme_classic() +
  scale_color_manual(values=c('#c0c0c0', '#ff7d00', '#78290f', '#15616d')) + labs(color='cluster')
ggsave('/grain/rad4/splitpea/out/plots/pca-all-nets-dbclust.pdf', units='in', width=6, height=5)

ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=as.factor(disease))) + geom_point() + theme_classic() +
  scale_color_manual(values=c('#E69F00', '#0072B2')) + labs(color='cancer') 
ggsave('/grain/rad4/splitpea/out/plots/pca-all-nets.pdf', units='in', width=6, height=5)

ggplot(combo.dim.red, aes(x=pc1, y=pc2, color=as.factor(cluster2))) + geom_point() + theme_minimal()



all_centralities = rbind(brca_centralities, paad_centralities)
all_centralities = merge(all_centralities, combo.dim.red, by.x="tcga_samp_id", by.y='sam')

top5_per_clust = all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster]
top5_per_clust = rbind(top5_per_clust,
                       all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster])

plot_top5 = all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
ggplot(plot_top5) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster)),stat="identity", width=0.5, position="dodge") +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12)


# combo 2
all_centralities = rbind(brca_centralities, paad_centralities)
all_centralities = merge(all_centralities, combo.dim.red, by.x="tcga_samp_id", by.y='sam')

top5_per_clust = all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster2,symbol)][order(mean_wdegree)][,head(.SD,5),by=cluster2]
top5_per_clust = rbind(top5_per_clust,
                       all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster2,symbol)][order(-mean_wdegree)][,head(.SD,5),by=cluster2])

plot_top5 = all_centralities[,.(mean_wdegree=mean(weighted_degree)),by=.(cluster2,symbol)][symbol %in% top5_per_clust$symbol]
plot_top5$symbol = factor(plot_top5$symbol, levels = plot_top5[,max(mean_wdegree),by=symbol][order(V1)]$symbol)
gain_genes = c('KRAS', 'CCL5', 'CDK4', 'MELK', 'RAB5A', 'NEDD8', 'AKT2', 'RPS3', 'IKBKB', 'PTPN18', 'P4HB', 'RHOF', 'KIF2A')
plot_top5[, direction := ifelse(symbol %in% gain_genes, "gain", 'lost')]

ggplot(plot_top5[direction == 'gain']) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster2)),stat="identity", position=position_dodge2(reverse = TRUE)) +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12) + theme_classic() + labs(fill='cluster') + 
  scale_fill_manual(values=c('#c0c0c0', '#696969', '#ff7d00', '#78290f', '#15616d'))
ggsave('/grain/rad4/splitpea/out/plots/net-clust-gains.pdf', units='in', width=6, height=5)


ggplot(plot_top5[direction == 'lost']) +
  geom_bar(aes(x=symbol, y=mean_wdegree,fill=as.factor(cluster2)),stat="identity", width=0.5, position=position_dodge2(reverse = TRUE)) +
  xlab("") + ylab("mean interactions gained / lost") +
  coord_flip() + theme_classic(base_size=12) + theme_classic() + labs(fill='cluster') + 
  scale_fill_manual(values=c('#c0c0c0', '#696969', '#ff7d00', '#78290f', '#15616d'))
ggsave('/grain/rad4/splitpea/out/plots/net-clust-losses.pdf', units='in', width=6, height=5)

