# ************************************************
# cluster the delta psi values directly
# to search for relationships between samples
# ************************************************
if (!require(data.table)) install.packages('data.table')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(gplots)) install.packages('gplots')
if (!require(matrixStats)) install.packages('matrixStats')
if (!require(viridis)) install.packages('viridis')
if (!require(rstudioapi)) install.packages('rstudioapi')

library(data.table)
library(ggplot2)
library(gplots)
library(matrixStats)
library(viridis)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../utils/heatmap3.R")


# PAAD heatmap
# read in the delta.psi values across files
# and store into a matrix
fns <- list.files("../../IRIS/PAAD-psi", pattern="*.txt", full.names=TRUE)
sam.names <- gsub("../../IRIS/PAAD-psi", '', fns)
sam.names <- gsub('-psi.txt', '', sam.names)
l <- lapply(fns, fread, sep="\t")
names(l) <- sam.names
dt <- rbindlist(l, idcol=T)
names(dt) <- c('sam', 'ensembl', 'symbol', 'chr', 'strand', 'exon.start', 
               'exon.end', 'psi.gtex', 'psi.tcga', 'delta.psi', 'pval')
dt[, pos := paste0(symbol,'-', exon.start, '-', exon.end)]
mat <- dcast(dt, pos~sam, fun.aggregate=mean, value.var = 'delta.psi')

# add NAN count
mat[, num_obs := Reduce(`+`, lapply(.SD,function(x) is.nan(x)))]

# remove rows with NaNs
mat.clean = mat[num_obs == 0]
mat.clean[,num_obs := NULL]
mat.clean$rvar = rowVars(as.matrix(mat.clean[,-1]))
mat.clean = mat.clean[rvar >= 0.01]
mat.clean[, rvar := NULL]
m <- as.matrix(mat.clean[,-1])

# base heatmap
heatmap.2(m, trace = 'none', density.info = 'none', scale='none',
          labRow = F, dendrogram='column', col=viridis, keysize=1)

# add meta data
tcga.meta = fread('../../examples/TCGA_metadata.csv', sep=',')
tcga.meta = tcga.meta[entity_submitter_id %in% names(mat.clean)]
tcga.meta = tcga.meta[, c('entity_submitter_id', 'sample_types', 'center', 'ajcc_pathologic_stage',
                          'tissue_or_organ_of_origin', 'days_to_last_follow_up', 'primary_diagnosis',
                          'year_of_diagnosis', 'vital_status', 'race', 'gender', 'ethnicity', 'age_at_index',
                          'site_of_resection_or_biopsy')]


c_annot = data.table(sam=names(mat.clean[,-1]))
c_annot = merge(c_annot, tcga.meta, by.x='sam', by.y = 'entity_submitter_id', all.x=T)

c_annot[, gender_c := ifelse(gender == 'male', "#1d3557", "#e63946")]
c_annot[, age_c := ifelse(age_at_index < 20, '#eae2b7', "#fcbf49")]
c_annot[, age_c := ifelse(age_at_index >= 40, '#f77f00', age_c)]
c_annot[, age_c := ifelse(age_at_index >= 60, '#d62828', age_c)]
c_annot[, age_c := ifelse(age_at_index >= 80, '#003049', age_c)]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage IIA' | ajcc_pathologic_stage == 'Stage IIB', '#7FCDBB', "#ECF7B0")]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage III', '#2091C0', stage_c)]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage IV', '#263494', stage_c)]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Adenocarcinoma with mixed subtypes', 'red', 'purple')]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Mucinous adenocarcinoma', 'pink', diag_c)]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Carcinoma, undifferentiated, NOS', 'green', diag_c)]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Infiltrating duct carcinoma, NOS', 'darkgreen', diag_c)]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Neuroendocrine carcinoma, NOS', 'blue', diag_c)]
c_annot[, site_c := ifelse(site_of_resection_or_biopsy == 'Head of pancreas', 'red', 'navy')]
c_annot[, site_c := ifelse(site_of_resection_or_biopsy == 'Tail of pancreas', 'green', site_c)]
c_annot[, site_c := ifelse(site_of_resection_or_biopsy == 'Overlapping lesion of pancreas', 'black', site_c)]
c_annot[, site_c := ifelse(site_of_resection_or_biopsy == 'Body of pancreas', 'blue', site_c)]


ann = c_annot[,c('gender_c', 'age_c', 'stage_c', 'diag_c', 'site_c')]
colnames(ann) <- c("sex", "age", 'stage', 'diagnosis', 'site')

# multicolumn version
heatmap.3(m, trace = 'none', density.info = 'none', scale='none', labCol = F,
          labRow = F, dendrogram='column', col=viridis,
          ColSideColors=as.matrix(ann))



# ----------------------------
#    MAKE HEATMAP FOR BRCA
# ----------------------------
# read in the delta.psi values across files
# and store into a matrix
fns <- list.files("../../IRIS/BRCA-psi", pattern="*.txt", full.names=TRUE)
sam.names <- gsub("../../IRIS/BRCA-psi", '', fns)
sam.names <- gsub('-psi.txt', '', sam.names)
l <- lapply(fns, fread, sep="\t")
names(l) <- sam.names
dt <- rbindlist(l, idcol=T)
names(dt) <- c('sam', 'ensembl', 'symbol', 'chr', 'strand', 'exon.start', 
               'exon.end', 'psi.gtex', 'psi.tcga', 'delta.psi', 'pval')
dt[, pos := paste0(symbol,'-', exon.start, '-', exon.end)]
mat <- dcast(dt, pos~sam, fun.aggregate=mean, value.var = 'delta.psi')

# add NAN count
mat[, num_obs := Reduce(`+`, lapply(.SD,function(x) is.nan(x)))]

# remove rows with NaNs
mat.clean = mat[num_obs == 0]
mat.clean[,num_obs := NULL]
mat.clean$rvar = rowVars(as.matrix(mat.clean[,-1]))
mat.clean = mat.clean[rvar >= 0.01]
mat.clean[, rvar := NULL]
m <- as.matrix(mat.clean[,-1])

# base heatmap
heatmap.2(m, trace = 'none', density.info = 'none', scale='none',
          labRow = F, dendrogram='column', col=viridis, keysize=1)

# add meta data
brca.meta = fread('../../examples/BRCA_pam50.txt', sep="\t")
brca.meta[, barcode := gsub('.{3}$', '', `Sample ID`)]
brca.meta[, c("PX", "PY", 'PZ', 'type') := tstrsplit(`Sample ID`, "-", fixed=TRUE)]
brca.meta[, PX := NULL]
brca.meta[, PY := NULL]
brca.meta[, PZ := NULL]
brca.meta = brca.meta[type != '11'] # remove the normal annotations


tcga.meta = fread('../../examples/TCGA_metadata.csv', sep=',')
tcga.meta = tcga.meta[entity_submitter_id %in% names(mat.clean)]
tcga.meta = tcga.meta[, c('entity_submitter_id', 'TCGABarcode', 'sample_types', 'center', 'ajcc_pathologic_stage',
                          'tissue_or_organ_of_origin', 'days_to_last_follow_up', 'primary_diagnosis',
                          'year_of_diagnosis', 'vital_status', 'race', 'gender', 'ethnicity', 'age_at_index',
                          'site_of_resection_or_biopsy')]

tcga.meta = merge(brca.meta, tcga.meta, by.y='TCGABarcode', by.x='barcode', all.y=T)

c_annot = data.table(sam=names(mat.clean[,-1]))
c_annot = merge(c_annot, tcga.meta, by.x='sam', by.y = 'entity_submitter_id', all.x=T)

c_annot[, age_c := ifelse(age_at_index < 20, '#eae2b7', "#fcbf49")]
c_annot[, age_c := ifelse(age_at_index >= 40, '#f77f00', age_c)]
c_annot[, age_c := ifelse(age_at_index >= 60, '#d62828', age_c)]
c_annot[, age_c := ifelse(age_at_index >= 80, '#003049', age_c)]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage IIA' | ajcc_pathologic_stage == 'Stage IIB' | ajcc_pathologic_stage == 'Stage II', '#7FCDBB', "#ECF7B0")]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage III' | ajcc_pathologic_stage == 'Stage IIIA' | ajcc_pathologic_stage == 'Stage IIIB', '#2091C0', stage_c)]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage IV', '#263494', stage_c)]
c_annot[, stage_c := ifelse(ajcc_pathologic_stage == 'Stage X', 'black', stage_c)]

c_annot[, diag_c := ifelse(primary_diagnosis == 'Lobular carcinoma, NOS', 'red', 'purple')]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Infiltrating duct carcinoma, NOS', 'blue', diag_c)]
c_annot[, diag_c := ifelse(primary_diagnosis == 'Infiltrating duct and lobular carcinoma', 'green', diag_c)]

c_annot[, pam := ifelse(PAM50 == 'Basal', 'red', 'navy')]
c_annot[, pam := ifelse(PAM50 == 'Her2', 'green', pam)]
c_annot[, pam := ifelse(PAM50 == 'LumB', 'blue', pam)]
c_annot[, pam := ifelse(PAM50 == 'Normal', 'grey', pam)]



ann = c_annot[,c('age_c', 'stage_c', 'diag_c', 'pam')]
colnames(ann) <- c("age", 'stage', 'diagnosis', 'pam50')

# multicolumn version
heatmap.3(m, trace = 'none', density.info = 'none', scale='none', labCol = F,
          labRow = F, dendrogram='column', col=viridis,
          ColSideColors=as.matrix(ann))

