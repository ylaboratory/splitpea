# ************************************************
# this script calculates delta psi values and
# related p-values for GTEx TCGA pairs
# ************************************************
library(data.table)
library(ggplot2)

# ------------
# mean
# ------------
can_type = 'PAAD'
outdir = paste0("examples/", can_type, "/")
gtex.fn = "examples/splicing_matrix.SE.cov10.GTEx_Pancreas.txt"
tcga.fn = paste0('examples/splicing_matrix.SE.cov10.TCGA_',can_type,'_T.txt')
summarized.prefix = 'examples/GTEx_Pancreas_combined'

gtex = fread(paste0(summarized.prefix,'_mean.txt'), sep='\t')
gtex.all = fread(gtex.fn, sep = '\t')
tcga.all = fread(tcga.fn, sep="\t")
#tcga.paad = fread('examples/TCGA/IRIS/TCGA_PAAD_T_combined_mean.txt', sep="\t")

# combine any remaining duplicate coordinates by aggregating via mean
gtex = gtex[, lapply(.SD, mean, na.rm=TRUE), by=.(ensembl_gene_id, symbol, chr, strand, exon_start, exon_end) ]
#tcga.paad = tcga.paad[, lapply(.SD, mean, na.rm=TRUE), by=.(`ensembl_gene_id`, symbol, chr, strand, exon_start, exon_end) ]
tcga.all[, upstreamEE := NULL]
tcga.all[, downstreamES := NULL]
tcga.all = tcga.all[, lapply(.SD, mean, na.rm=TRUE), by=.(AC, GeneName, chr, strand, exonStart, exonEnd) ]
tcga.all = tcga.all[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
tcga.all[, chr := gsub('chr', "", chr)]
tcga.all[, strand := as.character(strand)]
tcga.all[, strand := ifelse(strand == "-", -1, 1)]

gtex.all[, upstreamEE := NULL]
gtex.all[, downstreamES := NULL]
gtex.all = gtex.all[, lapply(.SD, mean, na.rm=TRUE), by=.(GeneID, geneSymbol, chr, strand, exonStart_0base, exonEnd) ]
gtex.all = gtex.all[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
gtex.all[, chr := gsub('chr', "", chr)]
gtex.all[, strand := as.character(strand)]
gtex.all[, strand := ifelse(strand == "-", -1, 1)]

gtex.sams = names(gtex.all)
gtex.sams = gtex.sams[7:length(gtex.sams)]
names(gtex.all) = append(c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'), gtex.sams)

# make the backgrounds for p-value calculation
getECDF <- function(vals) {
  min_valid = 10 # or 10% of samples
  v = as.numeric(vals)
  if ((length(v) *.1) > min_valid)
  {
    min_valid = length(v) * 0.1
  }
  
  # discard cases with less than min valid cases
  if (sum(!is.na(v)) < min_valid)
  {
    return(function(x){return(rep(NA,length(x)))})
  }
  
  d <- abs(outer(v,v,"-"))
  subtracts <- d[lower.tri(d, diag = FALSE)]
  
  if (all(subtracts[!is.na(subtracts)]==0)) # if all values the same
  {
    return(function(x){return(rep(NA,length(x)))})
    # return(function(x){
    #   p <- rep(NA,length(x))
    #   p[x==0] <- 0
    #   p[x>0] <- 1
    #   return(p)})
  } 
  else 
  {
    return(ecdf(subtracts))
  }
}
# make merged set so order of exons is preserved
# then precompute the empirical distributions
gtex.tmp.ecdf = apply(gtex.all[,gtex.sams, with=F], 1, getECDF)
gtex.ecdf = gtex.all
gtex.ecdf[, ecdf := gtex.tmp.ecdf]
gtex.ecdf = gtex.ecdf[,c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end', 'ecdf'), with=F]

# implement the directory structure
if (!dir.exists(outdir))
{
  dir.create(outdir)
}
if (!dir.exists(paste0(outdir, "mean/")))
{
  dir.create(paste0(outdir, "mean/"))
}
if (!dir.exists(paste0(outdir, "median/")))
{
  dir.create(paste0(outdir, "median/"))
}
if(!dir.exists(paste0('output/plots/mean/', can_type)))
{
  dir.create(paste0('output/plots/mean/', can_type))
}
if(!dir.exists(paste0('output/plots/median/', can_type)))
{
  dir.create(paste0('output/plots/median/', can_type))
}


# for each tcga sample calculate delta PSI
tcga.sams = names(tcga.all)
tcga.sams = tcga.sams[7:length(tcga.sams)] # remove all but the tcga sample identifiers

for (sam in tcga.sams)
{
  print(sam)
  # calculate delta PSI by subtraction
  tmp = tcga.all[,append(c(names(tcga.all)[1:6]), sam),with=F]
  names(tmp) = c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end', sam)
  joint = merge(gtex, tmp, by=c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'), all=T)
  names(joint) <- c('ensembl.id', 'symbol', 'chr', 'strand', 'exon.start', 'exon.end', 'psi.gtex', 'psi.tcga')
  joint = joint[!is.na(psi.gtex)]
  joint = joint[!is.na(psi.tcga)]
  joint[, delta.psi := psi.tcga - psi.gtex]
  
  # calculate the pvalue using empirical distribution
  joint.tmp = merge(joint, gtex.ecdf, 
                    by.x = c('ensembl.id', 'symbol', 'chr', 'strand', 'exon.start', 'exon.end'), 
                    by.y=c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'))
  
  
  joint[, pval := mapply(function(f,x){f(x)}, joint.tmp$ecdf, as.list(abs(joint.tmp$delta.psi)))]
  joint[, pval := 1 - pval]
  
  fwrite(joint,paste0(outdir,"mean/", sam,'-psi.txt'), sep='\t', quote=F)
  
  ggplot(joint, aes(x=delta.psi, y=-log10(pval + 0.001))) + geom_point() + theme_bw() +
    xlim(-1, 1) + ggtitle(sam)
  
  ggsave(paste0('out/plots/mean/', can_type, "/", sam,'-volcano.pdf'), width = 6, height = 6)
}

# ------------
# median
# ------------
gtex.all = fread(gtex.fn, sep = '\t')
tcga.all = fread(tcga.fn, sep="\t")
gtex = fread(paste0(summarized.prefix,'_median.txt'), sep='\t')
gtex = gtex[, lapply(.SD, median, na.rm=TRUE), by=.(ensembl_gene_id, symbol, chr, strand, exon_start, exon_end) ]


tcga.all[, upstreamEE := NULL]
tcga.all[, downstreamES := NULL]
tcga.all = tcga.all[, lapply(.SD, median, na.rm=TRUE), by=.(AC, GeneName, chr, strand, exonStart, exonEnd) ]
tcga.all = tcga.all[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
tcga.all[, chr := gsub('chr', "", chr)]
tcga.all[, strand := as.character(strand)]
tcga.all[, strand := ifelse(strand == "-", -1, 1)]

gtex.all[, upstreamEE := NULL]
gtex.all[, downstreamES := NULL]
gtex.all = gtex.all[, lapply(.SD, median, na.rm=TRUE), by=.(GeneID, geneSymbol, chr, strand, exonStart_0base, exonEnd) ]
gtex.all = gtex.all[, lapply(.SD, function(x) replace(x, is.nan(x), NA))]
gtex.all[, chr := gsub('chr', "", chr)]
gtex.all[, strand := as.character(strand)]
gtex.all[, strand := ifelse(strand == "-", -1, 1)]

gtex.sams = names(gtex.all)
gtex.sams = gtex.sams[7:length(gtex.sams)]
names(gtex.all) = append(c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'), gtex.sams)

gtex.tmp.ecdf = apply(gtex.all[,gtex.sams, with=F], 1, getECDF)
gtex.ecdf = gtex.all
gtex.ecdf[, ecdf := gtex.tmp.ecdf]
gtex.ecdf = gtex.ecdf[,c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end', 'ecdf'), with=F]


for (sam in tcga.sams)
{
  print(sam)
  # calculate delta PSI by subtraction
  tmp = tcga.all[,append(c(names(tcga.all)[1:6]), sam),with=F]
  names(tmp) = c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end', sam)
  joint = merge(gtex, tmp, by=c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'), all=T)
  names(joint) <- c('ensembl.id', 'symbol', 'chr', 'strand', 'exon.start', 'exon.end', 'psi.gtex', 'psi.tcga')
  joint = joint[!is.na(psi.gtex)]
  joint = joint[!is.na(psi.tcga)]
  joint[, delta.psi := psi.tcga - psi.gtex]
  
  # calculate the pvalue using empirical distribution
  joint.tmp = merge(joint, gtex.ecdf, 
                    by.x = c('ensembl.id', 'symbol', 'chr', 'strand', 'exon.start', 'exon.end'), 
                    by.y=c('ensembl_gene_id', 'symbol', 'chr', 'strand', 'exon_start', 'exon_end'))
  
  
  joint[, pval := mapply(function(f,x){f(x)}, joint.tmp$ecdf, as.list(abs(joint.tmp$delta.psi)))]
  joint[, pval := 1 - pval]
  
  fwrite(joint,paste0(outdir,"median/", sam,'-psi.txt'), sep='\t', quote=F)
  
  ggplot(joint, aes(x=delta.psi, y=-log10(pval + 0.001))) + geom_point() + theme_bw() +
    xlim(-1, 1) + ggtitle(sam)
  
  ggsave(paste0('output/plots/median/', can_type, '/', sam,'-volcano.pdf'), width = 6, height = 6)
}




# quick data checks
#hist(joint$delta.psi, breaks=500)
#length(unique(joint$symbol)) # how many genes covered
#length(unique(joint[abs(delta.psi) > 0.05]$symbol)) # how many if we make arbitrary 0.05 threshold for psi




