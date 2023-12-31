# ************************************************
# calculates delta psi values and
# related p-values for GTEx TCGA pairs
# ************************************************
if (!require(data.table)) install.packages('data.table')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(optparse)) install.packages('optparse')

library(data.table)
library(ggplot2)
library(optparse)


# get the file paths
option_list = list(
  make_option(c("-s", "--sumbackground"), type="character", default=NULL, 
              help="mean summarized background file", metavar="character"),
  make_option(c("-b", "--background"), type="character", default=NULL, 
              help="background spliced exon file", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=NULL, 
              help="cancer spliced exon file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output file directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

outdir = paste0(opt$outdir, "/")
gtex = opt$sumbackground
gtex.fn = opt$background
tcga.fn = opt$target

gtex = fread(gtex, sep='\t')
gtex.all = fread(gtex.fn, sep = '\t')
tcga.all = fread(tcga.fn, sep="\t")


# combine any remaining duplicate coordinates by aggregating via mean
gtex = gtex[, lapply(.SD, mean, na.rm=TRUE), by=.(ensembl_gene_id, symbol, chr, strand, exon_start, exon_end) ]
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

if(!dir.exists(paste0(outdir, '/plots/')))
{
  dir.create(paste0(outdir, '/plots/'))
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
  
  fwrite(joint,paste0(outdir, sam,'-psi.txt'), sep='\t', quote=F)
  
  ggplot(joint, aes(x=delta.psi, y=-log10(pval + 0.001))) + geom_point() + theme_bw() +
    xlim(-1, 1) + ggtitle(sam)
  
  ggsave(paste0(outdir, '/plots/', sam,'-volcano.pdf'), width = 6, height = 6)
}