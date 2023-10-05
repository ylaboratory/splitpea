'''
given a command line argument containing a directory of 
spliced exon files, collapse each spliced exon file into
one with coordinates summarized by mean

input files in the input directory must be structured as follows:
- tab delimited
- columns order: GeneID, geneSymbol, chr, strand, exonStart, exonEnd, 
- upstreamEE, downstreamES, [remaining columns one per sample]
(example spliced exon files are in the example directory)
'''
import os, sys
from pathlib import Path
from statistics import mean

BASE_DIR = Path(__file__).resolve().parent.parent

if len(sys.argv) > 1:
    in_dir = sys.argv[1]
    if not os.path.exists(in_dir):
        sys.exit('Error: Input directory does not exist.')
else:
    sys.exit('Error: Please include an input directory as a command line argument. See the README for more information.')


def get_nums_list(s):
    '''
    given a list of string numbers
    with possible NaNs return
    all non-NaN values
    '''
    l = []

    for elm in s:
        if elm != 'NaN':
            l.append(float(elm))

    return l


for fn in os.listdir(in_dir):
    if fn.endswith('combined_mean.txt'): # skip generated output files 
        continue
    print("processing file: " + fn)
    out_mean = open(in_dir + "/" + fn.replace('.txt', '') + '_combined_mean.txt', 'w')
    out_mean.write('ensembl_gene_id\tsymbol\tchr\tstrand\texon_start\texon_end\tpsi\n')

    header = True
    with open(in_dir + "/" + fn, 'r') as f:
        for line in f:
            if header:
                header = False
                continue
            words = line.strip().split('\t')
            strand = "1"
            if words[3] == '-':
                strand = "-1"
            s = words[0] + "\t" + words[1] + "\t" + words[2].replace("chr", "") + "\t" + strand + "\t" + words[4] + '\t' + words[5] + "\t"
            nums = get_nums_list(words[8:])
            a = "NaN"
            if len(nums) > 1:
                a = mean(nums)
            out_mean.write(s + str(a) + "\n")
    out_mean.close()


