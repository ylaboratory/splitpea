'''
collapse a spliced exon file into 2 new files
one with coordinates combined by median and one with average

TODO add a threshold for how many NaN values before
throwing out an exon
'''
import os
from pathlib import Path
from statistics import mean, median

BASE_DIR = Path(__file__).resolve().parent.parent
iris_path = str(BASE_DIR) + "/examples"


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


for fn in os.listdir(iris_path):
    if fn.startswith('splicing_matrix.SE.cov10'):
        print(fn)
        out_mean = open(iris_path + "/" + fn.replace('splicing_matrix.SE.cov10.', "").replace('.txt', '') + '_combined_mean.txt', 'w')
        out_med = open(iris_path + "/" + fn.replace('splicing_matrix.SE.cov10.', "").replace('.txt', '') + '_combined_median.txt', 'w')
        out_mean.write('ensembl_gene_id\tsymbol\tchr\tstrand\texon_start\texon_end\tpsi\n')
        out_med.write('ensembl_gene_id\tsymbol\tchr\tstrand\texon_start\texon_end\tpsi\n')
        header = True
        with open(iris_path + "/" + fn, 'r') as f:
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
                m = "NaN"
                if len(nums) > 1:
                    a = mean(nums)
                    m = median(nums)
                out_mean.write(s + str(a) + "\n")
                out_med.write(s + str(m) + "\n")
    out_mean.close()
    out_med.close()


