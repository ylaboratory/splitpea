'''
count up the number of tumor - normal pairs
per cancer
'''
import os

tcga_meta = {}
header = True
tmp_counts = {'primary': 0, 'norm': 0}

with open('examples/TCGA/TCGA_metadata.csv', 'r') as f:
    for line in f:
        if header:
            header = False
            continue
        words = line.replace('"', '').strip().split(',')
        barcode = words[4].replace('-', '_')
        if barcode not in tcga_meta:
            tcga_meta[barcode] = []
        tcga_meta[barcode].append((words[5], words[9]))
        if 'Primary Solid Tumor' in words[5]:
            tmp_counts['primary'] += 1
        else:
            tmp_counts['norm'] += 1

#print(tcga_meta.keys())
print('total count in metadata: ' + str(tmp_counts))

for fn in os.listdir('examples/TCGA/'):
    if fn == 'TCGA_metadata.csv':
        continue
    print(fn)
    samples = []
    mapped = 0
    norm_mapped = 0
    unmapped = set()
    with open('examples/TCGA/' + fn, 'r') as f:
        for line in f:
            words = line.strip().split("\t")
            samples = words[10:]
            break
    for s in samples:
        if s not in tcga_meta.keys():
            if '_Norm' in s:
                if s.replace('_Norm', '') in tcga_meta.keys():
                    norm_mapped += 1
                    continue
            unmapped.add(s)
        else:
            mapped += 1
    print('total mapped: ' + str(mapped) + ' norm mapped: ' + str(norm_mapped) + ' vs unmapped: ' + str(len(unmapped)))


            


