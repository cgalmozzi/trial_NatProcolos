############################
#                          #
#     Total_Enrichment     #
#                          #
############################

'''
Created on Tue Sep 25 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script calculates the total enrichment (TE) value for each gene comparing 
two biological replicates of total and selective samples. 

The output is a text file with the following columns: 
systematic gene name 
(trivial) gene name
for each of the four samples: 
    (normalized) gene expression value [RPKM]
    'included' or 'excluded'
ratio - replicate 1
ratio - replicate 2
average of both
log2 transform of average 
number of samples with 'included' 

Note: If the ratio, the average or the log2 transform cannot be calculated 
due to e.g. zero footprints for a gene, the respective value is given as 'n.a.'.


'''

import math
import os
import pickle


def calcTE(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name): 

    # read input data    
    data = {}
    spec = {}
    proc = {}
    
    for file in [file_total1, file_total2, file_selec1, file_selec2]: 
        with open(input_path + file + '_gene expression.txt', 'r') as f:
            for line in f:
                fields = line.split('\t')
                gene = fields[0]
                ge = float(fields[3])
                sp = fields[4].strip()
                data.setdefault(gene, []).append(ge)
                spec.setdefault(gene, []).append(sp)

    # reference file for gene names
    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = path_current + '/references_yeast/'
    dictGenes = pickle.load(open(path_ref + 'yeast_genes.pkl', 'rb'))

    # process data / calculate ratio and TE
    for (gene, values) in data.items():
        try:
            ratio1 = values[2] / values[0]
        except ZeroDivisionError: 
            ratio1 = 'n.a.'
        try:
            ratio2 = values[3] / values[1]
        except ZeroDivisionError: 
            ratio2 = 'n.a.'
        try:
            TE = (ratio1 + ratio2) / 2
        except (TypeError, ValueError):
            TE = 'n.a.'
        try:
            log2 = math.log(TE)
        except (TypeError, ValueError):
            log2 = 'n.a.'
        specNO = str(spec[gene].count('included')) + ' of 4'
        proc[gene] = [ratio1, ratio2, TE, log2, specNO]

    # write output file
    listGenes = list(data.keys())
    listGenes.sort()
    header = 'systematic gene name\tgene name\ttotal 1\t\ttotal2\t\tselective 1\
    \t\tselective 2\t\tratio 1\tratio 2\taverage\taverage [log2]\tno. of included values\n'

    with open(input_path + output_name + '_TE.txt', 'w') as f: 
        f.write(header)
        for gene in listGenes:
            f.write(gene + '\t' + dictGenes[gene][0] + '\t')
            for n in range(4):
                f.write(str(data[gene][n]) + '\t' + spec[gene][n] + '\t')
            for value in proc[gene]:
                f.write(str(value) + '\t')
            f.write('\n')

    return


if __name__ == '__main__':

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    file_total1 = 'sample name'                 # sample name total translatome 1 (no file extension)
    file_total2 = 'sample name'                 # sample name total translatome 2 (no file extension)
    file_selec1 = 'sample name'                 # sample name selective translatome 1 (no file extension)
    file_selec2 = 'sample name'                 # sample name selective translatome 2 (no file extension)
    output_name = 'experiment name'             # experiment name (no file extension)

    calcTE(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name)








