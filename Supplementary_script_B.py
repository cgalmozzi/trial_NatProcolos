##############################
#                            #
#      3_nt_Periodicity      #
#                            #
##############################

'''
Created on Mon Sep 24 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script calculates the average read density along all ORFs with nucleotide 
resolution. From -9nt to +30nt around the translation start side, data are
plotted as bar plot with every third nucleotide colored in green. This should 
help to visualize and analyse the expected 3nt periodicity of every (selective)
ribosome profiling sample. 

The graph is saved as png file and as pdf file. 


'''

import os
import re
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def calcPeriodicity(input_path, sample_name, threshold):

    # upload input data
    dictReads = pickle.load(open(input_path + sample_name + '_reads.pkl', 'rb'))

    # reference files
    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = path_current + '/references_yeast/'
    dictGenes = pickle.load(open(path_ref + 'yeast_genes.pkl', 'rb'))
    dictIntrons = pickle.load(open(path_ref + 'yeast_introns.pkl', 'rb'))

    # process data
    dictMeta = dict([i, []] for i in range(-50, 1501))

    for gene in dictGenes.keys():
        chrom = str(dictGenes[gene][1]).zfill(2)
        strand = dictGenes[gene][2]
        pos_list = dictIntrons[gene].copy()
        maxi = int(max(re.split('\W+', dictGenes[gene][3].strip(','))))
        mini = int(min(re.split('\W+', dictGenes[gene][3].strip(','))))
        length = len(pos_list)
        ind_raw = 0 if strand == '+' else 2
        ind_norm = 6 if strand == '+' else 8

        sum_raw = 0
        sum_norm = 0
                
        for pos in pos_list:
            sum_raw += dictReads[pos][ind_raw]                                  # for 64-threshold
            sum_norm += dictReads[pos][ind_norm]                                # for normalization

        norm_val = sum_norm / length
            
        if sum_raw < threshold:                                                      # exclude genes with < 64 raw reads
            continue 

        if strand == '+':                                                       # add 50nt before each ORF
            start_pos = mini
            for i in range(-1, -51, -1): 
                pos_list.insert(0, str(chrom) + '\t' + str(start_pos + i).zfill(7))
        elif strand == '-': 
            start_pos = maxi
            for i in range(1, 51, 1): 
                pos_list.insert(0, str(chrom) + '\t' + str(start_pos + i).zfill(7))

        pos_counter = -50
        for pos in pos_list:
            if pos_counter > 1500:
                break
            dictMeta[pos_counter].append(dictReads[pos][ind_norm] / norm_val)
            pos_counter += 1

    # generate graph
    mpl.rcParams['axes.linewidth'] = 0.64
    mpl.rcParams['font.sans-serif'] = "Arial"
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 8

    fig = plt.figure(1, figsize = (3.2, 1.6))
    plt.subplots_adjust(0,0,1,1)
    ax = fig.add_subplot(111)    
    plt.tick_params(axis='y', which='both',
                    left='on', right = 'off',
                    labelleft='on', labelright='off',
                    direction = 'out', width = 0.64, length = 2.0)
    plt.tick_params(axis='x', which='both',
                    bottom='on', top='off',
                    labelbottom='on', 
                    direction = 'out', width = 0.64, length = 2.0)

    x_label = ('position along ORF [nt]')
    y_label = ('average reads [a.u.]')
    ax.set_xlabel(x_label, fontsize = 8)
    ax.set_ylabel(y_label, fontsize = 8)

    y = []
    for n in range(-9,31):
        meta_v = np.mean(dictMeta[n])
        y.append(meta_v)
    y = np.array(y)
    x = np.arange(-9,31)
    ax.bar(x, y, align='center', color=['#999999', '#999999', '#00aa00'], width=0.6, linewidth=0)
    ax.set_xlim(-10, 31)

    # save graph
    plt.savefig(input_path + sample_name + '_3 nt periodicity.png', bbox_inches='tight', dpi = 600)
    plt.savefig(input_path + sample_name + '_3 nt periodicity.pdf', bbox_inches='tight')
    plt.close()

    return


if __name__ == '__main__':

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    sample_name = 'sample name'                 # e.g. 'S1' (no filename extension)
    threshold = 64.0                            # minimal number of reads per gene to be included

    calcPeriodicity(input_path, sample_name, threshold)

