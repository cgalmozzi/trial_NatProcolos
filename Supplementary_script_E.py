#######################
#                     #
#    Gene_Profiles    #
#                     #
#######################

'''
Created on Tue Sep 25 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script generates an enrichment profile for each gene of the given 
organism. All genes, even those with a read number lower than the given 
threshold (default: 64) in all samples are included. Introns are automatically 
removed from the ORF (= x-axis). 

The mean between two biological replicates is plotted as log2 transform (=black
line) and the range between the two replicates is highlighted as grey area. 
Gaps represent codons or regions of codons that have no reads in either the 
total translatome sample or the selective translatome sample and can therefore 
not give an enrichment between selective and total data set or cannot be trans-
formed to log2 scale. 


'''


import os
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def tranformToCodons(nt_list):
    ''' This module transforms a given list of nt-specific values in codon-
    specific values. Three adjacent values are summed up. 
    ''' 
    
    codon_list = []
    
    for i in range(0,len(nt_list)-1,3):
        first_nt = nt_list[i]
        second_nt = nt_list[i+1]
        third_nt = nt_list[i+2]
        codon = (first_nt + second_nt + third_nt) / 1
        codon_list.append(codon)
        
    return (codon_list)


def calculateRatio(ip, total):
    ''' This module calculates the ratio for each position of two given lists. 
    If the total is 0, the ratio is set to np.nan.
    output: ratio_list (array-type)
    '''
    
    ratio_list = []
    length = len(ip)
    for ind in range(length):
        try: 
            ratio = ip[ind] / total[ind]
        except: 
            ratio = np.nan
        ratio_list.append(ratio)

    ratio_list = np.array(ratio_list)        
    return (ratio_list)


def removeZero(input_list): 
    ''' This module replaces all zeros in a given list by np.nan to allow the 
    meaningful plotting on a log-scale. np.nan values produce gaps in the
    plotted line while zeros would decrease the line to the lower plotting area.
    ''' 
    
    output_list = []
    for elem in input_list: 
        if elem == 0:
            output_list.append(np.nan)
        else: 
            output_list.append(elem)
    
    output_list = np.array(output_list)    
    return (output_list)



def geneProfiles(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name): 

    # upload input data
    total_1 = pickle.load(open(input_path + file_total1 + '_reads.pkl', 'rb'))
    total_2 = pickle.load(open(input_path + file_total2 + '_reads.pkl', 'rb'))
    selec_1 = pickle.load(open(input_path + file_selec1 + '_reads.pkl', 'rb'))
    selec_2 = pickle.load(open(input_path + file_selec2 + '_reads.pkl', 'rb'))

    # reference files
    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = path_current + '/references_yeast/'
    dictGenes = pickle.load(open(path_ref + 'yeast_genes.pkl', 'rb'))
    dictIntrons = pickle.load(open(path_ref + 'yeast_introns.pkl', 'rb'))

    # generate output folder
    os.mkdir(input_path + output_name + '_gene profiles')
    output_path = input_path + output_name + '_gene profiles/'

    for gene in dictGenes.keys():

        # process data        
        nt_total1 = []
        nt_total2 = []
        nt_selec1 = []
        nt_selec2 = []
        pos_list = dictIntrons[gene].copy()
        strand = dictGenes[gene][2]
        value_pos = 6 if strand == '+' else 8

        for pos in pos_list: 
            nt_total1.append(total_1[pos][value_pos])
            nt_total2.append(total_2[pos][value_pos])
            nt_selec1.append(selec_1[pos][value_pos])
            nt_selec2.append(selec_2[pos][value_pos])

        codon_total1 = tranformToCodons(nt_total1)
        codon_total2 = tranformToCodons(nt_total1)
        codon_selec1 = tranformToCodons(nt_selec1)
        codon_selec2 = tranformToCodons(nt_selec2)

        ratio_1 = removeZero(calculateRatio(codon_selec1, codon_total1))
        ratio_2 = removeZero(calculateRatio(codon_selec2, codon_total2))
        average = removeZero(np.nanmean(np.array([ratio_1, ratio_2]), axis = 0))

        # generate graph
        mpl.rcParams['axes.linewidth'] = 0.64
        mpl.rcParams['font.sans-serif'] = "Arial"
        mpl.rcParams['xtick.labelsize'] = 6
        mpl.rcParams['ytick.labelsize'] = 6
        mpl.rcParams['legend.fontsize'] = 8

        triv = dictGenes[gene][0]
        filename = output_name + '_' + gene + ' ' + triv
        figuretitle = gene + ' - ' + triv

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

        x_label = 'position along ORF [codon / aa]'
        y_label = 'enrichment [a.u.]'
        ax.set_xlabel(x_label, fontsize = 8)
        ax.set_ylabel(y_label, fontsize = 8)
        gene_length = int(dictGenes[gene][7])
        ax.set_xlim(0, gene_length) 
        ax.set_title(figuretitle, fontsize = 10)

        x = np.arange(gene_length)
        plt.plot(x, average, color = 'black', linewidth = 0.8)
        plt.fill_between(x, ratio_1, ratio_2, color = 'black', linewidth=0.0, alpha= 0.2)
        plt.yscale('log', basey=2)

        plt.savefig(output_path + filename + '.png', bbox_inches='tight', dpi = 300)
        plt.savefig(output_path + filename + '.pdf', bbox_inches='tight')
        plt.close()

    return 


if __name__ == '__main__':

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    file_total1 = 'sample name'                 # sample name total translatome 1 (no file extension)
    file_total2 = 'sample name'                 # sample name total translatome 2 (no file extension)
    file_selec1 = 'sample name'                 # sample name selective translatome 1 (no file extension)
    file_selec2 = 'sample name'                 # sample name selective translatome 2 (no file extension)
    output_name = 'experiment name'             # experiment name (no file extension)

    geneProfiles(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name)

