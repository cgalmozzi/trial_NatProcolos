#######################################
#                                     #
#    Footprint_Length_Distribution    #
#                                     #
#######################################

'''
Created on Mon Sep 24 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script compares the read length distribution of a given deep sequencing 
data set (e.g. ribosome profiling or selective ribosome profiling) and plots 
the result as histogram. 

The graph is saved as png file and as pdf file. 


'''

import matplotlib.pyplot as plt
import matplotlib as mpl


def LengthDist(input_path, sample_name): 

    # read input data
    dict_length = {}
    list_length = []
    
    with open(input_path + sample_name + '_footprint length.txt', 'r') as f:
        for line in f: 
            fields = line.split('\t')
            length = int(fields[0])
            number = int(fields[1])
            dict_length[length] = number
            if number != 0:
                list_length.append(length)

    # generate graph
    mpl.rcParams['axes.linewidth'] = 0.64
    mpl.rcParams['font.sans-serif'] = "Arial"
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 8

    x_min = min(list_length) - 2
    x_max = max(list_length) + 2
    x_label = ('footprint length [nt]')
    y_label = ('reads [*1,000,000]')
    x = list(range(x_min + 2, x_max - 2))
    y = [dict_length[i]/1000000 for i in x]

    fig = plt.figure(1, figsize = (3.2, 2.4))
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

    ax.set_xlim(x_min, x_max-1)
    ax.set_xlabel(x_label, fontsize = 8)
    ax.set_ylabel(y_label, fontsize = 8)

    ax.bar(x, y, align='center', color='#999999', width=0.6, linewidth=0)
    plt.axvspan(xmin=27.5, xmax=31.5, color='green', zorder = 0, alpha = 0.3, linewidth = 0)

    # save graph
    plt.savefig(input_path + sample_name + '_footprint length distribution.png', bbox_inches='tight', dpi = 600)
    plt.savefig(input_path + sample_name + '_footprint length distribution.pdf', bbox_inches='tight')
    plt.close()

    return


if __name__ == '__main__':

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    sample_name = 'sample name'                 # e.g. 'S1' (no filename extension)

    LengthDist(input_path, sample_name)
