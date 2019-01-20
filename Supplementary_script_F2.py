
# Binding_Detection    
                         

import os
import numpy
import sys
import pickle
from time import localtime, strftime


def writeTempFiles(input_path, sample_name): 
    ''' This module generates the required text files based on the given 
    dictionary with reads per position for each data set.
    '''
    
    # generate output folder
    os.mkdir(input_path + sample_name + '_temp')
    output_path = input_path + sample_name + '_temp/'

    # upload input data    
    dictReads = pickle.load(open(input_path + sample_name + '_reads.pkl', 'rb'))

    # write output data
    listNorm = list(dictReads.keys())
    listNorm.sort()
    for chrom in range(1, 17):
        out_p = open(output_path + 'NormReads_Chr_' + str(chrom).zfill(2) + '_plus.txt', 'w')
        out_m = open(output_path + 'NormReads_Chr_' + str(chrom).zfill(2) + '_minus.txt', 'w')
        out_p.write('track type = Wiggle_0 \nvariableStep chrom=chr' + str(chrom).zfill(2) + '\n')
        out_m.write('track type = Wiggle_0 \nvariableStep chrom=chr' + str(chrom).zfill(2) + '\n')
        for pos in listNorm: 
            if pos[:2] == str(chrom).zfill(2):
                out_p.write(str(int(pos[3:])) + '\t' + str(dictReads[pos][6]) + '\n')
                out_m.write(str(int(pos[3:])) + '\t' + str(dictReads[pos][8]) + '\n')
        out_p.close()
        out_m.close()     

    return


def BindingDetection(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name): 

    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = path_current + '/references_yeast/'
#    dictGenes = pickle.load(open(path_ref + 'yeast_genes.pkl', 'rb')) #important: change name if used!! 
    
    # print used python version and start time
    print ('Running Python ',sys.version_info)
    print ('Start of peak detection script')
    print (strftime("%Y-%m-%d %H:%M:%S", localtime()))
    
    # read in genes with introns:
    inFile = open(path_ref + 'yeast_genes_with_introns.txt', 'r')
    dictIntrons = {}
    
    for line in inFile:
        fields = line.split()
        length = len(fields)
        gene = fields[0]
        dictIntrons[gene] = []
        for i in range(1, length):
            elem = fields[i]
            dictIntrons[gene].append(elem)
        
    # read in gene expression files to exclude genes that have less than 64 reads from the downstream analysis; write all genes that are
    # included in both samples into a list:

    inFile2 = open(input_path + file_selec1 + '_gene expression.txt', 'r')
    inFile3 = open(input_path + file_selec2 + '_gene expression.txt', 'r')
    inFile9 = open(input_path + file_total1 + '_gene expression.txt', 'r')
    inFile10 = open(input_path + file_total2 + '_gene expression.txt', 'r')
   
    listIncluded = [] 
    dictIncluded = {}
   
    for line in inFile2:
        fields = line.split()
        name = fields[0]
        status = fields[4].strip()        
        RPKM = float(fields[3])
        if status == 'included':
            listIncluded.append(name)
            dictIncluded[name] = [RPKM,'-', '-','-']
    
    for line in inFile3:
        fields = line.split()
        name = fields[0]
        status = fields[4].strip()        
        RPKM = float(fields[3])
        if status == 'excluded' and name in listIncluded:
            listIncluded.remove(name)
            del dictIncluded[name]
        if name in dictIncluded:
            dictIncluded[name][1] = RPKM
            
    for line in inFile9:
        fields = line.split()
        name = fields[0]
        status = fields[4].strip()        
        RPKM = float(fields[3])

        if name in dictIncluded and RPKM >= minRPKM:
            dictIncluded[name][2] = RPKM
        if name in dictIncluded and RPKM < minRPKM:
            del dictIncluded[name]
            listIncluded.remove(name)
    
    for line in inFile10:
        fields = line.split()
        name = fields[0]
        status = fields[4].strip()        
        RPKM = float(fields[3])

        if name in dictIncluded and RPKM >= minRPKM:
            dictIncluded[name][3] = RPKM
        if name in dictIncluded and RPKM < minRPKM:
            del dictIncluded[name]
            listIncluded.remove(name)
    
    inFile2.close()
    inFile3.close()
    inFile9.close()
    inFile10.close()
    
    print ('Included genes:' + '\t' + str(len(listIncluded)))

    dictReads = {}
    dictReads['rep1'] = {}
    dictReads['rep2'] = {}
    dictReads['rep1Trans'] = {}
    dictReads['rep2Trans'] = {}
    
    dictGenes = {}

    inFile1 = open(path_ref + 'yeast_cds_all.txt', 'r')        
    
    for line in inFile1:
        fields = line.split()
        name = fields[0]
        chrom = fields[1].zfill(2)
        start = fields[2]
        stop = fields[3]
        if name.startswith('Y') and name in listIncluded:
            if name in dictIntrons:
                try:
                    dictGenes[name] = [chrom, dictIntrons[name][0], dictIntrons[name][1], dictIntrons[name][2], dictIntrons[name][3], dictIntrons[name][4], dictIntrons[name][5]]
                except IndexError:
                    dictGenes[name] = [chrom, dictIntrons[name][0], dictIntrons[name][1], dictIntrons[name][2], dictIntrons[name][3]]
                    
                dictReads['rep1'][name] = []
                dictReads['rep2'][name] = []
                dictReads['rep1Trans'][name] = []
                dictReads['rep2Trans'][name] = []
            else:
                dictGenes[name] = [chrom, start, stop]
                dictReads['rep1'][name] = []
                dictReads['rep2'][name] = []
                dictReads['rep1Trans'][name] = []
                dictReads['rep2Trans'][name] = []
        else:
            continue
        
    inFile1.close()

    # generate required text files 
    writeTempFiles(input_path, file_total1)
    writeTempFiles(input_path, file_total2)
    writeTempFiles(input_path, file_selec1)
    writeTempFiles(input_path, file_selec2)

    # read and process generated text files       
    dictRep1 = {}
    dictRep1['plus'] = {}
    dictRep1['minus'] = {}
    
    dictRep2 = {}
    dictRep2['plus'] = {}
    dictRep2['minus'] = {}
    
    for chromosome in range(1,17):
        chromosome = str(chromosome).zfill(2)
        dictRep1['plus'][chromosome] = []
        inFile4 = open(input_path + file_selec1 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_plus.txt', 'r')
        inFile4.readline()
        inFile4.readline()
        for line in inFile4:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep1['plus'][chromosome].append(read)
            
        dictRep1['minus'][chromosome] = []
        inFile4 = open(input_path + file_selec1 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_minus.txt', 'r')
        inFile4.readline()
        inFile4.readline()
        for line in inFile4:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep1['minus'][chromosome].append(read)
            
        dictRep2['plus'][chromosome] = []
        inFile5 = open(input_path + file_selec2 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_plus.txt', 'r')
        inFile5.readline()
        inFile5.readline()
        for line in inFile5:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep2['plus'][chromosome].append(read)
            
        dictRep2['minus'][chromosome] = []
        inFile5 = open(input_path + file_selec2 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_minus.txt', 'r')
        inFile5.readline()
        inFile5.readline()
        for line in inFile5:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep2['minus'][chromosome].append(read)   
        
        inFile4.close()
        inFile5.close()

    dictRep1Trans = {}
    dictRep1Trans['plus'] = {}
    dictRep1Trans['minus'] = {}
    
    dictRep2Trans = {}
    dictRep2Trans['plus'] = {}
    dictRep2Trans['minus'] = {}
    
    for chromosome in range(1,17):
        chromosome = str(chromosome).zfill(2)
        dictRep1Trans['plus'][chromosome] = []
        inFile6 = open(input_path + file_total1 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_plus.txt', 'r')
        inFile6.readline()
        inFile6.readline()
        for line in inFile6:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep1Trans['plus'][chromosome].append(read)
            
        dictRep1Trans['minus'][chromosome] = []
        inFile6 = open(input_path + file_total1 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_minus.txt', 'r')
        inFile6.readline()
        inFile6.readline()
        for line in inFile6:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep1Trans['minus'][chromosome].append(read)
            
        dictRep2Trans['plus'][chromosome] = []
        inFile7 = open(input_path + file_total2 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_plus.txt', 'r')
        inFile7.readline()
        inFile7.readline()
        for line in inFile7:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep2Trans['plus'][chromosome].append(read)
            
        dictRep2Trans['minus'][chromosome] = []
        inFile7 = open(input_path + file_total2 + '_temp/NormReads_Chr_' + str(chrom).zfill(2) + '_minus.txt', 'r')
        inFile7.readline()
        inFile7.readline()
        for line in inFile7:
            fields = line.split()
            pos = fields[0]
            read = abs(float(fields[1]))
            dictRep2Trans['minus'][chromosome].append(read)
            
        inFile6.close()
        inFile7.close()
        
    print ('Reading of all input files finished')
    print (strftime("%Y-%m-%d %H:%M:%S", localtime()))     

    # for each gene that is included in both samples write reads as a list into the dictionary (dictReads) that was already created:
    for gene in dictGenes:
        if gene not in dictIntrons:
            chrom, start, stop = dictGenes[gene]
            start = int(start)
            stop = int(stop)
            if abs(stop - start) > 90:
                if gene[6] == 'W':
                    for pos in range(start, stop+1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                if gene[6] == 'C':
                    for pos in reversed(range(stop, start+1)):
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
            else:
                del dictReads['rep1'][gene]
                del dictReads['rep2'][gene]

        if gene in dictIntrons:
            columns = len(dictIntrons[gene])
            chrom = dictGenes[gene][0]
            if columns == 5:
                if gene[6] == 'W':
                    for pos in range(int(dictIntrons[gene][0]), int(dictIntrons[gene][1]) + 1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                    for pos in range(int(dictIntrons[gene][2]), int(dictIntrons[gene][3]) + 1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                if gene[6] == 'C':          
                    for pos in reversed(range(int(dictIntrons[gene][1]) , (int(dictIntrons[gene][0]) + 1))):    
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
                    for pos in reversed(range(int(dictIntrons[gene][3]), (int(dictIntrons[gene][2]) + 1))):
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
            if columns == 7:
                if gene[6] == 'W':
                    for pos in range(int(dictIntrons[gene][0]), int(dictIntrons[gene][1]) + 1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                    for pos in range(int(dictIntrons[gene][2]), int(dictIntrons[gene][3]) + 1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                    for pos in range(int(dictIntrons[gene][4]), int(dictIntrons[gene][5]) + 1):
                        dictReads['rep1'][gene].append(dictRep1['plus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['plus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['plus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['plus'][chrom][pos-1])
                if gene[6] == 'C':          
                    for pos in reversed(range(int(dictIntrons[gene][1]), (int(dictIntrons[gene][0]) + 1))):
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
                    for pos in reversed(range(int(dictIntrons[gene][3]), (int(dictIntrons[gene][2]) + 1))):
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
                    for pos in reversed(range(int(dictIntrons[gene][5]), (int(dictIntrons[gene][4]) + 1))):
                        dictReads['rep1'][gene].append(dictRep1['minus'][chrom][pos-1])
                        dictReads['rep2'][gene].append(dictRep2['minus'][chrom][pos-1])
                        dictReads['rep1Trans'][gene].append(dictRep1Trans['minus'][chrom][pos-1])
                        dictReads['rep2Trans'][gene].append(dictRep2Trans['minus'][chrom][pos-1])
         
    del dictRep1
    del dictRep2
    del dictRep1Trans
    del dictRep2Trans
                      
    dict90nt = {}     
    dictMax90nt = {}  
    dictMax2 = {}
        
    for gene in dictReads['rep1']:
        counter = 0
        counter0 = 0
        First90nt = 0
        dict90nt[gene] = ['-', '-']
        dictMax90nt[gene] = []
        for elem in dictReads['rep1'][gene]:
            counter += 1
            if counter < background_stretch:
                First90nt += elem
                dictMax90nt[gene].append(elem)
                if elem == 0:
                    counter0 += 1
            if counter == background_stretch:
                First90nt += elem
                dictMax90nt[gene].append(elem)
                if elem == 0:
                    counter0 += 1
                if First90nt > 0 and counter0 < background_stretch/2:
                    dict90nt[gene][0] = First90nt/background_stretch
                    max90 = max(dictMax90nt[gene])
                    dictMax90nt[gene] = [max90]
                    break
                if First90nt == 0:
                    dict90nt[gene][0] = dictIncluded[gene][2]/1000              
                    dictMax90nt[gene] = [dictIncluded[gene][2]/1000]
                    break
                if First90nt != 0 and counter0 >= background_stretch/2:
                    if dictIncluded[gene][2]/1000 > First90nt/background_stretch:
                        dict90nt[gene][0] = dictIncluded[gene][2]/1000 
                        max90 = max(dictMax90nt[gene])
                        dictMax90nt[gene] = [max90]
                        break
                    else:
                        dict90nt[gene][0] =  First90nt/background_stretch
                        max90 = max(dictMax90nt[gene])
                        dictMax90nt[gene] = [max90]
                        break
    
    for gene in dictReads['rep1']:
        dictMax2[gene] = []
        counter2 = 0
        counter0 = 0
        First90nt_2 = 0
        for elem in dictReads['rep2'][gene]:
            counter2 += 1
            if counter2 < background_stretch:
                First90nt_2 += elem
                dictMax2[gene].append(elem)
                if elem == 0:
                    counter0 += 1
            if counter2 == background_stretch:
                First90nt_2 += elem
                dictMax2[gene].append(elem)
                if elem == 0:
                    counter0 += 1
                if First90nt_2 > 0 and counter0 < background_stretch/2:
                    dict90nt[gene][1] = First90nt_2/background_stretch
                    max90 = max(dictMax2[gene])
                    dictMax90nt[gene].append(max90)
                    break
                if First90nt_2 == 0:
                    dict90nt[gene][1] = dictIncluded[gene][3]/1000 
                    dictMax90nt[gene].append(dictIncluded[gene][3]/1000)
                    break
                if First90nt_2 != 0 and counter0 >= background_stretch/2:
                    if dictIncluded[gene][3]/1000 > First90nt_2/background_stretch:
                        dict90nt[gene][1] = dictIncluded[gene][3]/1000 
                        max90 = max(dictMax2[gene])
                        dictMax90nt[gene].append(max90)
                        break
                    else:
                        dict90nt[gene][1] = First90nt_2/background_stretch
                        max90 = max(dictMax2[gene])
                        dictMax90nt[gene].append(max90)
                        break
            
    dictStrong = {}
    dictStrong['rep1'] = {}
    dictStrong['rep2'] = {} 

    dictBinders = {}
    dictBinders['rep1'] = {}
    dictBinders['rep2'] = {}   
    deleteGenes = []
    deleteGenesStrong = []
            
    for gene in dictReads['rep1']:
        length = len(dictReads['rep1'][gene])
        dictStrong['rep1'][gene] = []
        dictStrong['rep2'][gene] = []
        dictBinders['rep1'][gene] = []
        dictBinders['rep2'][gene] = []
        for index in range(0, length):
            try:
                dictStrong['rep1'][gene].append(dictReads['rep1'][gene][index]/dictReads['rep1Trans'][gene][index])
            except ZeroDivisionError:       
                dictStrong['rep1'][gene].append(dictReads['rep1'][gene][index]/float(dict90nt[gene][0]))
                
            try:
                dictStrong['rep2'][gene].append(dictReads['rep2'][gene][index]/dictReads['rep2Trans'][gene][index])
            except ZeroDivisionError:       
                dictStrong['rep2'][gene].append(dictReads['rep2'][gene][index]/float(dict90nt[gene][1]))
            
            dictBinders['rep1'][gene].append(dictReads['rep1'][gene][index]/float(dict90nt[gene][0]))
            dictBinders['rep2'][gene].append(dictReads['rep2'][gene][index]/float(dict90nt[gene][1]))
            
        if max(dictReads['rep1'][gene][background_stretch:]) < minHeightBinders*dictMax90nt[gene][0]:
            if gene not in deleteGenes:
                deleteGenes.append(gene)
        if max(dictReads['rep2'][gene][background_stretch:]) < minHeightBinders*dictMax90nt[gene][1]:
            if gene not in deleteGenes:
                deleteGenes.append(gene)
        if max(dictReads['rep1'][gene][background_stretch:]) < minHeightStrongBinders*dictMax90nt[gene][0]:
            if gene not in deleteGenesStrong:
                deleteGenesStrong.append(gene)
        if max(dictReads['rep2'][gene][background_stretch:]) < minHeightStrongBinders*dictMax90nt[gene][1]:
            if gene not in deleteGenesStrong:
                deleteGenesStrong.append(gene)
                
    print (str(len(deleteGenes)) + 'background_too_high')
   
    for gene in deleteGenes:
        del dictBinders['rep1'][gene]
        del dictBinders['rep2'][gene]
        
    for gene in deleteGenesStrong:
        del dictStrong['rep1'][gene]
        del dictStrong['rep2'][gene]  

    outFile3 = open(input_path + output_name + '_BD incl genes.txt', 'w')
    for gene in dictIncluded:
        if gene in deleteGenes:
            outFile3.write(gene + '\t' + str(dictIncluded[gene][0]) + '\t' + str(dictIncluded[gene][1])+ '\t' + str(dictIncluded[gene][2]) + '\t' + str(dictIncluded[gene][3]) + '\t' + 'excluded:max too low' + '\n')
        else:
            outFile3.write(gene + '\t' + str(dictIncluded[gene][0]) + '\t' + str(dictIncluded[gene][1]) + '\t' + str(dictIncluded[gene][2]) + '\t' + str(dictIncluded[gene][3]) + '\n')
            

    dictPeaks = {}          
    dict2fold = {}          
    
    dict2fold['rep1'] = {}
    dict2fold['rep2'] = {}
    dictPeaks['rep1'] = {}
    dictPeaks['rep2'] = {}

    #strong binders

    print ('Start of peak detection')
    print (strftime("%Y-%m-%d %H:%M:%S", localtime()))

    for gene in dictStrong['rep1']:
        length = len(dictStrong['rep1'][gene])
        dictPeaks['rep1'][gene] = {}
        dictPeaks['rep2'][gene] = {}
        dict2fold['rep1'][gene] = {}
        dict2fold['rep2'][gene] = {}
        counter = 0
        switch = 0
        for index in range(background_stretch, length+1):
            switch = 0
            ratio = float(dictStrong['rep1'][gene][index-1])
            if ratio >= threshold_strongBinders:
                counter += 1
                dict2fold['rep1'][gene][index] = ratio
                if index == length:
                    ratio = 1
                    counter -= 1
            if ratio < threshold_strongBinders:
                if counter >= peakwidth:
                    start = index-counter                         
                    dictPeaks['rep1'][gene][start] = []
                    dictPeaks['rep2'][gene][start] = []
                    counter2 = 0
                    counter = 0
                    for pos in range(start, index):
                        dictPeaks['rep1'][gene][start].append(dict2fold['rep1'][gene][pos])
                        dict2fold['rep2'][gene][pos] = float(dictStrong['rep2'][gene][pos-1])
                    for pos in range(start, index):
                        if switch == 1:
                            break
                        else:
                            ratio2 = float(dictStrong['rep2'][gene][pos-1])
                            if ratio2 >= threshold_strongBinders:
                                counter2 += 1
                                if pos == index-1:
                                    ratio2 = 1
                            if ratio2 < threshold_strongBinders and counter2 < peakwidth/2:
                                counter2 = 0
                                continue
                            if ratio2 < threshold_strongBinders and counter2 > peakwidth/2:
                                if counter2 >= peakwidth:
                                    for position in range(start, pos+1):
                                        dictPeaks['rep2'][gene][start].append(dict2fold['rep2'][gene][position])
                                    counter2 = 0
                                if peakwidth > counter2 >= peakwidth/2:         #peaks have to overlap for a minimum of half the minimal peakwidth
                                    if pos == index-1 and float(dictStrong['rep2'][gene][index-2]) >= threshold_strongBinders:
                                        counter4 = counter2
                                        for shifted2 in range(index, (index+peakwidth)):
                                            if shifted2 < length:
                                                ratio4 = float(dictStrong['rep2'][gene][shifted2-1])
                                                if ratio4 >= threshold_strongBinders:
                                                    dict2fold['rep2'][gene][shifted2] = ratio4
                                                    counter4 += 1
                                                    if shifted2 == index+peakwidth-1:
                                                        ratio4 = 1
                                                if ratio4 < threshold_strongBinders:
                                                    if counter4 >= peakwidth:
                                                        for position2 in range(start, (index-counter2+counter4)):
                                                            dictPeaks['rep2'][gene][start].append(dict2fold['rep2'][gene][position2])
                                                    break
                                            else:
                                                if counter4 >= peakwidth:
                                                    for position2 in range(start, (index-counter2+counter4)):
                                                        dictPeaks['rep2'][gene][start].append(dict2fold['rep2'][gene][position2])
                                                break
                                    elif float(dictStrong['rep2'][gene][(start-1)]) >= threshold_strongBinders and float(dictStrong['rep2'][gene][(start-2)]) >= threshold_strongBinders:
                                        counter3 = counter2
                                        for shifted1 in reversed(range(start-peakwidth, start)):
                                            ratio3 = float(dictStrong['rep2'][gene][(shifted1-1)])
                                            if ratio3 >= threshold_strongBinders:
                                                dict2fold['rep2'][gene][shifted1] = ratio3
                                                dict2fold['rep1'][gene][shifted1] = float(dictStrong['rep1'][gene][(shifted1-1)])
                                                counter3 += 1
                                                if shifted1 == start-peakwidth:
                                                    ratio3 = 1
                                            if ratio3 < threshold_strongBinders:
                                                if counter3 >= peakwidth:
                                                    start2 = shifted1 + 1
                                                    for position1 in range(start2, index):
                                                        if start2 in dictPeaks['rep2'][gene]:
                                                            dictPeaks['rep2'][gene][start2].append(dict2fold['rep2'][gene][position1])
                                                            dictPeaks['rep1'][gene][start2].append(dict2fold['rep1'][gene][position1])
                                                        else:
                                                            dictPeaks['rep2'][gene][start2] = []
                                                            dictPeaks['rep1'][gene][start2] = []
                                                            dictPeaks['rep2'][gene][start2].append(dict2fold['rep2'][gene][position1])
                                                            dictPeaks['rep1'][gene][start2].append(dict2fold['rep1'][gene][position1])
                                                            del dictPeaks['rep2'][gene][start]
                                                            del dictPeaks['rep1'][gene][start]
                                                            counter = 0
                                                            counter2 = 0
                                                            switch = 1
                                                    break
                                                else:
                                                    switch = 1
                                                    counter2 = 0
                                                    break
                                    else:
                                        break
                                if counter2 < peakwidth/2:
                                    break
                else:
                    counter = 0
                    continue

                              
    EmptyPeaks = {}
    
    for gene in dictPeaks['rep1']:
        for peaks in dictPeaks['rep1'][gene]:
            if len(dictPeaks['rep1'][gene][peaks]) < peakwidth or len(dictPeaks['rep2'][gene][peaks]) < peakwidth:
                if gene in EmptyPeaks:
                    EmptyPeaks[gene].append(peaks)
                else:
                    EmptyPeaks[gene] = [peaks]

       
    for gene in EmptyPeaks:
        for peaks in EmptyPeaks[gene]:
            try:
                del dictPeaks['rep1'][gene][peaks]
            except KeyError:
                continue
            try:
                del dictPeaks['rep2'][gene][peaks]
            except KeyError:
                continue
            
    EmptyList = []
                    
    for gene in dictPeaks['rep1']:
        if len(dictPeaks['rep1'][gene]) == 0:
            EmptyList.append(gene)
        elif len(dictPeaks['rep2'][gene]) == 0:
            EmptyList.append(gene)
    
    for name in EmptyList:
        del dictPeaks['rep1'][name]
        del dictPeaks['rep2'][name]

            
    dictIncrease = {}          
    dict2foldWeak = {}         
    
    dict2foldWeak['rep1'] = {}
    dict2foldWeak['rep2'] = {}
    dictIncrease['rep1'] = {}
    dictIncrease['rep2'] = {}
                       
    for gene in dictBinders['rep1']:
        if gene not in dictPeaks['rep1']:
            averageRatio1 = numpy.mean(dictStrong['rep1'][gene])
            averageRatio2 = numpy.mean(dictStrong['rep2'][gene])
            switch = 0
            length = len(dictBinders['rep1'][gene])
            dictIncrease['rep1'][gene] = {}
            dictIncrease['rep2'][gene] = {}
            dict2foldWeak['rep1'][gene] = {}
            dict2foldWeak['rep2'][gene] = {}
            counterWeak = 0
            for index in range(background_stretch, length+1):
                switch = 0
                ratio = float(dictBinders['rep1'][gene][index-1])
                if ratio >= threshold_Binders and dictStrong['rep1'][gene][index-1] >= threshold_average*averageRatio1:
                    counterWeak += 1
                    dict2foldWeak['rep1'][gene][index] = ratio
                    if index == length:
                        ratio = 0
                        counterWeak -= 1
                if ratio < threshold_Binders or dictStrong['rep1'][gene][index-1] < threshold_average*averageRatio1:
                    if counterWeak >= peakwidth:                           
                        start = index-counterWeak                          
                        dictIncrease['rep1'][gene][start] = []
                        dictIncrease['rep2'][gene][start] = []
                        counterWeak2 = 0
                        counterWeak = 0
                        for pos in range(start, index):
                            dictIncrease['rep1'][gene][start].append(dict2foldWeak['rep1'][gene][pos])
                            dict2foldWeak['rep2'][gene][pos] = float(dictBinders['rep2'][gene][pos-1])
                        for pos in range(start, index):
                            if switch == 1:
                                break
                            else:
                                ratio2 = float(dictBinders['rep2'][gene][pos-1])
                                if ratio2 >= threshold_Binders and dictStrong['rep2'][gene][pos-1] >= threshold_average*averageRatio2:
                                    counterWeak2 += 1
                                    if pos == index-1:
                                        ratio2 = 0
                                if ratio2 < threshold_Binders and counterWeak2 < peakwidth/2:
                                    counterWeak2 = 0
                                    continue
                                if (ratio2 < threshold_Binders and counterWeak2 > peakwidth/2) or (dictStrong['rep2'][gene][pos-1] < threshold_average*averageRatio2 and counterWeak2 > 0):
                                    if counterWeak2 >= peakwidth:
                                        for position in range(start, pos+1):
                                            dictIncrease['rep2'][gene][start].append(dict2foldWeak['rep2'][gene][position])
                                        counterWeak2 = 0
                                    if peakwidth > counterWeak2 >= peakwidth/2:             #peaks have to overlap for a minimum of half the minimal peakwidth
                                        if pos == index-1 and float(dictBinders['rep2'][gene][index-2]) >= threshold_Binders:
                                            counterWeak4 = counterWeak2
                                            for shifted2 in range(index, (index+peakwidth)):
                                                if shifted2 < length:
                                                    ratio4 = float(dictBinders['rep2'][gene][shifted2-1])
                                                    if ratio4 >= threshold_Binders and dictStrong['rep2'][gene][shifted2-1] >= threshold_average*averageRatio2:
                                                        dict2foldWeak['rep2'][gene][shifted2] = ratio4
                                                        counterWeak4 += 1
                                                        if shifted2 == index+peakwidth-1:
                                                            ratio4 = 0
                                                    if ratio4 < threshold_Binders or dictStrong['rep2'][gene][shifted2-1] < threshold_average*averageRatio2:
                                                        if counterWeak4 >= peakwidth:
                                                            for position2 in range(start, (index-counterWeak2+counterWeak4)):
                                                                dictIncrease['rep2'][gene][start].append(dict2foldWeak['rep2'][gene][position2])
                                                        break
                                                else:
                                                    if counterWeak4 >= peakwidth:
                                                        for position2 in range(start, (index-counterWeak2+counterWeak4)):
                                                            dictIncrease['rep2'][gene][start].append(dict2foldWeak['rep2'][gene][position2])
                                                    break
                                        elif float(dictBinders['rep2'][gene][start-1]) >= threshold_Binders:
                                            counterWeak3 = counterWeak2
                                            for shifted1 in reversed(range(start-peakwidth, start)):
                                                ratio3 = float(dictBinders['rep2'][gene][shifted1-1])
                                                if ratio3 >= threshold_Binders and dictStrong['rep2'][gene][shifted1-1] >= threshold_average*averageRatio2:
                                                    counterWeak3 += 1
                                                    dict2foldWeak['rep2'][gene][shifted1] = ratio3
                                                    dict2foldWeak['rep1'][gene][shifted1] = float(dictBinders['rep1'][gene][shifted1-1])
                                                    if shifted1 == start-peakwidth:
                                                        ratio3 = 0
                                                if ratio3 < threshold_Binders or dictStrong['rep2'][gene][shifted1-1] < threshold_average*averageRatio2:
                                                    if counterWeak3 >= peakwidth:
                                                        start2 = shifted1 + 1
                                                        for position1 in range(start2, index):
                                                            if start2 in dictIncrease['rep2'][gene]:
                                                                dictIncrease['rep2'][gene][start2].append(dict2foldWeak['rep2'][gene][position1])
                                                                dictIncrease['rep1'][gene][start2].append(dict2foldWeak['rep1'][gene][position1])
                                                            else:
                                                                dictIncrease['rep2'][gene][start2] = []
                                                                dictIncrease['rep1'][gene][start2] = []
                                                                dictIncrease['rep2'][gene][start2].append(dict2foldWeak['rep2'][gene][position1])
                                                                dictIncrease['rep1'][gene][start2].append(dict2foldWeak['rep1'][gene][position1])
                                                                del dictIncrease['rep2'][gene][start]
                                                                del dictIncrease['rep1'][gene][start]
                                                                counterWeak = 0
                                                                counterWeak2 = 0
                                                                switch = 1
                                                        break
                                                    else:
                                                        switch = 1
                                                        counterWeak2 = 0
                                                        break
                                        
                                        else:
                                            break  
                                        
                                    if counterWeak2 < peakwidth/2:
                                        break
                    else:
                        counterWeak = 0
                        counterWeak2 = 0
                        continue
                    
    del dict2fold
    del dict2foldWeak

        
    EmptyPeaks2 = {} 
    
    for gene in dictIncrease['rep1']:
        for peaks in dictIncrease['rep1'][gene]:
            if len(dictIncrease['rep2'][gene][peaks]) < peakwidth or len(dictIncrease['rep1'][gene][peaks]) < peakwidth:
                if gene in EmptyPeaks2:
                    EmptyPeaks2[gene].append(peaks)
                else:
                    EmptyPeaks2[gene] = [peaks]

    for gene in EmptyPeaks2:
        for peaks in EmptyPeaks2[gene]:
            del dictIncrease['rep1'][gene][peaks]
            del dictIncrease['rep2'][gene][peaks]
            
    EmptyList2 = []
                    
    for gene in dictIncrease['rep1']:
        if len(dictIncrease['rep1'][gene]) == 0:
            EmptyList2.append(gene)
        elif len(dictIncrease['rep2'][gene]) == 0:
            EmptyList2.append(gene)
    
    for name in EmptyList2:
        del dictIncrease['rep1'][name]
        del dictIncrease['rep2'][name]
            
    print ('End of peak detection')
    print ('Strong binders:' + '\t' + str(len(dictPeaks['rep1'])))
    print ('Binders:' + '\t' + str(len(dictIncrease['rep1'])))
    print (strftime("%Y-%m-%d %H:%M:%S", localtime()))
            
               
    dictCorrelation = {}
    
    for name in dictReads['rep1']:
        correlation = numpy.corrcoef(dictReads['rep1'][name], dictReads['rep2'][name])[0, 1]
        dictCorrelation[name] = correlation
        

    outFile1 = open(input_path + output_name + '_BD strong binders.txt', 'w')
    outFile2 = open(input_path + output_name + '_BD binders.txt', 'w')
    
    outFile1.write('systematic name' + '\t' + 'pearson correlation' + '\t' + 'average of first 90 nt' + '\t' + 'number of peaks' + '\t' + 'peakwidth' + '\t' + 'start of peak [nt]' + '\t' + 'maximal peak height' + '\t' + 'average peak height' + '\t' + '\n' )
    outFile2.write('systematic name' + '\t' + 'pearson correlation' + '\t' + 'average of first 90 nt' + '\t' + 'number of peaks' + '\t' + 'peakwidth' + '\t' + 'start of peak [nt]' + '\t' + 'maximal peak height' + '\t' + 'average peak height' + '\t' + '\n' )
    
    counterStrongBinder = 0
    counterBinder = 0   
    counterCorrelation = 0
    
    for gene in dictPeaks['rep1']:
        number_peaks1 = len(dictPeaks['rep1'][gene])
        number_peaks2 = len(dictPeaks['rep2'][gene])
        if dictCorrelation[gene] > minCorrelation:
            counterStrongBinder += 1
            outFile1.write(gene + '\t' + str(dictCorrelation[gene]) + '\t' + str(dict90nt[gene][0]) + '\t' + str(number_peaks1) + '\t')
            for peaks in sorted(dictPeaks['rep1'][gene]):
                peak_length = len(dictPeaks['rep1'][gene][peaks])
                max_height = max(dictPeaks['rep1'][gene][peaks])
                average_height = numpy.mean(dictPeaks['rep1'][gene][peaks])
                outFile1.write(str(peak_length) + '\t' + str(peaks) + '\t' + str(max_height) + '\t' + str(average_height) + '\t' + '\t')
            outFile1.write('\n')
            outFile1.write(gene + '\t' + str(dictCorrelation[gene]) + '\t' + str(dict90nt[gene][1]) + '\t'+ str(number_peaks2) + '\t')
            for peaks in sorted(dictPeaks['rep2'][gene]):
                peak_length = len(dictPeaks['rep2'][gene][peaks])
                max_height2 = max(dictPeaks['rep2'][gene][peaks])
                average_height2 = numpy.mean(dictPeaks['rep2'][gene][peaks])
                outFile1.write(str(peak_length) + '\t' + str(peaks) + '\t' + str(max_height2) + '\t' + str(average_height2) + '\t' + '\t')
            outFile1.write('\n')
        else:
            counterCorrelation += 1
               
    for gene in dictIncrease['rep1']:
        number_peaks1 = len(dictIncrease['rep1'][gene])
        number_peaks2 = len(dictIncrease['rep2'][gene])
        if dictCorrelation[gene] > minCorrelation:
            counterBinder += 1
            outFile2.write(gene + '\t' + str(dictCorrelation[gene]) + '\t' + str(dict90nt[gene][0]) + '\t'+ str(number_peaks1) + '\t')
            for peaks in sorted(dictIncrease['rep1'][gene]):
                peak_length = len(dictIncrease['rep1'][gene][peaks])
                max_height3 = max(dictIncrease['rep1'][gene][peaks])
                average_height3 = numpy.mean(dictIncrease['rep1'][gene][peaks])
                outFile2.write(str(peak_length) + '\t' + str(peaks) + '\t' + str(max_height3) + '\t' + str(average_height3) + '\t' + '\t')
            outFile2.write('\n')
            outFile2.write(gene + '\t' + str(dictCorrelation[gene]) + '\t' + str(dict90nt[gene][1]) + '\t'+ str(number_peaks2) + '\t')
            for peaks in sorted(dictIncrease['rep2'][gene]):
                peak_length = len(dictIncrease['rep2'][gene][peaks])
                max_height4 = max(dictIncrease['rep2'][gene][peaks])
                average_height4 = numpy.mean(dictIncrease['rep2'][gene][peaks])
                outFile2.write(str(peak_length) + '\t' + str(peaks) + '\t' + str(max_height4) + '\t' + str(average_height4) + '\t' + '\t')
            outFile2.write('\n')
        else:
            counterCorrelation += 1
    
    print (str(counterCorrelation) + 'below0.5')
    
    print ('End of script')
    print ('Strong Binders:' + '\t' + str(counterStrongBinder))
    print ('Binders:' + '\t' + str(counterBinder))
    print (strftime("%Y-%m-%d %H:%M:%S", localtime()))
    
    return
    
if __name__ == '__main__':

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    file_total1 = 'sample name'                 # sample name total translatome 1 (no file extension)
    file_total2 = 'sample name'                 # sample name total translatome 2 (no file extension)
    file_selec1 = 'sample name'                 # sample name selective translatome 1 (no file extension)
    file_selec2 = 'sample name'                 # sample name selective translatome 2 (no file extension)
    output_name = 'experiment name'             # experiment name (no file extension)

    # parameters to adjust    
    peakwidth = 15                  # minimal width of peak in nt
    minRPKM = 8                     # minimal coverage of included genes in RPKM
    background_stretch = 90         # first 90nt of each gene (= not exposed from ribosomal tunnel)
    minHeightBinders = 2            # minimal difference between highest signal in background stretch and later signals for binders
    minHeightStrongBinders = 2      # minimal difference between highest signal in background stretch and later signals for strong binders
    threshold_strongBinders = 1.5   # minimal ratio (interactome/translatome) at each position of a peak for strong binders
    threshold_Binders = 3           # minimal ratio (interactome/background stretch) at each position of a peak for binders
    threshold_average = 1.5         # minimal ratio ((interactome/translatome)/average ratio (interactome/translatome) along gene): exclude detection of peaks due to stalling peak in translatome
    minCorrelation = 0.5            # minimal Pearson correlation between replicates

    BindingDetection(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name)
