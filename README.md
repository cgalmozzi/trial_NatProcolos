# trial_NatProcolos
This is trial to create a repository webpage for Galmozzi &amp; Merker etat
Data Analysis Based on Provided Example: 

-	all python scripts can be executed either via a standard python environment (IDLE, Sypder etc.) or via shell commands
-	for details on usage and help, please type “python Supplementary_script_XXX.py –h”

-	assuming the following data structure: 
|- sample1.sam
|- sample2.sam
|- sample3.sam
|- sample4.sam
|- Supplementary_script_A.py
|- Supplementary_script_B.py
|- Supplementary_script_C.py
|- Supplementary_script_D.py
|- Supplementary_script_E.py
|- Supplementary_script_F.py
|- Supplementary_script_G.py
|- references_yeast
	|- offset.txt
	|- yeast_genes.pkl
	|- yeast_introns.pkl
	|- yeast_sequence.pkl
	|- yeast_tRNA.pkl
-	the example reference files contain information on transcripts/sequence/tRNA of yeast’s chromosome 1 only 
-	sample1-4.sam are Ssb1 selective and total translatome samples from Döring et al. 2017 containing only reads mapping to yeast’s chromosome 1 
-	set current working directory to: main folder containing samples and scripts

-	command lines:
o	python Supplementary_script_A.py sample1.sam sample1
o	python Supplementary_script_A.py sample2.sam sample2
o	python Supplementary_script_A.py sample3.sam sample3
o	python Supplementary_script_A.py sample4.sam sample4
o	python Supplementary_script_B.py sample1
o	python Supplementary_script_B.py sample2
o	python Supplementary_script_B.py sample3
o	python Supplementary_script_B.py sample4
o	python Supplementary_script_C.py sample1
o	python Supplementary_script_C.py sample2
o	python Supplementary_script_C.py sample3
o	python Supplementary_script_C.py sample4
o	python Supplementary_script_D.py sample1 sample3 sample2 sample4 experiment1
o	python Supplementary_script_E.py sample1 sample3 sample2 sample4 experiment1
o	python Supplementary_script_F.py sample1 sample3 sample2 sample4 experiment1
expected output: 
Start of binding detection script
2019-01-21 20:14:41
Reading of all input files finished 
2019-01-21 20:15:53
Start of binding detection
2019-01-21 20:15:55
Included based on raw reads: 	63	transcripts
Background too high:			5 	transcripts
Correlation too low:			0 	transcripts
Identified Strong Binders:		7	transcripts
Identified Binders:			0 	transcripts
End of binding detection script
2019-01-21 20:15:55
o	python Supplementary_script_G.py sample1 sample3 sample2 sample4 experiment1
-	the expected output files are provided


