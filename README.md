# get_interdistances
Python  : python 2.7.9

Matplotlib : 2.0.2

Numpy : 1.8.2

Biopython : 1.70

######################################################
Authors: Arnaud Stigliani, Adrien Bessy, Jeremy Lucas
######################################################

# How to use it 
Here is the command line with the mandatory arguments : 

$ python get_interdistances.py -pfm pfm.txt  -o output_directory -n figure_name.svg -maxInter window_size  -pos pos.fas -th threshold_1 threshold_2 threshold_n -neg neg_1.fas neg_2.fas neg_n.fas -points True  

################################################################################
Note : If you want to get a png, you can use the following command after the svg is generated:
$ inkscape -z -e figure_name.png -w 1800 -h 1024 figure_name.svg
################################################################################

The pfm (option -pfm) has to be in the following format :

68.0	182.0	65.0	289.0  
98.0	49.0	26.0	428.0  
1.0	1.0	44.0	557.0  
1.0	1.0	601.0	1.0  
1.0	1.0	1.0	601.0  
1.0	596.0	1.0	5.0  
2.0	4.0	463.0	134.0  
53.0	25.0	489.0	35.0  
131.0	173.0	140.0	158.0  
239.0	62.0	67.0	236.0  

The matrix has NO HEADER. The columns respectively stand for A, C, G, T frequencies. The decimal part ".0" is important for the frequencies to be interpreted as floats (This bug will soon be fixed). Each frequency must be different from 0.

The threshold list maximum size is 6 (option -th)

The list of sequences that you want to study is given in the fasta format (option -pos)

# About the control sets of sequences

(Option -neg)

The program should work without control sets, computing only the configuration frequencies. Nevertheless, It has been changed a lot and this function has not been used for a long time (not sure it still works).

You can give as much control sets as you want. A nice way to generate control sets is the BiasAway  tool :
https://github.com/wassermanlab/BiasAway
R. Worsley-Hunt et al., Improving detection and analysis of transcription factor binding sites within ChIP-Seq

###############################################################################
Reference:

This tool has been described in 
STIGLIANI, Arnaud, MARTIN-AREVALILLO, Raquel, LUCAS, Jeremy, et al. Capturing auxin response factors syntax using DNA binding models. Molecular plant, 2018.
Please cite this publication when using this program.
