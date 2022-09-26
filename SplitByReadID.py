#!/usr/bin/env python 

################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2022) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Version: 0.1
# Author: Hans-Hermann Wessels

################################################################# /COPYRIGHT ###
################################################################################

# Split Fastq based on Barcode in ReadID

################################################################################
### MODULES ####################################################################

from optparse import OptionParser
import sys
from itertools import islice



################################################################### /MODULES ###
################################################################################

################################################################################
### FUNCTIONS ##################################################################

def fqreader(filename):
	if filename != 'stdin':
		fq = open(filename, 'r')
	else:
		fq = sys.stdin
	while True:
		r = [line.rstrip() for line in islice(fq, 4)]
		if not r:
			break
		yield r[0], r[1], r[3]
	if filename != 'stdin':
		fq.close()

################################################################################
### ARGUMENTS,OPTIONS ##########################################################

parser = OptionParser(usage="\n%prog [options]", version="%prog v0.1")

parser.add_option(
	"-i",
	metavar = "FILE",
	type = "string",
	dest = "input_file",
	default = 'stdin',
	help = "Input R1 FASTQ file (default = 'stdin')"
	)

(opt, args) = parser.parse_args()

######################################################### /ARGUMENTS,OPTIONS ###
################################################################################

# Expected Barcodes
BC = 	[ 'AAGTAGAG' , 'ACACGATC' , 'CGCGCGGT' , 'CATGATCG' , 'CGTTACCA' , 'TCCTTGGT' ,	'AACGCATT' , 'ACAGGTAT' , 'AGGTAAGG' , 'AACAATGG' , 'ACTGTATC' , 'AGGTCGCA' ]


# BCs =	{
#   "INPUT": 'ACAGGTAT',
#   "BIN1": 'AGGTAAGG',
#   "BIN2": 'AACAATGG',
#   "BIN3": 'ACTGTATC',
#   "BIN4": 'AGGTCGCA'
# }


B1_out = open("Fwd_01.fastq", "w+") 
#B2_out = open("Fwd_02.fastq", "w+") 
B3_out = open("Fwd_03.fastq", "w+") 
#B4_out = open("Fwd_04.fastq", "w+") 
B5_out = open("Fwd_05.fastq", "w+") 
#B6_out = open("Fwd_06.fastq", "w+") 
B7_out = open("Fwd_07.fastq", "w+") 
B8_out = open("Fwd_08.fastq", "w+") 
B9_out = open("Fwd_09.fastq", "w+") 
B10_out = open("Fwd_10.fastq", "w+") 
B11_out = open("Fwd_11.fastq", "w+") 
B12_out = open("Fwd_12.fastq", "w+") 
Undetermined_out = open("Unassigned.fastq", "w+") 


################################################################################
### MAIN #######################################################################

if __name__ == "__main__":

	fastq = fqreader(opt.input_file)
	for rid, seq, qual in fastq:

		bc = rid.split(':')[0]

		if str(bc[1:]) == str(BC[0]):

			B1_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[1]):

			B1_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[2]):

			B3_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[3]):

			B3_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[4]):

			B5_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[5]):

			B5_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[6]):

			B7_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[7]):

			B8_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[8]):

			B9_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[9]):

			B10_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[10]):

			B11_out.write('\n'.join([rid, seq, '+', qual]) + '\n')

		if str(bc[1:]) == str(BC[11]):

			B12_out.write('\n'.join([rid, seq, '+', qual]) + '\n')
		
		if str(bc[1:]) not in  BC:

			Undetermined_out.write('\n'.join([rid, seq, '+', qual]) + '\n')
		

B1_out.close()
#B2_out.close()
B3_out.close()
#B4_out.close()
B5_out.close()
#B6_out.close()
B7_out.close()
B8_out.close()
B9_out.close()
B10_out.close()
B11_out.close()
B12_out.close()
Undetermined_out.close()


##################################################################### /MAIN ###
###############################################################################
