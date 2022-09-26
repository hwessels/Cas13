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



################################################################################
### MODULES ####################################################################

from optparse import OptionParser
import sys
from itertools import islice, chain
from Bio.Seq import Seq
from collections import defaultdict


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

def makedegenset(seq, mmalpha):
	degenset = set()
	for i in xrange(len(seq)):
		for mm in mmalpha:
			degenset.add(seq[:i] + mm + seq[i+1:])
	return degenset

def makebcdict(bclist, mmalpha):
	bcdict = {}
	for bc in bclist:
		if len(bc)>4:
			bcdegenset = makedegenset(bc, mmalpha)
			for mbc in bcdegenset:
				bcdict[mbc] = bc
		else:
			bcdict[bc] = bc
	return bcdict

def finddegenfeature(seq, degenset, flen, minpos, maxpos):
	for i in xrange(minpos, maxpos+1):
		pfeature = seq[i:i+flen]
		if pfeature in degenset:
			return i
	return -1

def getfeatures(seq, pos, offsets, flens):
	return (seq[pos+offsets[i]:pos+offsets[i]+flens[i]] for i in xrange(len(offsets)))

def bccorrect(bc, bcdict):
	if bc in bcdict:
		return bcdict[bc]
	return bc

def makedegenset_2mm(SeqList, mmalpha):
	for x in SeqList:
		degen = set()
		full = []
		for s in x:
			degen = makedegenset(s, mmalpha)	
			full.append(degen)
		ctrlsets_2mm = (set().union(*full))
		return ctrlsets_2mm


################################################################# /FUNCTIONS ###
################################################################################



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



################################################################################
### CONSTANTS ##################################################################

mmalpha = ['A', 'C', 'G', 'T', 'N']

# U6 Promoter and Cas13d Direct Repeat sequence
l1 = 'TCTTGTGGAAAGGACGAAACACCGAACCCCTACCAACTGGTCGGGGTTTGAAAC' 

l1len = len(l1)
l1set = makedegenset(l1, mmalpha)
l1List = []
l1List.append(list(l1set))
l1set_2mm = makedegenset_2mm(l1List, mmalpha)
l1List_2mm = []
l1List_2mm.append(list(l1set_2mm))
l1set_3mm = makedegenset_2mm(l1List_2mm, mmalpha)
l1List_3mm = []
l1List_3mm.append(list(l1set_3mm))
maxl1pos = 19
minl1pos = 9


l2= 'TTTTTGA'
l2len = len(l2)
l2set = makedegenset(l2, mmalpha)
l2List = []
l2List.append(list(l2set))
minl2pos = 21
maxl2pos = 25


l3= 'CAAGTA'
l3len = len(l3)
l3set = makedegenset(l3, mmalpha)
l3List = []
l3List.append(list(l3set))
minl3pos = 21
maxl3pos = 30


# offsets to retrieve adjacent sequences
# 1. i5 at -8
bcoffsets = [-8]
bclens = [8]

# guide start is actually at 58, but we trim the first nt
guideoffset = [54]
guideoffset2 = [1]

#max length of extension after removing the 5' anchor
extend = [40] 

# guide lengths is 19 or 20nt (20 when G is added for U6 start); However, we only record 17
guidelen = [23]




# in line i5 barcodes of fwd PCR2 primer
bclist_P5i = 	[ 'AAGTAGAG' , 'ACACGATC' , 'CGCGCGGT' , 'CATGATCG' , 'CGTTACCA' , 'TCCTTGGT' ,
				'AACGCATT' , 'ACAGGTAT' , 'AGGTAAGG' , 'AACAATGG' , 'ACTGTATC' , 'AGGTCGCA' ]



bcdict_P5i = makebcdict(bclist_P5i, mmalpha)


failseq = sum(bclens) * 'N'
failqual = sum(bclens) * '#'



################################################################# /CONSTANTS ###
################################################################################

ReadCnt = 0
PassCnt = 0
FailCnt = 0
BCFailCnt = 0
L1FailCnt = 0
L2FailCnt = 0
P5i_FailCnt = 0
U6count = 0
NoTermCount = 0
NoInsertCount = 0
ArrayCount = 0
TermCount = 0

# count barcode fail occurences
i5_set = set()
i5_faildict = defaultdict(int)
# count barcode match occurences
i5_dict = defaultdict(int)

RTbcCnt_out = open("Fwd_i5.txt", "w+")
RTbcCnt_failout = open("U6match_BCfail.txt", "w+")
U6_match_but_BCfail_out = open("U6match_BCfail.fastq", "w+") 
TerminatorFail_out = open("TerminatorFail.fastq", "w+")


################################################################################
### MAIN #######################################################################

if __name__ == "__main__":

	fastq = fqreader(opt.input_file)
	for rid, seq, qual in fastq:
		
		ReadCnt += 1

		# l1pos =  finddegenfeature(seq, l1set, l1len, minl1pos, maxl1pos)
		l1pos =  finddegenfeature(seq, l1set_3mm, l1len, minl1pos, maxl1pos)

		
		if l1pos >= 0:

			U6count += 1

			bcseqs = tuple(f for f in getfeatures(seq, l1pos, bcoffsets, bclens))
			
			bcfail_P5i = bcseqs[0] not in bcdict_P5i


			if bcfail_P5i:

				P5i_FailCnt += 1

				i5_set.add(bcseqs[0])
				i5_faildict[bcseqs[0]] += 1
				 
				concatBCseqs = ''.join(chain(bcseqs[0])) 

				guideSeq_gen = getfeatures(seq, l1pos, guideoffset, guidelen)
				guideQual_gen = getfeatures(qual, l1pos, guideoffset, guidelen)

				guideSeq = ''.join(chain(guideSeq_gen)) 
				guideQual = ''.join(chain(guideQual_gen)) 

				# Add barcode to the readID
				StripID=rid[1:]
				NewRid=''.join(['@',concatBCseqs,':',StripID])


				# keep failed Sequences
				U6_match_but_BCfail_out.write('\n'.join([NewRid, guideSeq, '+', guideQual]) + '\n')



			
			bcfail = bcfail_P5i
			

			if not bcfail:

				#get P5i barcode
				bcseqs_P5i = tuple(bccorrect(bc, bcdict_P5i) for bc in getfeatures(seq, l1pos, bcoffsets[0:1], bclens[0:1]))
				concatBCseqs = ''.join(chain(bcseqs_P5i)) 

				# Add barcode to the readID
				StripID=rid[1:]
				NewRid=''.join(['@',concatBCseqs,':',StripID])

				# remove 5' anchor sequence
				trimSeq_gen = getfeatures(seq, l1pos, guideoffset, extend)
				trimQual_gen = getfeatures(qual, l1pos, guideoffset, extend)
				trimSeq = ''.join(chain(trimSeq_gen)) 
				trimQual = ''.join(chain(trimQual_gen))

							

				# Test for presence of Terminator
				l2pos =  finddegenfeature(trimSeq, l2set, l2len, minl2pos, maxl2pos)
				# print l2pos	
				# if second anchor was not found correctly, count failure and report reads
				if l2pos >= 0:

					# count only if structure is correct
					i5_dict[bcseqs_P5i[0]] += 1

					# extract guide by trimming the 3' Terminator anchor
					guideSeq_gen = getfeatures(seq, l1pos, guideoffset, [l2pos])
					guideQual_gen = getfeatures(qual, l1pos, guideoffset, [l2pos])
					guideSeq = ''.join(chain(guideSeq_gen)) 
					guideQual = ''.join(chain(guideQual_gen)) 

					sys.stdout.write('\n'.join([NewRid, guideSeq, '+', guideQual]) + '\n')

					# quality score not needed here
					# bcqual_RTi = getfeatures(qual, l1pos, bcoffsets[0:1], bclens[0:1]) #would need to be reversed to match the rev comp of the RTi
					# bcqual_P5i = getfeatures(qual, l1pos, bcoffsets[1:2], bclens[1:2])
					# bcqual_P7i = getfeatures(qual, l1pos, bcoffsets[2:3], bclens[2:3])
					# concatBCquals = ''.join(chain(bcqual_P5i,bcqual_RTi,bcqual_P7i))
					
					PassCnt += 1
					TermCount +=1

				else:

					

					# Test for presence of early Terminator (if no guide inserted)
					l2pos =  finddegenfeature(trimSeq, l2set, l2len, minl2pos-5, maxl2pos-5)

					if l2pos >= 0:
						NoInsertCount += 1
						guideSeq_gen = getfeatures(seq, l1pos, guideoffset, [l2pos])
						guideQual_gen = getfeatures(qual, l1pos, guideoffset, [l2pos])
						guideSeq = ''.join(chain(guideSeq_gen)) 
						guideQual = ''.join(chain(guideQual_gen)) 

						sys.stdout.write('\n'.join([NewRid, guideSeq, '+', guideQual]) + '\n')

						L2FailCnt += 1
						FailCnt += 1
						TermCount +=1

					else: 

						# Test for presence of Array
						l3pos =  finddegenfeature(trimSeq, l3set, l3len, minl3pos, maxl3pos)

						if l3pos >= 0:
							
							
							guideSeq_gen = getfeatures(seq, l1pos, guideoffset, [l3pos])
							guideQual_gen = getfeatures(qual, l1pos, guideoffset, [l3pos])
							guideSeq = ''.join(chain(guideSeq_gen)) 
							guideQual = ''.join(chain(guideQual_gen)) 

							TerminatorFail_out.write('\n'.join([NewRid, guideSeq, '+', guideQual]) + '\n')

							PassCnt += 1
							ArrayCount += 1

						else:

							# Could still contain circRNA guides if > 23nt
							TerminatorFail_out.write('\n'.join([NewRid, trimSeq, '+', trimQual]) + '\n')

							NoTermCount += 1
							L2FailCnt += 1
							FailCnt += 1

						

				

		else:

			# sys.stdout.write('\n'.join([rid, failseq, '+', failqual]) + '\n')

			L1FailCnt += 1
			FailCnt += 1
		
	

for bc in i5_set: # type: str
	RTbcCnt_failout.write('\t'.join((bc, str(i5_faildict.get(bc, 0)))))
	RTbcCnt_failout.write('\n')
	RTbcCnt_failout.flush()
for bc in bclist_P5i: # type: str
	RTbcCnt_out.write('\t'.join((bc, str(i5_dict.get(bc, 0)))))
	RTbcCnt_out.write('\n')
	RTbcCnt_out.flush()


U6_match_but_BCfail_out.close()
TerminatorFail_out.close()


##################################################################### /MAIN ###
###############################################################################

sys.stderr.write("\n\n########################################################################\n")
sys.stderr.write("#  Barcode Cleanup finished.  \n")
sys.stderr.write(''.join(['#        ',str(ReadCnt),' reads have been processed']) + '\n')
sys.stderr.write(''.join(['#        ',str(U6count),' reads have the expected U6-DR anchor sequence']) + '\n')
sys.stderr.write(''.join(['#        ',str(PassCnt),' reads passed the barcode filter']) + '\n')
sys.stderr.write(''.join(['#        ',str(FailCnt),' reads did not pass the barcode or anchor sequence filters']) + '\n')
sys.stderr.write(''.join(['#        ',str(L1FailCnt),' reads did not pass the 5\' anchor sequence filter']) + '\n')
sys.stderr.write(''.join(['#        ',str(L2FailCnt),' reads did not pass the 3\' anchor sequence filter']) + '\n')
sys.stderr.write(''.join(['#        ',str(NoTermCount),' reads have no terminator or array sequence']) + '\n')
sys.stderr.write(''.join(['#        ',str(NoInsertCount),' have no gRNA insert']) + '\n')
sys.stderr.write(''.join(['#        ',str(TermCount),' have a terminator sequence']) + '\n')
sys.stderr.write(''.join(['#        ',str(ArrayCount),' have circRNA gRNA insert']) + '\n')
sys.stderr.write(''.join(['#        ',str(P5i_FailCnt),' P5i barcodes could not be assigned']) + '\n')
sys.stderr.write("########################################################################\n\n")
