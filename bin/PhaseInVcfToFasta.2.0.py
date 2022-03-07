#!/usr/bin/env python
#
#  Copyright 2014 CIRAD
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# -*- coding: utf-8 -*-
import sys
sys.stdout.write('loading modules\n')
import optparse
import os
import datetime
import gzip
from Bio import AlignIO

sys.stdout.write('modules loaded\n')

def countIUPAC(SEQ):
	total = len(SEQ)
	NonIUPAC = set("ATGCatgc")
	Val = 0
	for n in SEQ:
		if not(n in NonIUPAC):
			Val += 1
	return Val/total

def IdentIupacPOS(SEQ):
	NonIUPAC = set("ATGCatgc")
	setPos = set()
	for i in range(len(SEQ)):
		n = SEQ[i]
		if not(n in NonIUPAC):
			setPos.add(i)
	return setPos

def CodeIUPAC(GENO, ALLELE):
	
	TODECODE = set()
	for All in GENO:
		TODECODE.add(ALLELE[int(All)])
	
	if "A" in TODECODE or "a" in TODECODE:
		if "G" in TODECODE or "g" in TODECODE:
			if "C" in TODECODE or "c" in TODECODE:
				if "T" in TODECODE or "t" in TODECODE:
					return "N"
				else:
					return "V"
			elif "T" in TODECODE or "t" in TODECODE:
				return "D"
			else:
				return "R"
		elif "C" in TODECODE or "c" in TODECODE:
			if "T" in TODECODE or "t" in TODECODE:
				return "H"
			else:
				return "M"
		elif "T" in TODECODE or "t" in TODECODE:
			return "W"
		else:
			return "A"
	elif "G" in TODECODE or "g" in TODECODE:
		if "C" in TODECODE or "c" in TODECODE:
			if "T" in TODECODE or "t" in TODECODE:
				return "B"
			else:
				return "S"
		elif "T" in TODECODE or "t" in TODECODE:
			return "K"
		else:
			return "G"
	elif "C" in TODECODE or "c" in TODECODE:
		if "T" in TODECODE or "t" in TODECODE:
			return "Y"
		else:
			return "C"
	else:
		return "T"

def Phase(VCF, PCT, PREFIX, REGION, OUTGP, WIN, MAXIUP, MINLEN):
	
	"""
		Filter vcf file and output a filtered vcf
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param PCT: A dictionary containing with keys = F1 and values = [P1,P2]
		:type PCT: str
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
		:param REGION: Regions to work with
		:type REGION: str
		:param OUTGP: A path to a file containing accessions to add but for which no parent child trios are available
		:type OUTGP: Boolean
		:param WIN: Window size for sub fasta creation
		:type WIN: int
		:param MAXIUP: Maximal UIPAC letter proportion to keep the alignment
		:type MAXIUP: float
		:param MINLEN: Minimal size of the alignment to keep for analysis
		:type MINLEN: int
	"""
	
	# Initiating values
	CHROM, START, END = REGION.split(":")
	if WIN:
		window = int(START)
		while window < int(END):
			outfasta = open(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.fasta','w')
			outfasta.close()
			window += int(WIN)
	else:
		outfasta = open(PREFIX+'_Phase.fasta','w')
	DicoSeq = {}
	window = int(START) # For creating sub-fasta
	
	# Loading outgroup accessions
	file = open(OUTGP)
	dicoOUT = set()
	if OUTGP:
		for line in file:
			data = line.split()
			if data:
				DicoSeq[data[0]] = ""
				dicoOUT.add(data[0])
		file.close()
	
	# Fixing order in which accessions are treated
	F1ORDER = sorted(list(PCT.keys()))
	
	# Creating dictionary to count problems (impossible case, heterozygous, multiallelic)
	dicoImpossible = {}
	for n in F1ORDER:
		dicoImpossible[n] = [0,0,0]
	
	# Reading and phasing line by line
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
	PHYMLAnalysis = set()
	
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				Accession_start = header.index('FORMAT')+1
				ACCPOS = {}
				for n in header[Accession_start:]:
					ACCPOS[n] = header.index(n)
				for n in DicoSeq:
					ACCPOS[n] = header.index(n)
				FILTERpos = header.index('FILTER')
				REFpos = header.index('REF')
				ALTpos = header.index('ALT')
				FLAGformat = header.index('FORMAT')
				CHR = header.index('#CHROM')
				CHRpos = header.index('POS')
				
				for F1 in F1ORDER:
					# DicoSeq[''.join([F1,'-H1-','_X_'.join(PCT[F1])])] = ""
					# DicoSeq[''.join([F1,'-H2-','_X_'.join(PCT[F1])])] = ""
					DicoSeq[''.join([PCT[F1][0],'-H1-from-',F1])] = ""
					DicoSeq[''.join([PCT[F1][0],'-H2-from-',F1])] = ""
					DicoSeq[''.join([PCT[F1][1],'-H1-from-',F1])] = ""
					DicoSeq[''.join([PCT[F1][1],'-H2-from-',F1])] = ""
			
			# Printing headers
			elif data[0][0] == '#':
				pass
			# Working on variant
			else:
				if data[CHR] == CHROM:
					if int(START) <= int(data[CHRpos]) and int(END) >= int(data[CHRpos]):
						ref_allele = data[REFpos]
						alt_allele = data[ALTpos].split(',')
						AlleleList = [ref_allele]+alt_allele
						flag_format = data[FLAGformat].split(':')
						if 'AD' in flag_format and 'DP' in flag_format and 'GT' in flag_format:
							GTpos = flag_format.index('GT')
							DPpos = flag_format.index('DP')
							ADpos = flag_format.index('AD')
							
							for accession in F1ORDER:
								F1info = data[ACCPOS[accession]].split(':')
								P1info = data[ACCPOS[PCT[accession][0]]].split(':')
								P2info = data[ACCPOS[PCT[accession][1]]].split(':')
								F1GT = set(F1info[GTpos].split('/'))
								P1GT = set(P1info[GTpos].split('/'))
								P2GT = set(P2info[GTpos].split('/'))
								
								F1GTLi = list(F1info[GTpos].split('/'))
								P1GTLi = list(P1info[GTpos].split('/'))
								P2GTLi = list(P2info[GTpos].split('/'))
								OK = 0
								
								# Validating bi-allelic site
								CompleteSet = set(F1GTLi + P1GTLi + P2GTLi)
								if len(CompleteSet) > 2:
									dicoImpossible[accession][2] += 1
								
								# Looking for impossible case
								if len(F1GT) == 1: # The F1 is homozygous
									CommonAllSet = P1GT & P2GT
									if not(F1GTLi[0] in CommonAllSet):# impossible case
										# print('imp1', data[1], accession, F1GT, P1GT, P2GT)
										dicoImpossible[accession][0] += 1
									else:
										if len(CommonAllSet) == 0: # impossible case
											# print('imp2', data[1], accession, F1GT, P1GT, P2GT)
											dicoImpossible[accession][0] += 1
										else:
											CommonAll = list(CommonAllSet & F1GT)[0]
											F1L = int(CommonAll)
											F1R = int(CommonAll)
											P1L = int(CommonAll)
											P2L = int(CommonAll)
											P1R = P1GTLi[:]
											P2R = P2GTLi[:]
											P1R.remove(CommonAll)
											P2R.remove(CommonAll)
											P1R = int(P1R[0])
											P2R = int(P2R[0])
											OK = 1
								else:
									if len(P1GT) == 1:
										CommonP1F1 = P1GTLi[0]
										CommonP2F1Li = F1GTLi[:]
										CommonP2F1Li.remove(CommonP1F1)
										CommonP2F1 = CommonP2F1Li[0]
										if not(CommonP2F1 in P2GT):# impossible case
											# print('imp3', data[1], accession, F1GT, P1GT, P2GT)
											dicoImpossible[accession][0] += 1
										else:
											F1L = int(CommonP1F1)
											F1R = int(CommonP2F1)
											P1L = int(CommonP1F1)
											P2L = int(CommonP2F1)
											P1R = P1GTLi[:]
											P2R = P2GTLi[:]
											P1R.remove(CommonP1F1)
											P2R.remove(CommonP2F1)
											P1R = int(P1R[0])
											P2R = int(P2R[0])
											OK = 1
											
									elif len(P2GT) == 1:
										CommonP2F1 = P2GTLi[0]
										CommonP1F1Li = F1GTLi[:]
										CommonP1F1Li.remove(CommonP2F1)
										CommonP1F1 = CommonP1F1Li[0]
										if not(CommonP1F1 in P1GT):# impossible case
											# print('imp4', data[1], accession, F1GT, P1GT, P2GT)
											dicoImpossible[accession][0] += 1
										else:
											F1L = int(CommonP1F1)
											F1R = int(CommonP2F1)
											P1L = int(CommonP1F1)
											P2L = int(CommonP2F1)
											P1R = P1GTLi[:]
											P2R = P2GTLi[:]
											P1R.remove(CommonP1F1)
											P2R.remove(CommonP2F1)
											P1R = int(P1R[0])
											P2R = int(P2R[0])
											OK = 1
									else:
										# print('imp5', data[1], accession, F1GT, P1GT, P2GT)
										dicoImpossible[accession][1] += 1
								
								if WIN and int(data[CHRpos]) > window+int(WIN)-1:
									
									outfasta = open(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.fasta','w')
									
									# Verification of sequence size
									dicoNoATGC = set()
									SeqSize = 0
									for n in DicoSeq:
										SeqSize = max(SeqSize, len(DicoSeq[n]))
									for n in DicoSeq:
										if len(DicoSeq[n]) != SeqSize:
											sys.exit('Oups, there is a bug in the program\n')
									
									# Creating the fasta
									KeepForPhylo = True
									for n in DicoSeq:
										dicoNoATGC.update(IdentIupacPOS(DicoSeq[n]))
										outfasta.write('>'+n+'\n')
										outfasta.write(DicoSeq[n]+'\n')
										if countIUPAC(DicoSeq[n]) < MAXIUP:
											pass
										else:
											print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because too much IUPAC code in '+n+': '+str(countIUPAC(DicoSeq[n])))
											KeepForPhylo = False
									outfasta.close()
									
									if KeepForPhylo:
										if SeqSize >= MINLEN:
											PHYMLAnalysis.add('sparse -q agap_normal -N 1 -n 1 -m 10G -t 48:00:00 -j PHYML -c "phyml -d nt -m GTR -b -5 -v e -c 4 -a e -s BEST -i '+PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.phy"')
										else:
											print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because short alignment. Alignment length: '+str(SeqSize))
									
									outphylip = open(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.phy','w')
									record_dict = AlignIO.parse(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.fasta', "fasta")
									AlignIO.write(record_dict, outphylip, "phylip-relaxed")
									outphylip.close()
									
									
									while int(data[CHRpos]) > window+int(WIN)-1:
										window += int(WIN)
									for n in DicoSeq:
										DicoSeq[n] = ""
									
								if OK:
									# DicoSeq[''.join([accession,'-H1-','_X_'.join(PCT[accession])])] += AlleleList[F1L]
									# DicoSeq[''.join([accession,'-H2-','_X_'.join(PCT[accession])])] += AlleleList[F1R]
									DicoSeq[''.join([PCT[accession][0],'-H1-from-',accession])] += AlleleList[P1L]
									DicoSeq[''.join([PCT[accession][0],'-H2-from-',accession])] += AlleleList[P1R]
									DicoSeq[''.join([PCT[accession][1],'-H1-from-',accession])] += AlleleList[P2L]
									DicoSeq[''.join([PCT[accession][1],'-H2-from-',accession])] += AlleleList[P2R]
									
								else:
									# DicoSeq[''.join([accession,'-H1-','_X_'.join(PCT[accession])])] += CodeIUPAC(F1GT, AlleleList)
									# DicoSeq[''.join([accession,'-H2-','_X_'.join(PCT[accession])])] += CodeIUPAC(F1GT, AlleleList)
									DicoSeq[''.join([PCT[accession][0],'-H1-from-',accession])] += CodeIUPAC(P1GT, AlleleList)
									DicoSeq[''.join([PCT[accession][0],'-H2-from-',accession])] += CodeIUPAC(P1GT, AlleleList)
									DicoSeq[''.join([PCT[accession][1],'-H1-from-',accession])] += CodeIUPAC(P2GT, AlleleList)
									DicoSeq[''.join([PCT[accession][1],'-H2-from-',accession])] += CodeIUPAC(P2GT, AlleleList)
							
							for accession in dicoOUT:
								ACCinfo = data[ACCPOS[accession]].split(':')
								ACCGT = set(ACCinfo[GTpos].split('/'))
								if '.' in ACCGT:
									DicoSeq[accession] += "N"
								elif len(ACCGT) == 1:
									DicoSeq[accession] += AlleleList[int(list(ACCGT)[0])]
								else:
									DicoSeq[accession] += CodeIUPAC(ACCGT, AlleleList)
					elif int(END) <= int(data[CHRpos]):
						break
	if WIN:
		outfasta = open(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.fasta','w')
		
		# Verification of sequence size
		dicoNoATGC = set()
		SeqSize = 0
		for n in DicoSeq:
			SeqSize = max(SeqSize, len(DicoSeq[n]))
		for n in DicoSeq:
			if len(DicoSeq[n]) != SeqSize:
				sys.exit('Oups, there is a bug in the program\n')
		
		# Creating the fasta
		KeepForPhylo = True
		for n in DicoSeq:
			dicoNoATGC.update(IdentIupacPOS(DicoSeq[n]))
			outfasta.write('>'+n+'\n')
			outfasta.write(DicoSeq[n]+'\n')
			if countIUPAC(DicoSeq[n]) < MAXIUP:
				pass
			else:
				print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because too much IUPAC code in '+n+': '+str(countIUPAC(DicoSeq[n])))
				KeepForPhylo = False
		outfasta.close()
		
		if KeepForPhylo:
			if SeqSize >= MINLEN:
				PHYMLAnalysis.add('sparse -q agap_normal -N 1 -n 1 -m 10G -t 48:00:00 -j PHYML -c "phyml -d nt -m GTR -b -5 -v e -c 4 -a e -s BEST -i '+PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.phy"')
			else:
				print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because short alignment. Alignment length: '+str(SeqSize))
		
		outphylip = open(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.phy','w')
		record_dict = AlignIO.parse(PREFIX+'_Phase'+str(window)+'_'+str(window+int(WIN)-1)+'.fasta', "fasta")
		AlignIO.write(record_dict, outphylip, "phylip-relaxed")
		outphylip.close()
	else:
		
		# Verification of sequence size
		dicoNoATGC = set()
		SeqSize = 0
		for n in DicoSeq:
			SeqSize = max(SeqSize, len(DicoSeq[n]))
		for n in DicoSeq:
			if len(DicoSeq[n]) != SeqSize:
				sys.exit('Oups, there is a bug in the program\n')
		
		# Creating the fasta
		KeepForPhylo = True
		for n in DicoSeq:
			dicoNoATGC.update(IdentIupacPOS(DicoSeq[n]))
			outfasta.write('>'+n+'\n')
			outfasta.write(DicoSeq[n]+'\n')
			if countIUPAC(DicoSeq[n]) < MAXIUP:
				pass
			else:
				print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because too much IUPAC code in '+n+': '+str(countIUPAC(DicoSeq[n])))
				KeepForPhylo = False
		outfasta.close()
		
		if KeepForPhylo:
			if SeqSize >= MINLEN:
				PHYMLAnalysis.add('sparse -q agap_normal -N 1 -n 1 -m 10G -t 48:00:00 -j PHYML -c "phyml -d nt -m GTR -b -5 -v e -c 4 -a e -s BEST -i '+PREFIX+'_Phase.phy"')
			else:
				print('Not used region '+str(window)+' '+str(window+int(WIN)-1)+' because short alignment. Alignment length: '+str(SeqSize))
		
		outphylip = open(PREFIX+'_Phase.phy','w')
		record_dict = AlignIO.parse(PREFIX+'_Phase.fasta', "fasta")
		AlignIO.write(record_dict, outphylip, "phylip-relaxed")
		outphylip.close()
	
	# Running phyml analysis
	for n in PHYMLAnalysis:
		os.system(n)
	for n in PHYMLAnalysis:
		print(n)
	
	for n in dicoImpossible:
		print('Impossible','HetOnAll','NotBiallelic')
		print(n, dicoImpossible[n])

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '',	'--names',			dest='names',		default=None,			help='A 3 column file containing in this order F1 P1 P2. [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')
	parser.add_option( '',	'--region',			dest='region',		default=None, 			help='Region to work with. Region should be specified as follows : "chrX:start:end"')
	parser.add_option( '',	'--outgp',			dest='outgp',		default=False, 			help='Accessions for which parent child trio is not available but genotype should be obtained. Accessions should be passed in a one column file.')
	parser.add_option( '',	'--win',			dest='win',			default=False, 			help='Window size for each sub fasta')
	parser.add_option( '',	'--MinLength',		dest='MinLength',	default="1000", 		help='Minimal length of the alignment (in bp) to keep the alignment for PHYML.')
	parser.add_option( '',	'--maxIUP',			dest='maxIUP',		default="1",			help='Minimal proportion of IUPAC letters in a sequence to keep the sequence in the alignment. Values comprised between 0 and 1. [Default: %default]')
	
	(options, args) = parser.parse_args()
	
	# Filtering vcf file
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.names == None:
		sys.exit('Please provide a name file to --names argument')
	if options.region == None:
		sys.exit('Please provide argument to --region argument')
	
	# Recording information
	file = open(options.names)
	dicoPCT = {}
	for line in file:
		data = line.split()
		if data:
			dicoPCT[data[0]] = [data[1],data[2]]
	file.close()
	
	
	Phase(options.vcf, dicoPCT, options.prefix, options.region, options.outgp, options.win, float(options.maxIUP), int(options.MinLength))
		
if __name__ == "__main__": __main__()
