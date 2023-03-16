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

sys.stdout.write('modules loaded\n')


def Phase(VCF, PCT, PREFIX):
	
	"""
		Filter vcf file and output a filtered vcf
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param PCT: A dictionary containing with keys = F1 and values = [P1,P2]
		:type PCT: str
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
	"""
	
	# Fixing order in which accessions are treated
	F1ORDER = sorted(list(PCT.keys()))
	
	# Creating dictionary to count problems (impossible case, heterozygous, multiallelic)
	dicoImpossible = {}
	for n in F1ORDER:
		dicoImpossible[n] = [0,0,0,0]
	
	# Creating output file
	outfile = gzip.open(PREFIX+'_Phase.vcf.gz','wt')
	
	# Reading and phasing line by line
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
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
				FILTERpos = header.index('FILTER')
				REFpos = header.index('REF')
				ALTpos = header.index('ALT')
				FLAGformat = header.index('FORMAT')
				
				for F1 in F1ORDER:
					header.append(''.join([F1,'-H1-','_X_'.join(PCT[F1])]))
					header.append(''.join([F1,'-H2-','_X_'.join(PCT[F1])]))
					header.append
					header.append(''.join([PCT[F1][0],'-H1-from-',F1]))
					header.append(''.join([PCT[F1][0],'-H2-from-',F1]))
					header.append
					header.append(''.join([PCT[F1][1],'-H1-from-',F1]))
					header.append(''.join([PCT[F1][1],'-H2-from-',F1]))
				outfile.write('\t'.join(header))
				outfile.write('\n')
			
			# Printing headers
			elif data[0][0] == '#':
				outfile.write(line)
			# Working on variant
			else:
				FinalLine = data
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
						
						OK = 0 # By default we hypothesize that it is an impossible marker
						
						# Testing if missing data
						AllInPtc = F1GT | P1GT | P2GT
						if '.' in AllInPtc:
							dicoImpossible[accession][3] += 1
						else:
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
						
						if OK:
							genotypeF1_H1 = ['.']*len(flag_format)
							genotypeF1_H2 = ['.']*len(flag_format)
							genotypeP1_H1 = ['.']*len(flag_format)
							genotypeP1_H2 = ['.']*len(flag_format)
							genotypeP2_H1 = ['.']*len(flag_format)
							genotypeP2_H2 = ['.']*len(flag_format)
							
							ADF1_H1 = ['0']*len(AlleleList)
							ADF1_H2 = ['0']*len(AlleleList)
							ADP1_H1 = ['0']*len(AlleleList)
							ADP1_H2 = ['0']*len(AlleleList)
							ADP2_H1 = ['0']*len(AlleleList)
							ADP2_H2 = ['0']*len(AlleleList)
							
							
							ADF1_H1[F1L] = list(F1info[ADpos].split(','))[F1L]
							ADF1_H2[F1R] = list(F1info[ADpos].split(','))[F1R]
							ADP1_H1[P1L] = list(P1info[ADpos].split(','))[P1L]
							ADP1_H2[P1R] = list(P1info[ADpos].split(','))[P1R]
							ADP2_H1[P2L] = list(P2info[ADpos].split(','))[P2L]
							ADP2_H2[P2R] = list(P2info[ADpos].split(','))[P2R]
							
							genotypeF1_H1[DPpos] = ADF1_H1[F1L]
							genotypeF1_H2[DPpos] = ADF1_H2[F1R]
							genotypeP1_H1[DPpos] = ADP1_H1[P1L]
							genotypeP1_H2[DPpos] = ADP1_H2[P1R]
							genotypeP2_H1[DPpos] = ADP2_H1[P2L]
							genotypeP2_H2[DPpos] = ADP2_H2[P2R]
							
							genotypeF1_H1[ADpos] = ','.join(ADF1_H1)
							genotypeF1_H2[ADpos] = ','.join(ADF1_H2)
							genotypeP1_H1[ADpos] = ','.join(ADP1_H1)
							genotypeP1_H2[ADpos] = ','.join(ADP1_H2)
							genotypeP2_H1[ADpos] = ','.join(ADP2_H1)
							genotypeP2_H2[ADpos] = ','.join(ADP2_H2)
							
							genotypeF1_H1[GTpos] = str(F1L)
							genotypeF1_H2[GTpos] = str(F1R)
							genotypeP1_H1[GTpos] = str(P1L)
							genotypeP1_H2[GTpos] = str(P1R)
							genotypeP2_H1[GTpos] = str(P2L)
							genotypeP2_H2[GTpos] = str(P2R)
							
							FinalLine.append(':'.join(genotypeF1_H1))
							FinalLine.append(':'.join(genotypeF1_H2))
							FinalLine.append(':'.join(genotypeP1_H1))
							FinalLine.append(':'.join(genotypeP1_H2))
							FinalLine.append(':'.join(genotypeP2_H1))
							FinalLine.append(':'.join(genotypeP2_H2))
							
						else:
							genotype = ['.']*len(flag_format)
							genotype[DPpos] = '0'
							genotype[ADpos] = ','.join(['0']*len(AlleleList))
							FinalLine.append(':'.join(genotype))
							FinalLine.append(':'.join(genotype))
							FinalLine.append(':'.join(genotype))
							FinalLine.append(':'.join(genotype))
							FinalLine.append(':'.join(genotype))
							FinalLine.append(':'.join(genotype))
					
					outfile.write('\t'.join(FinalLine))
					outfile.write('\n')


	for n in dicoImpossible:
		print('Impossible','HetOnAll','NotBiallelic', 'Missing Data')
		print(n, dicoImpossible[n])

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '',	'--names',			dest='names',		default=None,			help='A 3 column file containing in this order F1 P1 P2. [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')
	
	(options, args) = parser.parse_args()
	
	# Filtering vcf file
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.names == None:
		sys.exit('Please provide a name file to --names argument')
	
	# Recording information
	file = open(options.names)
	dicoPCT = {}
	for line in file:
		data = line.split()
		if data:
			dicoPCT[data[0]] = [data[1],data[2]]
	file.close()
	
	
	Phase(options.vcf, dicoPCT, options.prefix)
		
if __name__ == "__main__": __main__()
