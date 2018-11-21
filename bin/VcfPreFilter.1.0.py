#!/usr/bin/env python
#
#  Copyright 2017 CIRAD
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
import time
import math
import gzip
import optparse
import itertools
import utilsSR.utilsSR as utils

def combin(n, k):
	"""Nombre de combinaisons de n objets pris k a k"""
	if k > n//2:
		k = n-k
	x = 1
	y = 1
	i = n-k+1
	while i <= n:
		x = (x*i)//y
		y += 1
		i += 1
	return x

def binom(k,n,p):
	x = combin(n,k)*(p**k)*((1-p)**(n-k))
	return x

def genotype_accession(COVERAGE, ALLELE, ERROR, PLOIDY):
	
	# WARNING The statistics will not be accurate on low coverage
	
	# To manage huge factorial
	nb_tirage = sum(COVERAGE)
	if nb_tirage >= 1000:
		coverage = list(map(int, [x/float(nb_tirage)*100 for x in COVERAGE]))
	else:
		coverage = list(COVERAGE)
	nb_tirage = sum(coverage)
	
	# Obtaining combinations
	ListCombPair = list(itertools.combinations_with_replacement(range(len(ALLELE)), 2))
	FinalComb = set()
	for n in ListCombPair:
		N = set(n)
		if len(N) == 2:
			IntermList = list(itertools.combinations_with_replacement(n, int(PLOIDY)))
			for k in IntermList:
				FinalComb.add(k)
	
	dico_proba = {}
	# calculating probabilities for homozygous state
	p = 1-ERROR
	for n in FinalComb:
		N = set(n)
		Key = '/'.join(list(map(str, n)))
		if len(N) == 1:
			Value = N.pop()
			Proba = binom(coverage[Value],nb_tirage,p)/binom(int(round(nb_tirage*p)),nb_tirage,p)
		else:
			Proba = 1
			for k in N:
				Ratio = n.count(k)/float(PLOIDY)
				Proba = Proba*(binom(coverage[k], nb_tirage, Ratio)/binom(int(round(nb_tirage*Ratio)),nb_tirage,Ratio))
		dico_proba[Key] = Proba
	
	# getting best probability
	best_genotype = '/'.join(['.']*int(PLOIDY))
	best_value = 0
	second_best = 0
	for genotype in dico_proba:
		if best_value < dico_proba[genotype]:
			if best_value != 0:
				second_best = best_value
			best_value = dico_proba[genotype]
			best_genotype = genotype
		elif best_value == dico_proba[genotype]:
			best_genotype = '/'.join(['.']*int(PLOIDY))
		elif second_best < dico_proba[genotype]:
			second_best = dico_proba[genotype]
	if best_value == 0:
		return best_genotype, (1)
	elif second_best == 0:
		return best_genotype, (9999)
	else:
		return best_genotype, best_value/second_best

def filter_on_read_cov(DATA, MINCOV, MAXCOV, MINALCOV, MINFREQ, ACCESSION_START, CHRpos, POSpos, REFpos, ALTpos, FORMATpos, DIAL):
	"""
		Convert VCF calling to a list of recoded variant

		param DATA: list corresponding to a vcf line
		type DATA: list
		param MINCOV: minimal site coverage per accessions
		type MINCOV: int
		param MAXCOV: maximal site coverage per accessions
		type MAXCOV: int
		param MINALCOV: Minimum allele coverage threshold
		type MINALCOV: int
		param MINFREQ: minimal allele frequency
		type MINFREQ: float
		param ACCESSION_START: Position of accession start
		type ACCESSION_START: int
		param CHRpos: Position of chromosome name in DATA
		type CHRpos: int
		param POSpos: Position of chromosome position in DATA
		type POSpos: int
		param REFpos: Position of reference allele in DATA
		type REFpos: int
		param ALTpos: Position of alternate alleles in DATA
		type ALTpos: int
		param FORMATpos: Position of FORMAT tags in DATA
		type FORMATpos: int
		param DIAL: Perform only bi-allelic calling
		type DIAL: str
		
		return: list
	"""
	
	if DIAL == 'y':
		ALLCOMB = False
	elif DIAL == 'n':
		ALLCOMB = True
	else:
		sys.exit('Unmanaged argument '+DIAL+'in --dial option. Possible values are "y" or "n"\n')
	
	# Getting line information
	CHR = DATA[CHRpos]
	POS = DATA[POSpos]
	REF = DATA[REFpos]
	ALT = DATA[ALTpos].split(",")
	FORMAT = DATA[FORMATpos].split(":")
	DPPOS = FORMAT.index('DP')
	ADPOS = FORMAT.index('AD')
	GTPOS = FORMAT.index('GT')
	
	dico_allele = set()
	dico_acc = set()
	
	# Cheking if Flags are present in the FORMAT tag
	if 'AD' in FORMAT and 'DP' in FORMAT and 'GT' in FORMAT:
	
		# Identification if the line pass coverages and minimum allele frequency ratio (MINFREQ and min_all_cov)
		for acc in range(len(DATA)):
			if acc >= ACCESSION_START:
				data = DATA[acc].split(':')
				# Cheking if the accession is in the same format of the FORMAT tag
				if len(FORMAT) == len(data):
					# Cheking if there is information in the DP and AD tag
					if '.' in data[DPPOS] or '.' in data[ADPOS]:
						pass
					elif int(data[DPPOS]) < MINCOV:
						pass
					elif int(data[DPPOS]) > MAXCOV:
						pass
					else:
						# additional step to identify between false and true multiallelic variant sites
						allele_cov_info = list(map(int, data[ADPOS].split(',')))
						# Checking for cases where the coverage is due to non A,T,G,C,* bases
						if sum(allele_cov_info) > 0:
							for n in range(len(allele_cov_info)):
								if allele_cov_info[n]/float(sum(allele_cov_info)) >= MINFREQ and allele_cov_info[n] >= MINALCOV:
									dico_allele.add(n)
									dico_acc.add(acc)
				else:
					pass
		# Now we are going to work on line with at least 1 variant site different from the reference
		if len(dico_allele) > 1 or (not(0 in dico_allele) and len(dico_allele) == 1):
			list_allele = list(dico_allele)
			list_allele.sort()
			# determination of allele to work with and print as alternative
			list_allele_select = [REF]
			allele2print = []
			for allele in list_allele:
				if allele != 0:
					allele2print.append(ALT[allele-1])
					list_allele_select.append(ALT[allele-1])
			genotype_coding = []
			for acc in range(len(DATA)):
				if acc >= ACCESSION_START:
					code = []
					data = DATA[acc].split(':')
					ploidy = len(data[GTPOS].split('/'))
					allele_cov_info = list(map(int, data[ADPOS].split(',')))
					# determination of allele to print coverage
					coverage2print = []
					if not (0 in list_allele):
						coverage2print.append(allele_cov_info[0])
					for allele in list_allele:
						coverage2print.append(allele_cov_info[allele])
					
					if acc in dico_acc:
						genotype = utils.genotype_accession(coverage2print, list_allele_select, 0.001, str(ploidy), True, ALLCOMB)
						code = ':'.join([genotype[0], ','. join(list(map(str, coverage2print))), str(sum(coverage2print)), str(genotype[1])])
					else:
						# No enough coverage to genotype
						code = ':'.join(['/'.join(['.']*ploidy), ','.join(list(map(str, coverage2print))), str(sum(coverage2print)), '.'])
					
					genotype_coding.append(code)
			list2print = [CHR, POS, '.', REF, ','.join(allele2print), '.', '.', '.', ':'.join(FORMAT+['GC'])] + genotype_coding
			return list2print
		else:
			return 0

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)\n"
	"This script filter VCF file generated by Process_RNAseq.py or process_reseq.py by removing homozygous sites for all accessions based on parameters passed")
	# Wrapper options.
	parser.add_option( '-v', '--vcf',		dest='vcf',			default=None,		help='The VCF file')
	parser.add_option( '-m', '--MinCov',	dest='MinCov',		default='10',		help='Minimal read coverage for site. [Default: %default]')
	parser.add_option( '-M', '--MaxCov',	dest='MaxCov',		default='1000',		help='Maximal read coverage for site. [Default: %default]')
	parser.add_option( '-f', '--minFreq',	dest='minFreq',		default='0.05',		help='Minimal allele frequency in an accession to keep the allele for calling in the row. [Default: %default]')
	parser.add_option( '-c', '--MinAlCov',	dest='MinAlCov',	default='3',		help='Minimal read number of minor allele to call variant heterozygous (between 1 and infinity). [Default: %default]')
	parser.add_option( '-d', '--dial',		dest='dial',		default='y',		help='Perform only a diallelic calling. i.e Only two allele are possible in a genotype if "y" is passed to this argument. Possible values "y" or "n". [Default: %default]')
	parser.add_option( '-o', '--out',		dest='out',			default='Pop',		help='Prefix for output files. [Default: %default]')
	parser.add_option( '-g', '--outgzip',	dest='outgzip',		default='n',		help='Output files in gzip format. [Default: %default]')
	(options, args) = parser.parse_args()
	
	if options.vcf == None:
		sys.exit('Please provide a vcf file to -v options')
	
	# Creating filtered output file
	if options.outgzip == 'n':
		outvcf = open(options.out,'w')
	elif options.outgzip == 'y':
		outvcf = gzip.open(options.out,'wb')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	# Working line by line
	PrintFilter = 1
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rb')
	else:
		file = open(options.vcf)
	for line in file:
		data = line.split()
		if data[0][0:8] == "##contig" and PrintFilter:
			outvcf.write('##FORMAT=<ID=GC,Number=1,Type=Float,Description="Ratio between best genotype probability and second best genotype probability">\n')
			PrintFilter = 0
		if data[0] == "#CHROM":
			outvcf.write(line)
			header = data
			Accession_start = header.index('FORMAT')+1
			# Accession_header = header[Accession_start:]
			CHRpos = header.index("#CHROM")
			POSpos = header.index("POS")
			REFpos = header.index("REF")
			ALTpos = header.index("ALT")
			FORMATpos = header.index("FORMAT")
			
		elif data[0][0] != "#":
			# Filtering on coverage and allele frequency
			liste = filter_on_read_cov(data, int(options.MinCov), int(options.MaxCov), int(options.MinAlCov), float(options.minFreq), Accession_start, CHRpos, POSpos, REFpos, ALTpos, FORMATpos, options.dial)
			if liste:
				outvcf.write('\t'.join(liste))
				outvcf.write('\n')
		else:
			outvcf.write(line)
			
if __name__ == "__main__": __main__()
