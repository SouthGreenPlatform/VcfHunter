
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
import optparse
import sys
import time
import math
import gzip

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
	
	dico_proba = {}
	# calculating probabilities for homozygous state
	p = 1-ERROR
	for n in range(len(ALLELE)):
		dico_proba['/'.join([str(n)]*int(PLOIDY))] = binom(coverage[n],nb_tirage,p)/binom(int(round(nb_tirage*p)),nb_tirage,p)
	
		if PLOIDY == '2':
			# calculating probabilities for heterozygous state
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n),str(k)])] = (binom(coverage[n], nb_tirage, 0.5)*binom(coverage[k], nb_tirage, 0.5))/(binom(int(round(nb_tirage*0.5)),nb_tirage,0.5)*binom(int(round(nb_tirage*0.5)),nb_tirage,0.5))
							comb_done.append(couple)
	
		if PLOIDY == '3':
			# calculating probabilities for heterozygous state for triploids
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 2.0/3.0)*binom(coverage[n], nb_tirage, 1.0/3.0))/(binom(int(round(nb_tirage*(2.0/3.0))),nb_tirage,(2.0/3.0))*binom(int(round(nb_tirage*(1.0/3.0))),nb_tirage,(1.0/3.0)))
							comb_done.append(couple)
						couple = sorted([n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 2.0/3.0)*binom(coverage[k], nb_tirage, 1.0/3.0))/(binom(int(round(nb_tirage*(2.0/3.0))),nb_tirage,(2.0/3.0))*binom(int(round(nb_tirage*(1.0/3.0))),nb_tirage,(1.0/3.0)))
							comb_done.append(couple)
		
		if PLOIDY == '4':
			# calculating probabilities for heterozygous state for tetraploid
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 3.0/4.0)*binom(coverage[n], nb_tirage, 1.0/4.0))/(binom(int(round(nb_tirage*(3.0/4.0))),nb_tirage,(3.0/4.0))*binom(int(round(nb_tirage*(1.0/4.0))),nb_tirage,(1.0/4.0)))
							comb_done.append(couple)
						couple = sorted([n,n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 2.0/4.0)*binom(coverage[k], nb_tirage, 2.0/4.0))/(binom(int(round(nb_tirage*(2.0/4.0))),nb_tirage,(2.0/4.0))*binom(int(round(nb_tirage*(2.0/4.0))),nb_tirage,(2.0/4.0)))
							comb_done.append(couple)
						couple = sorted([n,n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 3.0/4.0)*binom(coverage[k], nb_tirage, 1.0/4.0))/(binom(int(round(nb_tirage*(3.0/4.0))),nb_tirage,(3.0/4.0))*binom(int(round(nb_tirage*(1.0/4.0))),nb_tirage,(1.0/4.0)))
							comb_done.append(couple)
		
		if PLOIDY == '12':
			# calculating probabilities for heterozygous state for tetraploid
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 11.0/12.0)*binom(coverage[n], nb_tirage, 1.0/12.0))/(binom(int(round(nb_tirage*(11.0/12.0))),nb_tirage,(11.0/12.0))*binom(int(round(nb_tirage*(1.0/12.0))),nb_tirage,(1.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,k,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 10.0/12.0)*binom(coverage[n], nb_tirage, 2.0/12.0))/(binom(int(round(nb_tirage*(10.0/12.0))),nb_tirage,(10.0/12.0))*binom(int(round(nb_tirage*(2.0/12.0))),nb_tirage,(2.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 9.0/12.0)*binom(coverage[n], nb_tirage, 3.0/12.0))/(binom(int(round(nb_tirage*(9.0/12.0))),nb_tirage,(9.0/12.0))*binom(int(round(nb_tirage*(3.0/12.0))),nb_tirage,(3.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 8.0/12.0)*binom(coverage[n], nb_tirage, 4.0/12.0))/(binom(int(round(nb_tirage*(8.0/12.0))),nb_tirage,(8.0/12.0))*binom(int(round(nb_tirage*(4.0/12.0))),nb_tirage,(4.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 7.0/12.0)*binom(coverage[n], nb_tirage, 5.0/12.0))/(binom(int(round(nb_tirage*(7.0/12.0))),nb_tirage,(7.0/12.0))*binom(int(round(nb_tirage*(5.0/12.0))),nb_tirage,(5.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 6.0/12.0)*binom(coverage[n], nb_tirage, 6.0/12.0))/(binom(int(round(nb_tirage*(6.0/12.0))),nb_tirage,(6.0/12.0))*binom(int(round(nb_tirage*(6.0/12.0))),nb_tirage,(6.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 7.0/12.0)*binom(coverage[k], nb_tirage, 5.0/12.0))/(binom(int(round(nb_tirage*(7.0/12.0))),nb_tirage,(7.0/12.0))*binom(int(round(nb_tirage*(5.0/12.0))),nb_tirage,(5.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 8.0/12.0)*binom(coverage[k], nb_tirage, 4.0/12.0))/(binom(int(round(nb_tirage*(8.0/12.0))),nb_tirage,(8.0/12.0))*binom(int(round(nb_tirage*(4.0/12.0))),nb_tirage,(4.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 9.0/12.0)*binom(coverage[k], nb_tirage, 3.0/12.0))/(binom(int(round(nb_tirage*(9.0/12.0))),nb_tirage,(9.0/12.0))*binom(int(round(nb_tirage*(3.0/12.0))),nb_tirage,(3.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 10.0/12.0)*binom(coverage[k], nb_tirage, 2.0/12.0))/(binom(int(round(nb_tirage*(10.0/12.0))),nb_tirage,(10.0/12.0))*binom(int(round(nb_tirage*(2.0/12.0))),nb_tirage,(2.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 11.0/12.0)*binom(coverage[k], nb_tirage, 1.0/12.0))/(binom(int(round(nb_tirage*(11.0/12.0))),nb_tirage,(11.0/12.0))*binom(int(round(nb_tirage*(1.0/12.0))),nb_tirage,(1.0/12.0)))
							comb_done.append(couple)
	
	# getting best probability
	best_genotype = '/'.join(['.']*int(PLOIDY))
	best_value = 0
	second_best = 0
	for genotype in dico_proba:
		# print (dico_proba[genotype])
		if best_value < dico_proba[genotype]:
			if best_value != 0:
				second_best = best_value
			best_value = dico_proba[genotype]
			best_genotype = genotype
		elif best_value == dico_proba[genotype]:
			best_genotype = '/'.join(['.']*int(PLOIDY))
		elif second_best < dico_proba[genotype]:
			second_best = dico_proba[genotype]
		# print (best_value, second_best)
	# print (best_value, second_best, dico_proba)
	# print (-math.log(second_best,10)+math.log(best_value,10))
	if best_value == 0:
		return best_genotype, (1)
	elif second_best == 0:
		return best_genotype, (9999)
	else:
		return best_genotype, best_value/second_best

def filter_on_read_cov(DATA, MINCOV, MAXCOV, MINALCOV, MINFREQ, VCF_HEADER, ACCESSION_START):
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
		type MINFREQ: str
		param VCF_HEADER: VCF header
		type VCF_HEADER: list
		param ACCESSION_START: Position of accession start
		type ACCESSION_START: int
		return: list
	"""
	
	minfreq = float(MINFREQ)
	
	line_status = 0
	# 0 : OK
	# 1 : bad FORMAT tag
	
	# Getting line information
	CHR = DATA[VCF_HEADER.index("#CHROM")]
	POS = DATA[VCF_HEADER.index("POS")]
	REF = DATA[VCF_HEADER.index("REF")]
	ALT = DATA[VCF_HEADER.index("ALT")].split(",")
	FORMAT = DATA[VCF_HEADER.index("FORMAT")].split(":")
	
	dico_allele = set()
	dico_acc = set()
	
	# Cheking if Flags are present in the FORMAT tag
	if 'AD' in FORMAT and 'DP' in FORMAT and 'GT' in FORMAT:
	
		# Identification if the line pass coverages and minimum allele frequency ratio (MINFREQ and min_all_cov)
		for acc in range(len(DATA)):
			if acc >= ACCESSION_START:
				data_point_status = 0
				# 0 : OK
				# 1 : bad format
				# 2 : no coverage
				data = DATA[acc].split(':')
				# Cheking if the accession is in the same format of the FORMAT tag
				if len(FORMAT) == len(data):
					# Cheking if there is information in the DP and AD tag
					if '.' in data[FORMAT.index('DP')] or '.' in data[FORMAT.index('AD')]:
						data_point_status = 2
					elif int(data[FORMAT.index('DP')]) < MINCOV:
						data_point_status = 2
					elif int(data[FORMAT.index('DP')]) > MAXCOV:
						data_point_status = 2
					else:
						# additional step to identify between false and true multiallelic variant sites
						allele_cov_info = list(map(int, data[FORMAT.index('AD')].split(',')))
						# Checking for cases where the coverage is due to non A,T,G,C bases
						if sum(allele_cov_info) > 0:
							for n in range(len(allele_cov_info)):
								if allele_cov_info[n]/float(sum(allele_cov_info)) >= minfreq and allele_cov_info[n] >= MINALCOV:
									dico_allele.add(n)
									dico_acc.add(acc)
				else:
					data_point_status = 1
		# Now we are going to work on line with at least 1 variant site different from the reference
		if len(dico_allele) > 1 or (not(0 in dico_allele) and len(dico_allele) == 1):
			# print (dico_allele, allele_cov_info, DATA)
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
					ploidy = len(data[FORMAT.index('GT')].split('/'))
					allele_cov_info = list(map(int, data[FORMAT.index('AD')].split(',')))
					# determination of allele to print coverage
					coverage2print = []
					if not (0 in list_allele):
						coverage2print.append(allele_cov_info[0])
					for allele in list_allele:
						coverage2print.append(allele_cov_info[allele])
					
					if acc in dico_acc:
						genotype = genotype_accession(coverage2print, list_allele_select, 0.001, str(ploidy))
						code = ':'.join([genotype[0], ','. join(list(map(str, coverage2print))), str(sum(coverage2print)), str(genotype[1])])
					else:
						# No enough coverage to genotype
						code = ':'.join(['/'.join(['.']*ploidy), ','.join(list(map(str, coverage2print))), str(sum(coverage2print)), '.'])
					
					genotype_coding.append(code)
			list2print = [CHR, POS, '.', REF, ','.join(allele2print), '.', '.', '.', ':'.join(FORMAT+['GC'])] + genotype_coding
			# print (list2print)
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
	parser.add_option( '-o', '--out',		dest='out',			default='Pop',		help='Prefix for output files. [Default: %default]')
	parser.add_option( '-g', '--outgzip',	dest='outgzip',		default='n',		help='Output files in gzip format. [Default: %default]')
	(options, args) = parser.parse_args()
	
	if options.vcf == None:
		sys.exit('Please provide a vcf file to -v options')
	
	# Creating filtered output file
	if options.outgzip == 'n':
		outvcf = open(options.out,'w')
	elif options.outgzip == 'y':
		outvcf = gzip.open(options.out+'.gz','w')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	
	# GLOBAL VARIABLES
	GLOB_allele_ratio = []
	GLOB_coverage = []
	
	
	TO_EXCLUDE = set()
	NO_READ, NO_COV, TOO_COV, NO_FREQ, NO_COV_PLUS, TRONCATED_FORMAT = 0, 0, 0, 0, 0, 0
	FLTR_MISS = 0
	FLTR_PVAL = 0
	PASSED = 0
	TO_MUCH_ALLELE = 0
	
	MORE_THAN_DI_ALLELE, ALL_MISSING, MONOMORPHOUS, DIMORPHOUS_SITES, MORE_THAN_DIMORPHOUS = 0, 0, 0, 0, 0
	
	
	
	# Working line by line
	PrintFilter = 1
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
			Accession_header = header[Accession_start:]
			
		elif data[0][0] != "#":
			alleles = [data[header.index('REF')]]
			alleles = alleles + data[header.index('ALT')].split(',')
			chrom = data[header.index('#CHROM')]
			pos = data[header.index('POS')]
			format = data[header.index('FORMAT')].split(':')
			# Filtering on coverage and allele frequency
			liste = filter_on_read_cov(data, int(options.MinCov), int(options.MaxCov), int(options.MinAlCov), options.minFreq, header, Accession_start)
			if liste:
				outvcf.write('\t'.join(liste))
				outvcf.write('\n')
		else:
			outvcf.write(line)
			
if __name__ == "__main__": __main__()
