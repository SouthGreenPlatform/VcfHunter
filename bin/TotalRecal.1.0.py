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

def Recal(DATA, MINCOV, MAXCOV, MINALCOV, MINFREQ, ACCESSION_START, CHRpos, POSpos, REFpos, ALTpos, FORMATpos, DIAL, DICOREG, DICOACCPOS):
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
		param DICOREG: A dictionary containing accessions and given ploidy
		type DICOREG: dictionary
		param DICOACCPOS: A dictionary containing key: position in header and value: accession name
		type DICOACCPOS: dictionary
		
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
	
	dico_allele = set(range(len([REF]+ALT)))
	# dico_allele = set() #####
	dico_acc = set()
	
	# Cheking if Flags are present in the FORMAT tag
	if 'AD' in FORMAT and 'DP' in FORMAT and 'GT' in FORMAT:
	
		# Identification of accessions that pass MINCOV an MAXCOV (This is not optimized but it use same bloc as VcfPrefilter. Modified section of this bloc was identified with "#####")
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
									dico_acc.add(acc)
									# dico_allele.add(n) #####
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
			AltAllele = ','.join(allele2print)
			
			if AltAllele != DATA[ALTpos]:
				sys.exit('Oups, there is a bug1\n')
			
			for acc in range(len(DATA)):
				if acc >= ACCESSION_START:
					accName = DICOACCPOS[acc]
					if accName in DICOREG:
						code = []
						data = DATA[acc].split(':')
						ploidy = DICOREG[accName]
						allele_cov_info = list(map(int, data[ADPOS].split(',')))
						# determination of allele to print coverage
						coverage2print = []
						if not (0 in list_allele):
							coverage2print.append(allele_cov_info[0])
						for allele in list_allele:
							coverage2print.append(allele_cov_info[allele])
						
						if acc in dico_acc: # Testing if accession pass min and max coverage threshold
							genotype = utils.genotype_accession(coverage2print, list_allele_select, 0.001, str(ploidy), True, ALLCOMB)
							code = ':'.join([genotype[0], ','. join(list(map(str, coverage2print))), str(sum(coverage2print)), str(genotype[1])])
						else:
							# No enough coverage to genotype
							code = ':'.join(['/'.join(['.']*ploidy), ','.join(list(map(str, coverage2print))), str(sum(coverage2print)), '.'])
					else:
						code = DATA[acc]
					
					genotype_coding.append(code)
			list2print = [CHR, POS, '.', REF, AltAllele, '.', '.', '.', ':'.join(FORMAT+['GC'])] + genotype_coding
			return list2print
		else:
			return 0

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)\n"
	"This script recode genotype calling in specified regions of specified accessions according to given ploidy")
	# Wrapper options.
	parser.add_option( '-v', '--vcf',		dest='vcf',			default=None,		help='The VCF file')
	parser.add_option( '-t', '--table',		dest='table',		default=None,		help='A table file with column1: accession name, column2: chromosome name; column3: start region, column4: end region, column 5: ploidy')
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
	if options.table == None:
		sys.exit('Please provide a file to -t options')
	
	# Creating filtered output file
	if options.outgzip == 'n':
		outvcf = open(options.out,'w')
	elif options.outgzip == 'y':
		outvcf = gzip.open(options.out,'wt')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	# Getting options
	MINCOV = int(options.MinCov)
	MAXCOV = int(options.MaxCov)
	MINALCOV = int(options.MinAlCov)
	MINFREQ = float(options.minFreq)
	
	# Recording regions to work with
	dicoReg = {}
	dicoAcc = set()
	file = open(options.table)
	for line in file:
		data = line.split()
		if data:
			acc = data[0]
			chr = data[1]
			start = int(data[2])
			end = int(data[3])
			plo = int(data[4])
			dicoAcc.add(acc)
			if not(chr in dicoReg):
				dicoReg[chr] = {}
			for i in range((end-start) + 1):
				pos = i + start
				if not(pos in dicoReg[chr]):
					dicoReg[chr][pos] = {}
				dicoReg[chr][pos][acc] = plo
	file.close()
	
	# Working line by line
	PrintFilter = 1
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
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
			LineSize = len(header)
			Accession_start = header.index('FORMAT')+1
			CHRpos = header.index("#CHROM")
			POSpos = header.index("POS")
			REFpos = header.index("REF")
			ALTpos = header.index("ALT")
			FORMATpos = header.index("FORMAT")
			
			# Recording accessions position
			dicoAccPos = {}
			for acc in header:
				dicoAccPos[header.index(acc)] = acc
			
			# Just a piece of verification
			for acc in dicoAcc:
				if not(acc in header):
					sys.stdout.write("Warning, accession "+acc+" was not found in the vcf. This could be a problem according to what you want to do...\n")
			
		elif data[0][0] != "#":
			# Selecting variant line
			CHR = data[CHRpos]
			POS = int(data[POSpos])
			if CHR in dicoReg:
				if POS in dicoReg[CHR]:
					liste = Recal(data, MINCOV, MAXCOV, MINALCOV, MINFREQ, Accession_start, CHRpos, POSpos, REFpos, ALTpos, FORMATpos, options.dial, dicoReg[CHR][POS], dicoAccPos)
					if len(liste) != LineSize:
						sys.exit('Oups, there is a bug2\n')
					if liste:
						outvcf.write('\t'.join(liste))
						outvcf.write('\n')
				else:
					outvcf.write(line)
			else:
				outvcf.write(line)
		else:
			outvcf.write(line)
			
if __name__ == "__main__": __main__()
