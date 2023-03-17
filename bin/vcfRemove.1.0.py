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
import gzip
import optparse

def recup_geno(LINE, ACCESSION, FORMAT, REF, ALT, CHROM, POS):
	
	"""
		Return the genotype of an accession from a vcf line
		
		:param LINE: A list corresponding to a vcf line
		:type LINE: list
		:param ACCESSION: Acession position in LINE to return the genotype
		:type ACCESSION: int
		:param FORMAT: Position of the FORMAT column in LINE
		:type FORMAT: int
		:param REF: Position of the REF column in LINE
		:type REF: int
		:param ALT: Position of the ALT column in LINE
		:type ALT: int
		:param CHROM: Position of the CHROM column in LINE
		:type CHROM: int
		:param POS: Position of the POS column in LINE
		:type POS: int
		:return: A list with [0] --> chromosome, [1] --> position (integer), [2] --> allele 1, [3] --> allele 2, ... [n] --> allele [n-1]
		:rtype: list
	"""
	
	# Getting possible alleles
	liste_allele = []
	# Getting reference allele
	liste_allele.append(LINE[REF])
	# Getting alternative alleles
	alt_all = LINE[ALT].split(',')
	for n in alt_all:
		liste_allele.append(n)
	
	# Getting accession format
	format_col = LINE[FORMAT].split(':')
	
	# Getting genotypes
	geno = LINE[ACCESSION].split(':')[format_col.index('GT')].replace('/', '|').split('|')
	
	# Preparing the list to return
	list2return = [LINE[CHROM], LINE[POS],[], []]
	for n in geno:
		if n == '.':
			list2return[2].append('.')
		else:
			list2return[2].append(int(n))
	
	# Obtaining allele coverage
	AllDp = list(map(int, LINE[ACCESSION].split(':')[format_col.index('AD')].split(',')))
	list2return[3] = AllDp
	
	# returning list
	return list2return
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(
		description="This program remove alleles in a given accession if they are shared by a "
		"second accession. For Example, if an accession in which we will remove allele is A/A/T "
		"and a second one which is the one that will be used to remove is A/A, the result for "
		"the first accession will be T in the final vcf. If the accession used to remove is T/T, "
		"the result will be A/A as only T is common between the two accessions.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-a',	'--acc',			dest='acc',			default=None,			help='Accession name in which alleles will be removed. [Default: %default]')
	parser.add_option( '-r',	'--remove',			dest='remove',		default=None,			help='Accession used to remove alleles. [Default: %default]')
	parser.add_option( '-A',	'--onlyAcc',		dest='onlyAcc',		default='n',			help='Output only accession passed to --acc argument. Possible values: "y" or "n". [Default: %default]')
	parser.add_option( '-o',	'--out',			dest='out',			default=None,			help='Prefix for output files names. [Default: %default]')

	(options, args) = parser.parse_args()
	
	
	ACC = options.acc
	REMOVE = options.remove
	
	OK = 0
	NoOk = 0
	Missing = 0
	
	# recording informations in a dictionary
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
	else:
		file = open(options.vcf)
	outfile = gzip.open(options.out+'.vcf.gz','wt')
	outfile1 = gzip.open(options.out+'_haplo.vcf.gz','wt')
	outfile2 = gzip.open(options.out+'_MoreThanHaplo.vcf.gz','wt')
	for line in file:
		data = line.split()
		if data:
			if data[0] == '#CHROM':
				header = data
				CHRPOS = data.index('#CHROM')
				POSPOS = data.index('POS')
				IDPOS = data.index('ID')
				REFPOS = data.index('REF')
				ALTPOS = data.index('ALT')
				QUALPOS = data.index('QUAL')
				FILTERPOS = data.index('FILTER')
				INFOPOS = data.index('INFO')
				FORMATPOS = data.index('FORMAT')
				ACCPOS = data.index(ACC)
				REMOVEPOS = data.index(REMOVE)
				if options.onlyAcc == 'y':
					outfile.write('\t'.join([data[CHRPOS], data[POSPOS], data[IDPOS], data[REFPOS], data[ALTPOS], data[QUALPOS], data[FILTERPOS], data[INFOPOS], data[FORMATPOS], data[ACCPOS]])+'\n')
					outfile1.write('\t'.join([data[CHRPOS], data[POSPOS], data[IDPOS], data[REFPOS], data[ALTPOS], data[QUALPOS], data[FILTERPOS], data[INFOPOS], data[FORMATPOS], data[ACCPOS]])+'\n')
					outfile2.write('\t'.join([data[CHRPOS], data[POSPOS], data[IDPOS], data[REFPOS], data[ALTPOS], data[QUALPOS], data[FILTERPOS], data[INFOPOS], data[FORMATPOS], data[ACCPOS]])+'\n')
				else:
					outfile.write(line)
					outfile1.write(line)
					outfile2.write(line)
			elif data[0][0] != "#":
				geno_acc = recup_geno(data, ACCPOS, FORMATPOS, REFPOS, ALTPOS, CHRPOS, POSPOS)
				geno_remove = recup_geno(data, REMOVEPOS, FORMATPOS, REFPOS, ALTPOS, CHRPOS, POSPOS)
				ExpectedPlo = len(geno_acc[2])-len(geno_remove[2])
				if '.' in geno_acc or '.' in geno_remove:
					Missing += 1
				else:
					genotype = geno_acc[2][:]
					removed = []
					for allele in geno_remove[2]:
						if allele in genotype:
							genotype.remove(allele)
							removed.append(allele)
					list_to_print = []
					
					for n in set(geno_acc[2]):
						if n != '.':
							geno_acc[3][n] = round(geno_acc[3][n] * (1-(removed.count(n)/geno_acc[2].count(n))))
					TotalDP = sum(geno_acc[3])
					
					FORMAT = data[FORMATPOS].split(':')
					GENOACC = data[ACCPOS].split(':')
					GENOACC[FORMAT.index('GT')] = '/'.join(list(map(str, genotype)))
					GENOACC[FORMAT.index('AD')] = ','.join(list(map(str, geno_acc[3])))
					GENOACC[FORMAT.index('DP')] = str(TotalDP)
					
					if options.onlyAcc == 'y':
						list_to_print += [data[CHRPOS], data[POSPOS], data[IDPOS], data[REFPOS], data[ALTPOS], data[QUALPOS], data[FILTERPOS], data[INFOPOS], data[FORMATPOS]]
						list_to_print.append(':'.join(GENOACC))
					else:
						list_to_print = data[:]
						list_to_print[ACCPOS] = ':'.join(GENOACC)
					outfile.write('\t'.join(list_to_print+['\n']))
					if len(genotype) == ExpectedPlo:
						outfile1.write('\t'.join(list_to_print+['\n']))
						OK += 1
					else:
						outfile2.write('\t'.join(list_to_print+['\n']))
						NoOk += 1
			else:
				outfile.write(line)
				outfile1.write(line)
				outfile2.write(line)
	print('Total sites without missing:', OK+NoOk)
	print('OK sites without missing:', OK, str(OK/(OK+NoOk)))
	print('NoOK sites without missing:', NoOk, str(NoOk/(OK+NoOk)))
	print('Missing sites:', Missing, str(Missing/(OK+NoOk+Missing)))
	
	
					
if __name__ == "__main__": __main__()