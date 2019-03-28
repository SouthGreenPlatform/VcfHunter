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

def recup_geno(LINE, HEADER, ACCESSION):
	
	"""
		Return the genotype of an accession from a vcf line
		
		:param LINE: A list corresponding to a vcf line
		:type LINE: list
		:param HEADER: A list containing the header of the vcf
		:type HEADER: list
		:param ACCESSION: ACCESSION name to return the genotype
		:type ACCESSION: str
		:return: A list with [0] --> chromosome, [1] --> position (integer), [2] --> allele 1, [3] --> allele 2, ... [n] --> allele [n-1]
		:rtype: list
	"""
	
	# Getting possible alleles
	liste_allele = []
	# Getting reference allele
	liste_allele.append(LINE[HEADER.index('REF')])
	# Getting alternative alleles
	alt_all = LINE[HEADER.index('ALT')].split(',')
	for n in alt_all:
		liste_allele.append(n)
	
	# Getting accession format
	format_col = LINE[HEADER.index('FORMAT')].split(':')
	
	# Getting genotypes
	geno = LINE[HEADER.index(ACCESSION)].split(':')[format_col.index('GT')].replace('/', '|').split('|')
	
	# Preparing the list to return
	list2return = [LINE[HEADER.index('#CHROM')], LINE[HEADER.index('POS')]]
	for n in geno:
		if n == '.':
			list2return.append('NA')
		else:
			list2return.append(int(n))
	
	# returning list
	return list2return
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-a',	'--acc',			dest='acc',			default=None,			help='Accession name in which alleles will be removed. [Default: %default]')
	parser.add_option( '-r',	'--remove',			dest='remove',		default=None,			help='Accession used to remove alleles. [Default: %default]')
	parser.add_option( '-o',	'--out',			dest='out',			default=None,			help='Prefix for output files names. [Default: %default]')

	(options, args) = parser.parse_args()
	
	
	ACC = options.acc
	REMOVE = options.remove
	
	# recording informations in a dictionnary
	file = open(options.vcf)
	outfile = open(options.out+'.vcf','w')
	outfile1 = open(options.out+'_haplo.vcf','w')
	outfile2 = open(options.out+'_MoreThanHaplo.vcf','w')
	for line in file:
		data = line.split()
		if data:
			if data[0] == '#CHROM':
				header = data
				outfile.write(line)
				outfile1.write(line)
				outfile2.write(line)
			elif data[0][0] != "#":
				geno_acc = recup_geno(data, header, ACC)
				geno_remove = recup_geno(data, header, REMOVE)
				genotype = geno_acc[2:]
				for allele in geno_remove[2:]:
					if allele in genotype:
						genotype.remove(allele)
				list_to_print = []
				for n in range(len(header)):
					if header[n] == ACC:
						list_to_print.append(':'.join(['/'.join(list(map(str, genotype)))]+data[n].split(':')[1:]))
					else:
						list_to_print.append(data[n])
				outfile.write('\t'.join(list_to_print+['\n']))
				if len(genotype) == 1:
					outfile1.write('\t'.join(list_to_print+['\n']))
				else:
					outfile2.write('\t'.join(list_to_print+['\n']))
			else:
				outfile.write(line)
				outfile1.write(line)
				outfile2.write(line)
	
	
					
if __name__ == "__main__": __main__()