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


def Format(VCF, PREFIX, OUTGZIP):
	
	"""
		Filter vcf file and output a filtered vcf
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
		:param OUTGZIP: value to tell to gzip output
		:type RMTYPE: str
	"""
	
	# Creating output file
	if OUTGZIP == 'n':
		outfile = open(PREFIX,'w')
	elif OUTGZIP == 'y':
		outfile = gzip.open(PREFIX+'.gz','wt')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
	RmBecauseHomo = 0
	
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				FILTERpos = header.index('FILTER')
				REFpos = header.index('REF')
				ALTpos = header.index('ALT')
				FLAGformat = header.index('FORMAT')
				Accession_start = header.index('FORMAT')+1
				dicoAccPos = {}
				for n in header[Accession_start:]:
					dicoAccPos[n] = header.index(n)
				outfile.write(line)
			# Printing headers
			elif data[0][0] == '#':
				outfile.write(line)
			# Working on variant
			else:
				# Identification of alleles present in accessions
				ref_allele = data[REFpos]
				alt_allele = data[ALTpos].split(',')
				# print(data[0], data[1], ref_allele, alt_allele)
				AlleleList = [ref_allele]+alt_allele
				flag_format = data[FLAGformat].split(':')
				GTpos = flag_format.index('GT')
				DPpos = flag_format.index('DP')
				ADpos = flag_format.index('AD')
				if 'AD' in flag_format and 'DP' in flag_format and 'GT' in flag_format:
					DicoAllele = set()
					for accession in data[Accession_start:]:
						GTAcc = set(accession.split(':')[GTpos].replace('|','/').split('/'))
						DicoAllele.update(GTAcc)
					if '.' in DicoAllele:
						DicoAllele.remove('.')
					ListAllele = sorted(list(DicoAllele))
					
					if len(ListAllele) > 1:
						dicAllele = {}
						for i in range(len(ListAllele)):
							dicAllele[ListAllele[i]] = i
						ref_allele = int(ListAllele[0])
						alt_allele = ListAllele[1:]
						
						data[REFpos] = AlleleList[ref_allele]
						data[ALTpos] = ','.join([AlleleList[int(x)] for x in alt_allele])
						
						for acc in header[Accession_start:]:
							ACCINFO = data[dicoAccPos[acc]].split(':')
							GTAcc = ACCINFO[GTpos].replace('|','/').split('/')
							if not('.' in GTAcc):
								FGTAcc = '/'.join([str(dicAllele[x]) for x in GTAcc])
							else:
								FGTAcc = ACCINFO[GTpos]
							
							ADcov = []
							ADacc = ACCINFO[ADpos].split(',')
							for val in ListAllele:
								ADcov.append(ADacc[int(val)])
							
							ACCINFO[GTpos] = FGTAcc
							ACCINFO[ADpos] = ','.join(ADcov)
							
							data[dicoAccPos[acc]] = ':'.join(ACCINFO)
						outfile.write('\t'.join(data))
						outfile.write('\n')
					else:
						RmBecauseHomo += 1
	file.close()
	outfile.close()
	
	sys.stdout.write('Removed sites because homozygous: '+str(RmBecauseHomo)+'\n')


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',		default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-o',	'--out',			dest='out',		default='WorkOnVcf', 	help='Name for output files. [Default: %default]')
	parser.add_option( '-g',	'--outgzip',		dest='outgzip',	default='n',			help='Output files in gzip format. [Default: %default]')
	
	(options, args) = parser.parse_args()
	
	# Filtering vcf file
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	Format(options.vcf, options.out, options.outgzip)
		
if __name__ == "__main__": __main__()
