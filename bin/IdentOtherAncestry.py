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
import numpy as np

sys.stdout.write('modules loaded\n')


def CalcRatio(VCF, NAME, ANCESTOR, NEWANCESTOR, OUT):
	
	# Recording ancestral accessions
	dicoAnc = set()
	file = open(ANCESTOR)
	for line in file:
		data = line.split()
		if data:
			dicoAnc.add(data[0])
	file.close()
	
	# Recording potentially ancestral accessions
	ListeNewAnc = []
	file = open(NEWANCESTOR)
	for line in file:
		data = line.split()
		if data:
			if not(data[0] in ListeNewAnc):
				ListeNewAnc.append(data[0])
	ListeAccStat = [0]*len(ListeNewAnc)
	file.close()
	
	# Working on the vcf
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	outfile = gzip.open(OUT+'For_Density.tab.gz', 'wt')
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				list_to_print_in_output = []
				## Looking for accessions position
				AccPos = header.index(NAME)
				dicoPosAnc = {}
				for n in dicoAnc:
					dicoPosAnc[n] = header.index(n)
				dicoNewAnc = {}
				for i in range(len(ListeNewAnc)):
					dicoNewAnc[i] = header.index(ListeNewAnc[i])
				FILTERpos = header.index('FILTER')
				REFpos = header.index('REF')
				ALTpos = header.index('ALT')
				CHRpos = header.index('#CHROM')
				POSpos = header.index('POS')
				FLAGformat = header.index('FORMAT')
				TOTAL = 0
				TOTALOTHERORIGIN = 0
				TOTALALLELE = 0
			# Printing headers
			elif data[0][0] == '#':
				pass
			# Working on variant
			else:
				flag_format = data[FLAGformat].split(':')
				if 'AD' in flag_format and 'DP' in flag_format and 'GT' in flag_format:
					GTpos = flag_format.index('GT')
					DPpos = flag_format.index('DP')
					ADpos = flag_format.index('AD')
					# Recording accession genotype
					ACCGeno = set(data[AccPos].split(':')[GTpos].replace('|', '/').split('/'))
					ACCGeno.discard('.')
					TOTALALLELE += len(ACCGeno)
					# Recording ancestral alleles
					SetAnc = set()
					for n in dicoPosAnc:
						ANCGeno = set(data[dicoPosAnc[n]].split(':')[GTpos].replace('|', '/').split('/'))
						SetAnc.update(ANCGeno)
					remainingAlleles = list(ACCGeno-SetAnc)
					if remainingAlleles:
						ALLELE = [data[REFpos]]+data[ALTpos].split(',')
						for allele in remainingAlleles:
							outfile.write("\t".join([data[CHRpos], data[POSpos], ALLELE[int(allele)], "Unknown", "1"]))
							outfile.write("\n")
						NewOr = set()
						TOTAL += len(remainingAlleles)
						for acc in dicoNewAnc:
							NANC = set(data[dicoNewAnc[acc]].split(':')[GTpos].replace('|', '/').split('/'))
							for allele in remainingAlleles:
								if allele in NANC:
									ListeAccStat[acc] += 1
									NewOr.add(allele)
						TOTALOTHERORIGIN += len(NewOr)
	file.close()
	outfile.close()
	outfile = open(OUT, 'w')
	outfile.write('\t'.join(['Accession', 'TOTALALLELE', 'TOTAL', 'TOTALNewOrigin']+ListeNewAnc+ListeNewAnc))
	outfile.write('\n')
	outfile.write('\t'.join(map(str, [NAME, str(TOTALALLELE),str(TOTAL), str(TOTALOTHERORIGIN)] + ListeAccStat + list(np.array(ListeAccStat)/TOTAL))))
	outfile.write('\n')
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program look for the proportion of allele present on a studied accession"
	"\nthat are not present in set of accessions identified as ancestor. It also look"
	"\nfor the proportion of these alleles not present in defined ancestors in potential"
	"\nother (new) ancestors.")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file.')
	parser.add_option( '-n',	'--name',			dest='name',		default=None,			help='Accession name in which alleles will be compared')
	parser.add_option( '-a',	'--ancestor',		dest='ancestor',	default=None,			help='A tabulated file that contain on column 1 the name of accessions used to look for unspecific alleles')
	parser.add_option( '-l',	'--newancestor',	dest='newancestor',	default=None,			help='A tabulated file that contain on column 1 the name of accessions in which specific alleles from the searched accession will be compared')
	parser.add_option( '-o',	'--out',			dest='out',			default='Specific.tab', help='The output file. [Default: %default]')
	
	(options, args) = parser.parse_args()
	
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.name == None:
		sys.exit('Please provide a name file to --name argument')
	CalcRatio(options.vcf, options.name, options.ancestor, options.newancestor, options.out)
		
if __name__ == "__main__": __main__()
