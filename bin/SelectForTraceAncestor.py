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
import gzip



def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program identify sepecifi allele and format data so that they can be"
	"\nloaded to perform chromosome painting with TraceAncestor.")
	# Wrapper options. 
	parser.add_option( '',	'--conf',			dest='conf',		default=None,			help='Conf file containing vcf location (one per chromosome). [Default: %default]')
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='A 2 column file containing accession name (col1) origin/group (Col2). [Default: %default]')
	parser.add_option( '',	'--NoMiss',			dest='NoMiss',		default='n',			help='No missing data are allowed in accessions used to group alleles. Possible value: "y" (no missing allowed), "n" (missing allowed, but a group with only missing will remove the SNP). [Default: %default]')
	parser.add_option( '',	'--all',			dest='all',			default='n',			help='Allele should be present in all accessions of the group. Possible value: "y", "n". [Default: %default]')
	parser.add_option( '',	'--hom',			dest='hom',			default='n',			help='Allele should be present at homozygous state in all accessions of the group containing it (if not, it is removed). Possible value: "y", "n". [Default: %default]')
	parser.add_option( '',	'--out',			dest='out',			default='GpAllele.tab',	help='Output file name. [Default: %default]')
	(options, args) = parser.parse_args()
	
	if options.conf == None and options.vcf == None:
		sys.exit('Please provide either a conf file to --conf argument or a vcf file to --vcf argument')
	if options.conf != None and options.vcf != None:
		sys.exit('Both --conf and --vcf argument were provided. They are mutually exclusive... Please provide either a conf file to --conf argument or a vcf file to --vcf argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	
	# recording accession group
	dico_group = {}
	file = open(options.origin)
	for line in file:
		data = line.split()
		if data:
			if not(data[1] in dico_group):
				dico_group[data[1]] = []
			dico_group[data[1]].append(data[0])
	file.close()
	
	# recording vcf file to work with
	dico_vcf = set()
	if options.conf != None:
		file = open(options.conf)
		for line in file:
			data = line.split()
			if data:
				dico_vcf.add(data[0])
		file.close()
	else:
		dico_vcf.add(options.vcf)
	
	outfile = open(options.out, 'w')
	outfile.write('\t'.join(('ancestor', 'chromosome', 'position', 'allele')))
	outfile.write('\n')
	
	# For linear drawing
	dico_chr = {}
	total = []
	for vcf in dico_vcf:
		if vcf[-3:] == '.gz':
			file = gzip.open(vcf,'rt')
		else:
			file = open(vcf)
		for line in file:
			data = line.split()
			if data:
				if data[0] == "#CHROM":
					header = data
				elif data[0][0] != "#":
					FORMAT = data[header.index("FORMAT")].split(':')
					CHR = data[header.index("#CHROM")]
					POS = data[header.index("POS")]
					REF = [data[header.index("REF")]]
					ALT = data[header.index("ALT")].split(',')
					VARIANT = REF + ALT
					
					# recording group1 and group2 alleles
					dico_allele = {}
					for gp in dico_group:
						dico_allele[gp] = set()
						for ACC in dico_group[gp]:
							ACCESSION = data[header.index(ACC)].split(':')
							GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
							for geno in GENOTYPE:
								dico_allele[gp].add(geno)
					
					# to manage missing data in groups
					Good = 1
					for gp in dico_allele:
						if options.NoMiss == 'n':
							if '.' in dico_allele[gp] and len(dico_allele[gp]) == 1:
								Good = 0
						elif options.NoMiss == 'y':
							if '.' in dico_allele[gp]:
								Good = 0
						else:
							sys.exit('Oups, their is a bug... either the vcf has a problem or I made a mistake in the programing\n')
					if Good:
						# Selecting potential group specific alleles
						DicoSpecAll = {}
						for n in dico_allele:
							for all in dico_allele[n]:
								if all != '.':
									IsSpec = True
									for k in dico_allele:
										if n != k:
											if all in dico_allele[k]:
												IsSpec = False
									if IsSpec:
										if not (n in DicoSpecAll):
											DicoSpecAll[n] = set()
										DicoSpecAll[n].add(all)
						# Testing homozygosity and presence in accession
						for gp in DicoSpecAll:
							for all in DicoSpecAll[gp]:
								InAllAcc = True
								AlwaysHomo = True
								for ACC in dico_group[gp]:
									ACCESSION = data[header.index(ACC)].split(':')
									GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
									if not('.' in  GENOTYPE):
										# Testing homozygosity
										if all in GENOTYPE:
											if len(GENOTYPE) != 1:
												AlwaysHomo = False
										# Testing presence in all accessions
										if not(all in GENOTYPE):
											InAllAcc = False
								if options.all == 'y' and options.hom == 'y':
									if InAllAcc and AlwaysHomo:
										outfile.write('\t'.join((gp, CHR, POS, VARIANT[int(all)])))
										outfile.write('\n')
								elif options.all == 'y' and options.hom == 'n':
									if InAllAcc:
										outfile.write('\t'.join((gp, CHR, POS, VARIANT[int(all)])))
										outfile.write('\n')
								elif options.all == 'n' and options.hom == 'y':
									if AlwaysHomo:
										outfile.write('\t'.join((gp, CHR, POS, VARIANT[int(all)])))
										outfile.write('\n')
								elif options.all == 'n' and options.hom == 'n':
									outfile.write('\t'.join((gp, CHR, POS, VARIANT[int(all)])))
									outfile.write('\n')
								else:
									sys.exit('Unmannaged options... please review options --all and --hom. Only possible values are "n" and "y"')
						

		file.close()
	
if __name__ == "__main__": __main__()