#!/usr/bin/env python
#
#  Copyright 2023 CIRAD
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
import gzip


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(
		description="This program select sites in a vcf according to a list of positions.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser.add_option( '-v',  '--vcf',		dest='vcf',		default=None,	help='Vcf file. [Default: %default]')
	parser.add_option( '-s',  '--sites',	dest='sites',	default=None,	help='Sites to keep in a 2 column file (c1 = chromosome, c2 = position). [Default: %default]')
	parser.add_option( '-o',  '--out',		dest='out',		default=None,	help='Output file name. Output will be gzipped. [Default: %default]')
	(options, args) = parser.parse_args()

	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.sites == None:
		sys.exit('Please provide a file to --sites argument')
	if options.out == None:
		sys.exit('Please provide a output name to --out argument')
	
	# Recording sites to keep
	dico = {}
	file = open(options.sites)
	for line in file:
		data = line.split()
		if data:
			if not data[0] in dico:
				dico[data[0]] = set()
			dico[data[0]].add(data[1])
	file.close()
	
	# Parsing the vcf file
	VCF = options.vcf
	
	outfile = gzip.open(options.out,'wt')
	
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == "#":
				if data[0] == "#CHROM":
					CHROMPOS = data.index("#CHROM")
					POPPOS = data.index("POS")
				outfile.write(line)
			else:
				chr = data[CHROMPOS]
				pos = data[POPPOS]
				if chr in dico:
					if pos in dico[chr]:
						outfile.write(line)
	

if __name__ == "__main__": __main__()
