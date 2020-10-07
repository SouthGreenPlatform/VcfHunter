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

	
import gzip
import optparse

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n")
	# Wrapper options.
	parser.add_option( '-c',	'--haplo',	dest='haplo',	default=None,				help='A three column file containing on column 1: chromosome, column 2 : position, column 3 : allele')
	parser.add_option( '-v',	'--vcf',	dest='vcf',		default=None,				help='A vcf file containing accessions to treat')
	parser.add_option( '-o',	'--out',	dest='out',		default='HaploProp.tab',	help='Output file name')
	
	(options, args) = parser.parse_args()
	
	if options.haplo == None:
		sys.exit('Please provide a haplotype file to --haplo argument')
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	
	DicoMarker = {}
	file = open(options.haplo, 'r')
	for line in file:
		data = line.split()
		if data:
			if not(data[0]) in DicoMarker:
				DicoMarker[data[0]] = {}
			DicoMarker[data[0]][data[1]] = data[2]
	
	file.close()
	
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
	else:
		file = open(options.vcf)
	
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				dicoAccPos = {}
				dicoAccStat = {}
				ChrPos = data.index("#CHROM")
				PosPos = data.index("POS")
				RefPos = data.index("REF")
				AltPos = data.index("ALT")
				FormatPos = data.index("FORMAT")
				for acc in data[FormatPos+1:]:
					dicoAccPos[acc] = data.index(acc)
					dicoAccStat[acc] = [0,0]
			elif data[0][0] != "#":
				Chr = data[ChrPos]
				if Chr in DicoMarker:
					Pos = data[PosPos]
					if Pos in DicoMarker[Chr]:
						Allele = [data[RefPos]]+data[AltPos].split(',')
						GoodAllele = str(Allele.index(DicoMarker[Chr][Pos]))
						FORMAT = data[FormatPos].split(':')
						GTPos = FORMAT.index('GT')
						for acc in dicoAccPos:
							AccGeno = data[dicoAccPos[acc]].split(':')[GTPos].split('/')
							if not('.' in AccGeno):
								if GoodAllele in AccGeno:
									dicoAccStat[acc][0] += 1
								dicoAccStat[acc][1] += 1
	file.close()
	
	outfile = open(options.out,'w')
	for acc in dicoAccStat:
		if dicoAccStat[acc][1] != 0:
			outfile.write('\t'.join([acc, str(dicoAccStat[acc][0]/dicoAccStat[acc][1]), str(dicoAccStat[acc][0]), str(dicoAccStat[acc][1])]))
			outfile.write('\n')
	outfile.close()
	exit()



if __name__ == "__main__": __main__()
