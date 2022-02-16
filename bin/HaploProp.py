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
	"This program calculate the proportion of alleles of an haplotype found in an accession genotyped in a vcf file.")
	# Wrapper options.
	parser.add_option( '-c',	'--haplo',	dest='haplo',	default=None,				help='A three column file containing on column 1: chromosome, column 2 : position, column 3 : allele')
	parser.add_option( '-v',	'--vcf',	dest='vcf',		default=None,				help='A vcf file containing accessions to treat')
	parser.add_option( '-l',	'--loc',	dest='loc',		default=None,				help='Optional: locate specific region in which alleles will be searched. Should be formated as follows: chr,start,end')
	parser.add_option( '-o',	'--out',	dest='out',		default='HaploProp.tab',	help='Output file name')
	
	(options, args) = parser.parse_args()
	
	if options.haplo == None:
		sys.exit('Please provide a haplotype file to --haplo argument')
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	
	
	
	DicoMarker = {}
	file = open(options.haplo, 'r')
	if options.loc != None:
		CHR, Start, End = options.loc.split(',')
		Start = int(Start)
		End = int(End)
		for line in file:
			data = line.split()
			if data:
				if data[0] == CHR:
					position = int(data[1])
					if Start <= position and position <= End:
						if not(data[0]) in DicoMarker:
							DicoMarker[data[0]] = {}
						DicoMarker[data[0]][data[1]] = data[2]
						# print(data[0], data[1])
	else:
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
				dicoAccReadStat = {}
				ChrPos = data.index("#CHROM")
				PosPos = data.index("POS")
				RefPos = data.index("REF")
				AltPos = data.index("ALT")
				FormatPos = data.index("FORMAT")
				for acc in data[FormatPos+1:]:
					dicoAccPos[acc] = data.index(acc)
					dicoAccStat[acc] = [0,0]
					dicoAccReadStat[acc] = [0,0]
			elif data[0][0] != "#":
				Chr = data[ChrPos]
				if Chr in DicoMarker:
					Pos = data[PosPos]
					if Pos in DicoMarker[Chr]:
						Allele = [data[RefPos]]+data[AltPos].split(',')
						if DicoMarker[Chr][Pos] in Allele: ## If the diagnostic allele is in the vcf
							GoodAllele = str(Allele.index(DicoMarker[Chr][Pos]))
							FORMAT = data[FormatPos].split(':')
							GTPos = FORMAT.index('GT')
							ADPos = FORMAT.index('AD')
							DPPos = FORMAT.index('DP')
							for acc in dicoAccPos:
								AccGeno = data[dicoAccPos[acc]].split(':')[GTPos].split('/')
								ADGeno = data[dicoAccPos[acc]].split(':')[ADPos].split(',')
								DPGeno = data[dicoAccPos[acc]].split(':')[DPPos]
								if not('.' in AccGeno):
									if GoodAllele in AccGeno:
										dicoAccStat[acc][0] += 1
									dicoAccStat[acc][1] += 1
									dicoAccReadStat[acc][0] += int(ADGeno[int(GoodAllele)])
									dicoAccReadStat[acc][1] += int(DPGeno)
						else: # The diagnostic allele is not in vcf, thus individuals do not have this allele
							for acc in dicoAccPos:
								AccGeno = data[dicoAccPos[acc]].split(':')[GTPos].split('/')
								if not('.' in AccGeno):
									dicoAccStat[acc][1] += 1
									dicoAccReadStat[acc][0] += int(ADGeno[int(GoodAllele)])
									dicoAccReadStat[acc][1] += int(DPGeno)
							
	file.close()
	
	outfile = open(options.out,'w')
	outfile.write('\t'.join(["name", "SiteRatio", "SiteWithSpecAllele", "TotalSites", "ReadRatio", "ReadWithSpecAllele", "TotalReads"]))
	outfile.write('\n')
	for acc in dicoAccStat:
		if dicoAccStat[acc][1] != 0:
			outfile.write('\t'.join([acc, str(dicoAccStat[acc][0]/dicoAccStat[acc][1]), str(dicoAccStat[acc][0]), str(dicoAccStat[acc][1]), str(dicoAccReadStat[acc][0]/dicoAccReadStat[acc][1]), str(dicoAccReadStat[acc][0]), str(dicoAccReadStat[acc][1])]))
			outfile.write('\n')
	outfile.close()
	exit()



if __name__ == "__main__": __main__()
