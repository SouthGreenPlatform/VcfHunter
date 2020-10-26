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
import gzip
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
import optparse

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\nThis program looks for specific alleles of a consensus haplotype found in one specific group of accessions and absent "
	"from another group of accessions. The read allelic ratio of these specific alleles was then searched and plotted along "
	"chromosome for each studied accession. This program only work on one chromosome and several chromosomes should not be "
	"passed in the consensus.")
	# Wrapper options.
	parser.add_option( '-H', '--haplo',		dest='haplo',	default=None, help='A three column file containing on column 1: chromosome, column 2 : position, column 3 : allele')
	parser.add_option( '-W', '--With',		dest='With',	default=None, help='A one column file containing accessions with the haplotype')
	parser.add_option( '-w', '--Without',	dest='Without',	default=None, help='A one column file containing accessions without the haplotype')
	parser.add_option( '-a', '--acc',		dest='acc',		default=None, help='A one column file containing accessions to treat')
	parser.add_option( '-v', '--vcf',		dest='vcf',		default=None, help='A vcf file containing accessions passed to "--With" and "--Without" options')
	parser.add_option( '-V', '--VCF',		dest='VCF',		default=None, help='A vcf file containing accessions to treat')
	parser.add_option( '-o', '--out',		dest='out',		default=None, help='Prefix for output file names')
	
	(options, args) = parser.parse_args()
	
	if options.haplo == None:
		sys.exit('Please provide a haplotype file to --haplo argument')
	if options.With == None:
		sys.exit('Please provide a file to --With argument')
	if options.Without == None:
		sys.exit('Please provide a file to --Without argument')
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.VCF == None:
		sys.exit('Please provide a vcf file to --VCF argument')
	
	
	# Loading markers to treat
	DicoMarker = {}
	file = open(options.haplo)
	for line in file:
		data = line.split()
		if data:
			if not(data[0]) in DicoMarker:
				DicoMarker[data[0]] = {}
			DicoMarker[data[0]][data[1]] = data[2]
	file.close()
	
	# Loading accessions to work with
	ACCTOTREAT = []
	if options.acc != None:
		file = open(options.acc)
		for line in file:
			data = line.split()
			if data:
				ACCTOTREAT.append(data[0])
		file.close()
	
	# Selecting marker present in accessions from the group, and absent from the other group
	ListAvec = set()
	file = open(options.With)
	for line in file:
		data = line.split()
		if data:
			ListAvec.add(data[0])
	file.close()
	
	ListeSans = set()
	file = open(options.Without)
	for line in file:
		data = line.split()
		if data:
			ListeSans.add(data[0])
	file.close()
	
	
	
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
	else:
		file = open(options.vcf)
	NewDicoMarker = {}
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				dicoAccPos = {}
				ChrPos = data.index("#CHROM")
				PosPos = data.index("POS")
				RefPos = data.index("REF")
				AltPos = data.index("ALT")
				FormatPos = data.index("FORMAT")
				for acc in data[FormatPos+1:]:
					dicoAccPos[acc] = data.index(acc)
			elif data[0][0] != "#":
				Chr = data[ChrPos]
				if Chr in DicoMarker:
					Pos = data[PosPos]
					if Pos in DicoMarker[Chr]:
						Allele = [data[RefPos]]+data[AltPos].split(',')
						GoodAllele = str(Allele.index(DicoMarker[Chr][Pos]))
						FORMAT = data[FormatPos].split(':')
						GTPos = FORMAT.index('GT')
						DicoAvec = set()
						DicoSans = set()
						for acc in dicoAccPos:
							AccGeno = set(data[dicoAccPos[acc]].split(':')[GTPos].split('/'))
							if acc in ListAvec:
								DicoAvec.update(AccGeno)
							elif acc in ListeSans:
								DicoSans.update(AccGeno)
						if DicoAvec and DicoSans:
							if GoodAllele in DicoAvec-DicoSans:
								if not(Chr in NewDicoMarker):
									NewDicoMarker[Chr] = {}
								NewDicoMarker[Chr][Pos] = DicoMarker[Chr][Pos]
	
	file.close()
	
	if len(NewDicoMarker) != 1:
		sys.exit('Warning, their is more than one chromosome in the haplotype file. This is not allowed... The program stopped before finishing\n')
	for Chr in NewDicoMarker:
		sys.stdout.write('Number of specific alleles on '+Chr+': '+str(len(NewDicoMarker[Chr]))+'\n')
		Chromosome = Chr
	
	if options.VCF[-3:] == '.gz':
		file = gzip.open(options.VCF,'rt')
	else:
		file = open(options.VCF)
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				dicoAccPos = {}
				dicoAccStat = {}
				dicoPos = {}
				ChrPos = data.index("#CHROM")
				PosPos = data.index("POS")
				RefPos = data.index("REF")
				AltPos = data.index("ALT")
				FormatPos = data.index("FORMAT")
				if ACCTOTREAT:
					for acc in data[FormatPos+1:]:
						if acc in ACCTOTREAT:
							dicoAccPos[acc] = data.index(acc)
							dicoAccStat[acc] = [[],[],[]]
				else:
					for acc in data[FormatPos+1:]:
						dicoAccPos[acc] = data.index(acc)
						dicoAccStat[acc] = [[],[],[]]
						ACCTOTREAT.append(acc)
			elif data[0][0] != "#":
				Chr = data[ChrPos]
				if Chr in NewDicoMarker:
					Pos = data[PosPos]
					if Pos in NewDicoMarker[Chr]:
						if not(Chr in dicoPos):
							dicoPos[Chr] = []
						dicoPos[Chr].append(Pos)
						Allele = [data[RefPos]]+data[AltPos].split(',')
						GoodAllele = Allele.index(NewDicoMarker[Chr][Pos])
						GOODAllele = set()
						GOODAllele.add(GoodAllele)
						FORMAT = data[FormatPos].split(':')
						GTPos = FORMAT.index('GT')
						ADPos = FORMAT.index('AD')
						DPPos = FORMAT.index('DP')
						for acc in ACCTOTREAT:
							GENO = data[dicoAccPos[acc]].split(':')[GTPos].split('/')
							if '.' in GENO:
								dicoAccStat[acc][0].append(numpy.nan)
							else:
								AccGeno = list(map(int, GENO))
								ADGeno = list(map(int, data[dicoAccPos[acc]].split(':')[ADPos].split(',')))
								DPGeno = int(data[dicoAccPos[acc]].split(':')[DPPos])
								if not('.' in AccGeno):
									if GoodAllele in AccGeno:
										dicoAccStat[acc][0].append(ADGeno[GoodAllele]/DPGeno)
										NewSet = list(set(AccGeno) - GOODAllele)
										if len(NewSet) == 1:
											dicoAccStat[acc][1].append(ADGeno[NewSet[0]]/DPGeno)
											dicoAccStat[acc][2].append(0)
										elif len(NewSet) == 0:
											dicoAccStat[acc][1].append(0)
											dicoAccStat[acc][2].append(0)
										else:
											sys.exit('There is a bug')
									else:
										NewSet = list(set(AccGeno))
										if len(NewSet) == 1:
											dicoAccStat[acc][0].append(0)
											dicoAccStat[acc][1].append(ADGeno[NewSet[0]]/DPGeno)
											dicoAccStat[acc][2].append(0)
										elif len(NewSet) == 2:
											dicoAccStat[acc][0].append(0)
											dicoAccStat[acc][1].append(ADGeno[NewSet[0]]/DPGeno)
											dicoAccStat[acc][2].append(ADGeno[NewSet[1]]/DPGeno)
										else:
											sys.exit('There is a bug')
		
	# Color definition
	color = ((1,0,0),(0,0.8,0))
	
	NB=15
	
	# Calcul de la taille max
	MAXI = 29070452
	
	# Calculating ticks
	i = 5000000
	ticks_pos = [0]
	ticks_labels = [0]
	
	while i < MAXI:
		ticks_pos.append(i)
		ticks_labels.append(str(int(i/1000000)))
		i += 5000000
	
	i = 1000000
	
	minor_ticks = []
	while i < MAXI:
		minor_ticks.append(i)
		i += 1000000
	
	######Drawing allele ratio######
	
	# On fait la figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)
	
	# AccToPlot=list(sorted(dicoAccStat.keys()))
	
	outfile = open(options.out+'.tab','w')
	
	for i, acc in enumerate(ACCTOTREAT):
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1)
		ax.set_xlim(0, MAXI)
		ax.axhline(y=0.33, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.66, linewidth=0.1, color = (0,0,0.5))
		ax.plot(dicoPos[Chromosome], dicoAccStat[acc][0], 'o', ms=5, mew=0, mfc=color[0])
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off') 
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, acc, size=12, va='center', fontweight='bold')
		POSSPAN += 1 
		if (i+1)%15==0:   
			OUT2=options.out+'_'+Chromosome+'_'+str(int(i/15)+1)+'_Ratio.png'
			fig.savefig(OUT2)
			plt.close(fig)
			POSSPAN = 0
			fig = plt.figure(figsize=(10.5, 14.85))
			fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)
		j=i
		outfile.write(acc+'\t'+str(numpy.nanmean(dicoAccStat[acc][0]))+'\t'+str(numpy.nanmedian(dicoAccStat[acc][0]))+'\t'+str(numpy.count_nonzero(~numpy.isnan(dicoAccStat[acc][0])))+'\n')
	
	if (i+1)%15 != 0:
		fig.savefig(options.out+'_'+Chromosome+'_'+str(int(j/15)+1)+'_Ratio.png')
		plt.close(fig)
	
	outfile.close()


if __name__ == "__main__": __main__()
