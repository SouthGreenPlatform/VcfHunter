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

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap


def draw_chr(DICO_INFO, CHR_INFO, OUT, MEAN_COV, PLOIDY, halfWindow, PSIZE, LSIZE, VERT, VERTREG):

	# Obtaining chromosomes list
	chr_list = sorted(list(DICO_INFO.keys()))

	# Calculating chromosome numbre
	NB = len(chr_list)

	# Calculating max chromosome size
	MAXI = 0
	for n in CHR_INFO:
		if MAXI < CHR_INFO[n]:
			MAXI = CHR_INFO[n]

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


		######Drawing coverage######

	# Initiating Figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)

	# Drawing per chromosomes
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 2*MEAN_COV)
		ax.set_xlim(0, MAXI)
		position = sorted(list(DICO_INFO[chr].keys()))
		value =  []
		final_pos = []
		fill_plus1 = [(PLOIDY+0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_moins1 = [(PLOIDY-0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_plus2 = [(PLOIDY+1.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_moins2 = [(PLOIDY-1.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		
		if chr in VERT:
			for pos in VERT[chr]:
				ax.axvline(x=pos, ymin=0, ymax = 2*MEAN_COV, linewidth=2, color='black')
		
		if chr in VERTREG:
			for pos in VERTREG[chr]:
				ax.fill_betweenx([0,2*MEAN_COV], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)

		for pos in position:
			value.append(DICO_INFO[chr][pos])
			final_pos.append(pos)

		ValueNorm = []
		PosNorm = []
		Window = 1+2*halfWindow
		for n in range(len(value)-Window):
			IntermVal = value[n:n+Window]
			ValueNorm.append(sum(IntermVal)/len(IntermVal))
			PosNorm.append(final_pos[n+halfWindow])

		ax.fill_between([0]+final_pos+[CHR_INFO[chr]], fill_plus1, fill_moins1, color=(1,0,0), alpha=0.20) ## For band
		ax.axhline(y=MEAN_COV, linewidth=0.5, color = (0,0,0))
		ax.axhline(y=(PLOIDY+1)/PLOIDY*MEAN_COV, linewidth=0.5, color = 'chartreuse')
		ax.axhline(y=(PLOIDY-1)/PLOIDY*MEAN_COV, linewidth=0.5, color = 'chartreuse')
		ax.axhline(y=(PLOIDY+2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.axhline(y=(PLOIDY-2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.plot(final_pos, value, 'o', ms=PSIZE, mew=0, mfc='blue')
		ax.plot(PosNorm, ValueNorm, color='red', lw=LSIZE)
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])

		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)

	# Adding chromosome name
	POSSPAN = 0
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', fontweight='bold')

		POSSPAN += 1

	fig.savefig(OUT)
	plt.close(fig)

def get_chr_size(LINE):

		"""
			Get reference chromosome size from vcf file

			:param LINE: The ##contig line of the vcf file splitted on spaces
			:type LINE: list
			:return: A list with [0] --> chromosome name, [1] --> chromosome size
			:rtype: list
		"""

		split_on_chev = LINE[0].split('>')[0]

		split_on_eq = split_on_chev.split('=')

		return [split_on_eq[2].split(',')[0], int(split_on_eq[3])]

def __main__():
		#Parse Command Line
		parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\nThis program draw coverage from accessions found in a vcf file and output figures in png format. Several options are available to custom figures.")
		# Wrapper options.
		parser.add_option( '',  '--vcf',	dest='vcf',		default=None,	help='Vcf file. [Default: %default]')
		parser.add_option( '',  '--acc',	dest='acc',		default=None,	help='Accession list in a 2 column file (c1 = acc name, c2 = ploidy). [Default: %default]')
		parser.add_option( '',  '--minCov',	dest='minCov',	default='10',	help='Minimal coverage to keep the data point [Default: %default]')
		parser.add_option( '',  '--maxCov',	dest='maxCov',	default='100',	help='Maximal coverage to keep the data point [Default: %default]')
		parser.add_option( '',  '--psize',	dest='psize',	default='1',	help='Dot size in graph. [Default: %default]')
		parser.add_option( '',  '--lsize',	dest='lsize',	default='1',	help='Size of the line of the mean value curve. [Default: %default]')
		parser.add_option( '',  '--win',	dest='halfwin',	default='10',	help='Size of half sliding window that allow to draw mean value curve [Default: %default]')
		parser.add_option( '',  '--loc',	dest='loc',		default='',		help='Regions to locate by vertical line. This should be formated this way: Chromosome_name:position,chromosome_name:position, ... [Default: %default]')
		parser.add_option( '',  '--out',	dest='out',		default=None,	help='Prefix for file name output. [Default: %default]')
		(options, args) = parser.parse_args()

		if options.vcf == None:
				sys.exit('Please provide a vcf file to --vcf argument')
		if options.acc == None:
				sys.exit('Please provide a accession name to --acc argument')
		if options.out == None:
				sys.exit('Please provide a output name to --out argument')

		MINCOV = int(options.minCov)
		MAXCOV = int(options.maxCov)
		
		VERT = {}
		VERTREG = {}
		for n in options.loc.split(':'):
			info = n.split(',')
			if len(info) == 1:
				pass
			else:
				chr = info[0]
				if not(chr in VERT):
					VERT[chr] = set()
					VERTREG[chr] = []
				if len(info) == 2:
					VERT[chr].add(int(info[1]))
				elif len(info) == 3:
					VERTREG[chr].append(info[1:])
				else:
					sys.stdout.write('Wrong format passed to --loc argument\n')
			

		ACCLIST = open(options.acc)
		for AccLine in ACCLIST:
			data = AccLine.split()
			if data:
				# Recording accession name
				ACCESS = data[0]
				PLOIDY = int(data[1])

				# For linear drawing
				dico_draw = {}
				dico_chr = {}
				total = []
				if options.vcf[-3:] == '.gz':
					file = gzip.open(options.vcf,'rt')
				else:
					file = open(options.vcf)
				for line in file:
					data = line.split()
					if data:
						if data[0] == "#CHROM":
								header = data
						elif data[0][0] != "#":
							FORMAT = data[header.index("FORMAT")].split(':')
							CHR = data[header.index("#CHROM")]
							POS = int(data[header.index("POS")])
							REF = [data[header.index("REF")]]
							ALT = data[header.index("ALT")].split(',')
							VARIANT = REF + ALT
							ACCESSION = data[header.index(ACCESS)].split(':')
							COVERAGE = int(ACCESSION[FORMAT.index("DP")])
							if COVERAGE >= MINCOV and COVERAGE <= MAXCOV:
								if not(CHR in dico_draw):
									dico_draw[CHR] = {}
								dico_draw[CHR][POS] = COVERAGE
								total.append(COVERAGE)
						elif data[0][0:9] == '##contig=':
							chr_info = get_chr_size(data)
							dico_chr[chr_info[0]] = chr_info[1]
				file.close()
				sys.stdout.write('Average coverage : '+str(sum(total)/len(total))+'\n')
				sys.stdout.write('Median coverage : '+str(np.median(total))+'\n')
				# It's time to draw
				draw_chr(dico_draw, dico_chr, options.out+'_'+ACCESS+'.png', sum(total)/len(total), PLOIDY, int(options.halfwin), float(options.psize), float(options.lsize), VERT, VERTREG)

if __name__ == "__main__": __main__()
