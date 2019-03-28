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

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.stdout.write('modules loaded\n')


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
			list2return.append(liste_allele[int(n)])
	
	# returning list
	return list2return

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

def calcul_ident(LISTE1, LISTE2):
	
	"""
		Calculate identity between two genotypes
	"""
	
	min_ploidy = min(len(LISTE1), len(LISTE2))
	
	ident = 0
	if len(LISTE1) < len(LISTE2):
		for n in LISTE1:
			if n in LISTE2:
				ident += 1
				LISTE2.remove(n)
	else:
		for n in LISTE2:
			if n in LISTE1:
				ident += 1
				LISTE1.remove(n)
	return float(ident)/float(min_ploidy)
		
def draw_density_plot(CHR, SIZE, DICO, PDF):

	fig = plt.figure(figsize=(14.85, 10.5))
	taille = len(DICO)
	i = 0
	COLOR = ((0,0,1),(0,1,0),(1,0,0))
	
	for n in DICO:
		if i%5 == 0 and i != 0:
			PDF.savefig()
			plt.close(fig)
			fig = plt.figure(figsize=(14.85, 10.5))
			i = 0
		i += 1
		XCOORD = numpy.array(DICO[n]['pos'])
		YCOORD = numpy.array(DICO[n]['ident'])
		ax = plt.subplot(5,1,i)
		ax.scatter(XCOORD, YCOORD, linewidth=1, linestyle='-')
		
		j = 0
		liste_ploidy = sorted(DICO[n]['values'].keys())
		for ploidy in liste_ploidy:
			YLINE = DICO[n]['values'][ploidy]
			ax.plot(XCOORD, YLINE, linewidth=1, linestyle='-', color=COLOR[j], label=str(ploidy))
			j += 1
		
		ax.set_xlim(0, SIZE)
		ax.set_ylim(-0.1, 1.1)
		ax.set_ylabel(n, fontsize=12, fontweight='bold')
		ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		plt.title(CHR)
	
	PDF.savefig()
	plt.close(fig)
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-a',	'--acc',			dest='acc',			default=None,			help='Accession to compare to others. [Default: %default]')
	parser.add_option( '-c',	'--comp',			dest='comp',		default=None,			help='Accessions to be compared separated by ":". [Default: %default]')
	parser.add_option( '-w',	'--win',			dest='win',			default=0,				help='Half window (unit = variant site) to compute density along chromosomes. [Default: %default]')
	parser.add_option( '-o',	'--out',			dest='out',			default=None,			help='Prefix for output file names. [Default: %default]')

	(options, args) = parser.parse_args()
	
	WINDOW = int(options.win)
	
	liste_acc = options.comp.split(':')
	liste_acc.append(options.acc)
	
	dico_acc_geno = {}
	for accession in liste_acc:
		dico_acc_geno[accession] = {}
	
	# recording informations in a dictionnary
	dico_chr = {}
	file = open(options.vcf)
	for line in file:
		data = line.split()
		if data:
			if data[0][0:9] == '##contig=':
				info = get_chr_size(data)
				dico_chr[info[0]] = info[1]
				for accession in dico_acc_geno:
					dico_acc_geno[accession][info[0]] = {}
			elif data[0] == '#CHROM':
				header = data
			elif data[0][0] != "#":
				for accession in liste_acc:
					genotype = recup_geno(data, header, accession)
					if not('NA' in genotype[2:]):
						dico_acc_geno[accession][genotype[0]][genotype[1]] = list(genotype[2:])
	file.close()
	
	# Now it's time to calculate identity
	dico_ident = {}
	dico_values = set()
	for accession in dico_acc_geno:
		if accession != options.acc:
			dico_ident[accession] = {}
			for chr in dico_acc_geno[accession]:
				dico_ident[accession][chr] = {}
				for pos in dico_acc_geno[options.acc][chr]:
					if pos in dico_acc_geno[accession][chr]:
						identity = calcul_ident(list(dico_acc_geno[options.acc][chr][pos]), list(dico_acc_geno[accession][chr][pos]))
						dico_ident[accession][chr][pos] = identity
						dico_values.add(identity)
	
	# Now it's time to plot identity along chromosomes
	for acc in liste_acc:
		if acc != options.acc:
			outfile = open(options.out+'-'+acc+'.scatter.txt','w')
			outfile.close()
			for i in dico_values:
				outfile = open(options.out+'-'+acc+'.density'+str(i)+'.txt','w')
				outfile.close()
	
	
	list_chromosomes = sorted(dico_chr.keys())
	with PdfPages(options.out+'.pdf') as pdf:
		for chr in list_chromosomes:
			if len(dico_acc_geno[options.acc][chr].keys()) < WINDOW*2 + 1:
				sys.stdout.write("Chromosome "+chr+" is removed from statistics because of insufficient marker number\n")
			else:
				dico_final = {}
				for accession in dico_ident:
					dico_final[accession] = {}
					dico_final[accession]['pos'] = list(map(str, sorted(list(map(int, dico_ident[accession][chr].keys())))))
					dico_final[accession]['ident'] = []
					dico_final[accession]['values'] = {}
					for values in dico_values:
						dico_final[accession]['values'][values] = []
					liste_value = []
					for pos in dico_final[accession]['pos']:
						liste_value.append(dico_ident[accession][chr][pos])
						dico_final[accession]['ident'].append(dico_ident[accession][chr][pos])
						if len(liste_value) == (WINDOW*2 + 1):
							# dico_final[accession]['ident'].append(sum(liste_value)/float(len(liste_value)))
							for values in dico_values:
								dico_final[accession]['values'][values].append(liste_value.count(values)/float(len(liste_value)))
							
							del liste_value[0]
				
					# Finalizing data
					for values in dico_values:
						debut_value = dico_final[accession]['values'][values][0]
						fin_value = dico_final[accession]['values'][values][-1]
						dico_final[accession]['values'][values] = [debut_value]*WINDOW + dico_final[accession]['values'][values] + [fin_value]*WINDOW
					
				draw_density_plot(chr, dico_chr[chr], dico_final, pdf)
				
				for acc in options.comp.split(':'):
					outfile = open(options.out+'-'+acc+'.scatter.txt','a')
					for i in range(len(dico_final[acc]['pos'])):
						outfile.write('\t'.join(map(str,[chr, dico_final[acc]['pos'][i], dico_final[acc]['pos'][i], dico_final[acc]['ident'][i]])))
						outfile.write('\n')
						for k in dico_values:
							outfile2 = open(options.out+'-'+acc+'.density'+str(k)+'.txt','a')
							outfile2.write('\t'.join(map(str,[chr, dico_final[acc]['pos'][i], dico_final[acc]['pos'][i], dico_final[acc]['values'][k][i]])))
							outfile2.write('\n')
							outfile2.close()
					outfile.close()
					
if __name__ == "__main__": __main__()