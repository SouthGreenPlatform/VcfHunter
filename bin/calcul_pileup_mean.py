#!/usr/bin/env python
#
#  Copyright 2017 CIRAD
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
import optparse, os, shutil, subprocess, sys, tempfile, fileinput, operator, time, gzip

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def draw_plot(X, Y, OUT, TITLE, FILLUNDER, YLIM, NEGATIVE, V):
	fig = plt.figure(figsize=(15, 5))
	ax = plt.subplot2grid((1,1),(0,0))
	
	if NEGATIVE == 'y':
		ax.plot(X, Y, linewidth=2, color=(0.75, 0.31, 0.30))
		ax.plot(X, Y*-1, linewidth=2, color=(0.75, 0.31, 0.30))
		if FILLUNDER == 'y':
			ax.fill_between(X,Y*-1, Y, color=(0.75, 0.31, 0.30))
		ax.set_xlim(0, X[-1])
		if YLIM != None:
			ax.set_ylim(float(YLIM*-1), float(YLIM))
		else:
			ax.set_ylim(max(Y)*-1, max(Y))
	else:
		ax.plot(X, Y, linewidth=2, color=(0.75, 0.31, 0.30))
		if FILLUNDER == 'y':
			ax.fill_between(X,0, Y, color=(0.75, 0.31, 0.30))
		ax.set_xlim(0, X[-1])
		if YLIM != None:
			ax.set_ylim(0, float(YLIM))
		else:
			ax.set_ylim(0, max(Y))
	ax.set_title(TITLE, fontweight='bold', position=(0.5, 1))
	
	
	fig.savefig(OUT)
	plt.close(fig)
	
	os.system("pmap %s | grep total | sed 's/ total/RAM used:/' | sed 's:  \+: :'" % V)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program count from a pileup like file the average coverage of covered"
	"\npositions in a sliding window.")
	# Wrapper options. 
	parser.add_option( '-p', '--pileup',		dest='pileup',			default=None,		help='The pileup file (tabulated)')
	parser.add_option( '-w', '--window',		dest='window',			default='100000',	help='The window to calculate proportion')
	parser.add_option( '-f', '--fasta',			dest='fasta',			default=None,		help='The multifasta reference file')
	parser.add_option( '-F', '--FillUnder',		dest='FillUnder',		default='n',		help='Fill the curve under, (y or n)')
	parser.add_option( '-y', '--Ylim',			dest='Ylim',			default=None,		help='Set Y limit. Can be omitted')
	parser.add_option( '-c', '--chr',			dest='chr',				default=None,		help='Chromosome to work with. They should be separated by ":". If omitted, all chromosomes will be used')
	parser.add_option( '-T', '--title',			dest='title',			default='n',		help='Add title to figure: y or n')
	parser.add_option( '-n', '--negative',		dest='negative',		default='n',		help='Draw also negative curve: y or n')
	parser.add_option( '-o', '--out',			dest='out',				default='',			help='Prefix for the output files')
	parser.add_option( '-t', '--outtype',		dest='outtype',			default='svg',		help='Output file type: png or svg')
	parser.add_option( '-d', '--draw',			dest='draw',			default='y',		help='Draw a figure. Possible values "y" or "n".')
	(options, args) = parser.parse_args()
	V = os.getpid()
	
	WIN = int(options.window)
	
	if options.chr == None:
		ChrToWorkWith = []
	else:
		ChrToWorkWith = options.chr.split(':')
	
	# Getting sequence statistics
	sys.stdout.write("Loading fasta sequence\n")
	dico_chr = {}
	sequence_dict = SeqIO.index(options.fasta, "fasta")
	if ChrToWorkWith:
		for n in ChrToWorkWith:
			dico_chr[n] = len(str(sequence_dict[n].seq))
	else:
		for n in sequence_dict:
			dico_chr[n] = len(str(sequence_dict[n].seq))
	os.system("pmap %s | grep total | sed 's/ total/RAM used:/' | sed 's:  \+: :'" % V)
	del sequence_dict
	
	# Working chromosome per chromosome
	listChr = sorted(dico_chr.keys())
	outfile = gzip.open(options.out+'.tab.gz', 'wt')
	for ch in listChr:
		sys.stdout.write("Working on chromosome "+str(ch)+" ...\n")
		sys.stdout.flush()
		ListVal = [0]*dico_chr[ch]
		ListPos = [0]*dico_chr[ch]
		# Filling the list
		file = open(options.pileup)
		for line in file:
			data = line.split()
			if data:
				chro = data[0]
				if chro == ch:
					ListVal[int(data[1])-1] = float(data[2])
					ListPos[int(data[1])-1] = 1
		# Calculating proportion on windows
		debut = 0
		fin = WIN
		chr_fin = len(ListVal)-1
		list_pos = [1]
		list_cov = []
		while debut < chr_fin:
			temp_fin = fin
			if chr_fin < fin:
				temp_fin = chr_fin + 1
			Value = float(sum(ListVal[debut:temp_fin]))
			PosNum = sum(ListPos[debut:temp_fin])
			if PosNum > 0:
				MeanCov = Value/PosNum
			else:
				MeanCov = 0
			outfile.write('\t'.join([ch, str(debut+1), str(temp_fin), str(MeanCov), str(PosNum)])+'\n')
			if len(list_cov) == 0:
				list_cov.append(MeanCov)
			list_cov.append(MeanCov)
			list_pos.append(temp_fin)
			debut += WIN
			fin += WIN
			sys.stdout.flush()
		if options.draw == 'y':
			if options.title == 'y':
				draw_plot(np.array(list_pos), np.array(list_cov), options.out+ch+'.'+options.outtype, options.out+'-'+ch, options.FillUnder, options.Ylim, options.negative, V)
			else:
				draw_plot(np.array(list_pos), np.array(list_cov), options.out+ch+'.'+options.outtype, ch, options.FillUnder, options.Ylim, options.negative, V)
		os.system("pmap %s | grep total | sed 's/ total/RAM used:/' | sed 's:  \+: :'" % V)
		

if __name__ == "__main__": __main__()



