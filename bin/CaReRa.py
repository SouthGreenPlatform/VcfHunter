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
# import os
# import shutil
# import subprocess
# import tempfile
# import fileinput
# import time
# import random
# import math
# import datetime
# import traceback
# import multiprocessing as mp
# from inspect import currentframe, getframeinfo
# from operator import itemgetter

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.patches as mpatches
from scipy.interpolate import spline
# from textwrap import wrap


def draw_chr(DICO_INFO, CHR_INFO, DICOREGIONHIGH, MAXx, OUT):


	# getting chromosomes list
	chr_list = sorted(list(DICO_INFO.keys()))

	# On calcule le nombre de fenetres
	NB = len(chr_list)

	# Calcule de la taille max des chromosomes
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
	
	Yticks_pos = []
	Yticks_labels = []
	Val = MAXx / 4
	count = 0
	for n in str(Val):
		if n != '.':
			count += 1
			if n != '0':
				break
		else:
			count = 0
	PAS = round(Val, count)
	Yticks_pos = [round(PAS, count), round(PAS*2, count), round(PAS*3, count)]
	Yticks_labels = list(map(str, Yticks_pos))


		######Drawing coverage######

	# On fait la figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)

	# On dessine par chromosomes
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1*MAXx)
		ax.set_xlim(0, MAXI)
		
		if chr in DICOREGIONHIGH:
			for bloc in DICOREGIONHIGH[chr]:
				ax.fill_betweenx([0,1.1*MAXx], bloc[0], bloc[1], color=(1,0,0), alpha=0.20)
		
		ax.plot(DICO_INFO[chr][1], DICO_INFO[chr][0], color='red', lw=2)
		ax.plot(DICO_INFO[chr][2], DICO_INFO[chr][3], color='blue', lw=2)
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])
		
		ax.set_yticks(Yticks_pos)
		ax.set_yticklabels(Yticks_labels)


		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)

	# On mets le nom des chromosomes
	POSSPAN = 0
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', ha='right', fontweight='bold')

		POSSPAN += 1

	fig.savefig(OUT+'.png')
	fig.savefig(OUT+'.svg')
	plt.close(fig)

def recode_matrix(LISTE):
	
	name = LISTE[0]
	coding = LISTE[1].split(',')
	if coding == ['nn','np']:
		LISTE_recode = recode_nn_np(LISTE[4:])
	elif coding == ['ll','lm']:
		LISTE_recode = recode_ll_lm(LISTE[4:])
	elif coding == ['hh','hk', 'kk']:
		LISTE_recode = recode_hh_hk_kk(LISTE[4:])
	elif coding == ['hh','k-']:
		LISTE_recode = recode_hh_k(LISTE[4:])
	else:
		sys.exit('Bad file format: The matrix file should contain 3 or 4 header columns\n'+LINE)
	return (name, numpy.array(LISTE_recode), coding)

def recode_hh_k(LISTE):
	liste = []
	for x in LISTE:
		if x == 'hh':
			liste.append(numpy.nan)
			liste.append(0)
		elif x == 'k-':
			liste.append(numpy.nan)
			liste.append(1)
		else:
			liste.append(numpy.nan)
			liste.append(numpy.nan)
	return liste

def recode_nn_np(LISTE):
	liste = []
	for x in LISTE:
		if x == 'nn':
			liste.append(numpy.nan)
			liste.append(0)
		elif x == 'np':
			liste.append(0)
			liste.append(1)
		else:
			liste.append(numpy.nan)
			liste.append(numpy.nan)
	return liste

def recode_ll_lm(LISTE):
	liste = []
	for x in LISTE:
		if x == 'll':
			liste.append(numpy.nan)
			liste.append(0)
		elif x == 'lm':
			liste.append(0)
			liste.append(1)
		else:
			liste.append(numpy.nan)
			liste.append(numpy.nan)
	return liste
	
def recode_hh_hk_kk(LISTE):
	liste = []
	for x in LISTE:
		if x == 'hh':
			liste.append(0)
			liste.append(0)
		elif x == 'hk':
			liste.append(0)
			liste.append(1)
		elif x == 'kk':
			liste.append(1)
			liste.append(1)
		else:
			liste.append(numpy.nan)
			liste.append(numpy.nan)
	return liste

def get_chr_size(FASTA):

	"""
		Get reference chromosome size from multifasta file

		:param FASTA: A multifasta file
		:type FASTA: fasta
		:return: A dictionary containing with key = sequence name and value = sequence size
		:rtype: list
	"""

	dico = {}
	SeqSize = 0
	file = open(FASTA)
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '>':
				if SeqSize == 0:
					SeqName = data[0][1:]
					SeqSize = 0
				else:
					dico[SeqName] = SeqSize
					SeqName = data[0][1:]
					SeqSize = 0
			else:
				if len(data) != 1:
					sys.exit('Oups, a space or a tabulation has been found in the multifasta file... This is not allowed and the program has crashed.')
				SeqSize += len(data[0])
	return dico

def get_win(DICO, WIN):

	"""
		Generate windows in which recombination rate will be calculated

		:param DICO: A dictionary containing file 
		:type DICO: dictionary
		:param WIN: Window size
		:type WIN: integer
		:return: A dictionary : {chr}->{pos}->[start, end, count]
		:rtype: list
	"""

	dicowin = {}
	for n in DICO:
		start = 1
		dicowin[n] = {}
		while start + WIN < DICO[n]:
			dicowin[n][start] = [start, (start+WIN)-1, 0, False]
			start += WIN
		dicowin[n][start] = [start, DICO[n], 0, False]
		
	return dicowin

def __main__():
		#Parse Command Line
		parser = optparse.OptionParser(usage='python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n This program calculate, on a given window, the recombination rate observed\n')
		# Wrapper options.
		parser.add_option( '',  '--fasta',	dest='fasta',	default=None,			help='The multifasta file containing reference sequence on which tags were located.')
		parser.add_option( '',  '--mat',	dest='mat',		default=None,			help='Phased and corrected matrix file. col 1 : markers name, col marker coding ("nn,np" or "hh,hk,kk" or "hh,k-" or "ll,lm"), 3, 4 necessary but not used by the program,'
		'col 5 to end : individual genotypes. First line contain header. Redundant names (markers and individuals) are not allowed')
		parser.add_option( '',  '--win',	dest='win',		default='500000',		help='The windows (in bases) in which the number of recombinations will be calculated. [Default: %default]')
		parser.add_option( '',  '--chr',	dest='chr',		default='all',			help='List of chromosomes to draw separated by ",". If omitted, all chromosomes will be drawn. [Default: %default]')
		parser.add_option( '',  '--excl',	dest='excl',	default=None,			help='List of marker couple to exclude. A couple of marker should have this structure (Name1:Name2) '
		'and each couple should be separated by ",".')
		parser.add_option( '',  '--rmun',	dest='rmun',	default='n',			help='Remove regions in which recombination rate has been estimated based on markers of surrounding windows. '
		'Values: "y" or "n" [Default: %default]')
		parser.add_option( '',  '--prefix',	dest='prefix',	default='RecombRate',	help='Prefix for the output files. [Default: %default]')
		(options, args) = parser.parse_args()

		if options.fasta == None:
			sys.exit('Please provide a fasta file to --fasta argument')
		if options.mat == None:
			sys.exit('Please provide a marker matrix file to --mat argument')
		
		EXCLCouple = {}
		if options.excl != None:
			for n in options.excl.split(','):
				m1, m2 = n.split(':')
				EXCLCouple[m1] = m2
		
		WIN = int(options.win)
		
		# Obtaining chromosome size
		chr_info = get_chr_size(options.fasta)
		
		# Calculating windows
		DicoWin = get_win(chr_info, int(WIN))
		
		# Calculating recombination rate
		file = open(options.mat)
		header = file.readline().split()
		IndNum = len(header[4:])
		print('Individuals in the population :', IndNum)
		
		newLine = file.readline().split()
		for line in file:
			preLine = newLine[:]
			newLine = line.split()
			
			if len(preLine) != len(newLine):
				sys.exit('Oups, the program has crashed... Lines in the file do not have the same length\n')
			
			preChr, prePos = preLine[0].split('M')
			newChr, newPos = newLine[0].split('M')
			
			if preChr == newChr:
				# If th couple in not excluded one
				OK = True
				if preLine[0] in EXCLCouple:
					if newLine[0] in EXCLCouple[preLine[0]]:
						OK = False
				if newLine[0] in EXCLCouple:
					if preLine[0] in EXCLCouple[newLine[0]]:
						OK = False
				
				if OK:
					
					############################
					# We can calculate recombination rate
				
					## First we need to identify impacted windows
					WindowsTofill = set()
					wBe = int(int(prePos)/WIN)
					wEn = int(int(newPos)/WIN)
					
					wToAdd = []
					while wBe <= wEn:
						wToAdd.append((wBe*WIN))
						wBe += 1
					
					# Identification of windows with markers
					DicoWin[preChr][wToAdd[0]+1][3] = True
					DicoWin[preChr][wToAdd[-1]+1][3] = True
					
					## We calculate the number of recombination between markers
					preRE = recode_matrix(preLine)
					newRE = recode_matrix(newLine)
					NumRec = numpy.nanmean(abs(preRE[1]-newRE[1]))
					if preRE[2] != newRE[2]:
						sys.exit('Oups, the program has crashed... There is marker of distinct type in the file and this is not allowed\n')
					
					## Filling the dictionary
					Value = NumRec/len(wToAdd)
					for n in wToAdd:
						DicoWin[preChr][n+1][2] += Value
		
		# Preparing result for drawing
		if options.chr == 'all':
			listChr = sorted(list(DicoWin.keys()))
		else:
			listChr = options.chr.split(',')
		# Collecting information in a first dictionary
		DicoValToPrint = {}
		mxVal = 0
		for n in listChr:
			listPos = sorted(list(DicoWin[n].keys()))
			DicoValToPrint[n] = [[],[]]
			for j in listPos:
				Value = (DicoWin[n][j][2]/(DicoWin[n][j][1]-DicoWin[n][j][0]+1))*WIN
				DicoValToPrint[n][0].append(Value)
				DicoValToPrint[n][1].append((DicoWin[n][j][1]+DicoWin[n][j][0])/2)
				if Value > mxVal:
					mxVal = Value
		
		# Smoothing data
		for n in listChr:
			DicoValToPrint[n][0].insert(0, DicoValToPrint[n][0][0])
			DicoValToPrint[n][1].insert(0, 1)
			DicoValToPrint[n][0].append(DicoValToPrint[n][0][-1])
			DicoValToPrint[n][1].append(chr_info[n])
			XNEW = numpy.linspace(DicoValToPrint[n][1][0],DicoValToPrint[n][1][-1],3000)
			POWER_SMOOTH = spline(DicoValToPrint[n][1],DicoValToPrint[n][0],XNEW)
			DicoValToPrint[n].append(XNEW)
			DicoValToPrint[n].append(POWER_SMOOTH)
		
		# Identification of regions with no data
		DicoRegionHigh = {}
		for n in listChr:
			listPos = sorted(list(DicoWin[n].keys()))
			start = False
			end = False
			for j in listPos:
				if DicoWin[n][j][3]: # The region contains markers
					if start: # The preceding region does not contain markers
						if not(n in DicoRegionHigh):
							DicoRegionHigh[n] = []
						DicoRegionHigh[n].append([start, end])
						start = False
						end = False
					else: # The preceding region contained markers
						pass
				else: # The region does not contain markers
					if start: # The preceding region does not contain markers
						end = DicoWin[n][j][1]
					else: # The preceding region contained markers
						start = DicoWin[n][j][0]
						end = DicoWin[n][j][1]
			if start:
				if not(n in DicoRegionHigh):
					DicoRegionHigh[n] = []
				DicoRegionHigh[n].append([start, end])
		
		# Printing results
		outfile = open(options.prefix+'.tab','w')
		if options.rmun == 'n':
			for n in listChr:
				for i in range(len(DicoValToPrint[n][2])):
					outfile.write('\t'.join([n, str(int(DicoValToPrint[n][2][i]-1)), str(int(DicoValToPrint[n][2][i])), str(max(0,DicoValToPrint[n][3][i]))]))
					outfile.write('\n')
		else: # Removing recombination rate estimated from surrounding windows
			for n in listChr:
				for i in range(len(DicoValToPrint[n][2])):
					NotComprised = True
					if n in DicoRegionHigh:
						for k in DicoRegionHigh[n]:
							value = int(DicoValToPrint[n][2][i])
							if k[0] <= value and value <= k[1]:
								NotComprised = False
					if NotComprised:
						outfile.write('\t'.join([n, str(int(DicoValToPrint[n][2][i]-1)), str(int(DicoValToPrint[n][2][i])), str(max(0,DicoValToPrint[n][3][i]))]))
						outfile.write('\n')
		outfile.close()
		
		outfile = open(options.prefix+'_high.tab','w')
		for n in listChr:
			if n in DicoRegionHigh:
				for k in DicoRegionHigh[n]:
					outfile.write('\t'.join(list(map(str, [n]+k))))
					outfile.write('\n')
		outfile.close()
		
		draw_chr(DicoValToPrint, chr_info, DicoRegionHigh, mxVal, options.prefix)

if __name__ == "__main__": __main__()
