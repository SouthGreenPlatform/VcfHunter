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
import os
import sys
import numpy
from scipy import stats
import math
import optparse
import datetime
import tempfile

def recode_matrix(LINE, RATIO):
	
	LISTE = LINE.split()
	# recuperation du nom de marker
	name = LISTE[0]
	coding = LISTE[1].split(',')
	if RATIO:
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
	else:
		if coding == ['nn','np']:
			LISTE_recode = recode_nn_np(LISTE[3:])
		elif coding == ['ll','lm']:
			LISTE_recode = recode_ll_lm(LISTE[3:])
		elif coding == ['hh','hk', 'kk']:
			LISTE_recode = recode_hh_hk_kk(LISTE[3:])
		elif coding == ['hh','k-']:
			LISTE_recode = recode_hh_k(LISTE[3:])
		else:
			sys.exit('Bad file format: The matrix file should contain 3 or 4 header columns\n'+LINE)
	return (name, LISTE_recode, coding)

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

def convert(MATRIX):
	file = open(MATRIX)
	# recording header
	header = file.readline().split()
	
	# recording matrix and reverse matrix
	marker, matrix, reverse_matrix, marker_coding = [], [], [], []
	for line in file:
		RECODE = recode_matrix(line, 'ratio' in header)
		marker.append(RECODE[0])
		matrix.append(RECODE[1])
		marker_coding.append(RECODE[2])
	file.close()
	
	# converting matrix into array
	return (marker, numpy.array(matrix), header, marker_coding)

def calculate_pairwise_recomb_rate(MATRIX, MATRIX2):
	"""
		Retrun a matrix with pairwise distance calculated
	"""
	NB_MARKER = MATRIX.shape[0]
	
	t0 = datetime.datetime.now()
	
	list2return = []
	
	define_to_stop = NB_MARKER/2
	
	NEW_MATRIX = numpy.array(MATRIX2)
	i = 0
	
	if len(MATRIX.shape) == 2:
		sub_list = numpy.nanmean(abs(MATRIX-NEW_MATRIX), axis=1)
		list2return.append(sub_list)
		i += 1
		sys.stdout.write(str(i)+" "+str((datetime.datetime.now()-t0)*(define_to_stop-i))+"\n")
		sys.stdout.flush()
		
		
		while i <= (NB_MARKER/2):
			t0 = datetime.datetime.now()
			i += 1
			NEW_MATRIX = numpy.concatenate((NEW_MATRIX[1:NB_MARKER,],NEW_MATRIX[0:1,]))
			sub_list = numpy.nanmean(abs(MATRIX-NEW_MATRIX), axis=1)
			list2return.append(sub_list)
			if i%100 == 0:
				sys.stdout.write(str(i)+" "+str((datetime.datetime.now()-t0)*(define_to_stop-i+1))+"\n")
				sys.stdout.flush()
			# print i, (datetime.datetime.now()-t0)*(define_to_stop-i+1)
	else:
		sub_list = numpy.nanmean(abs(MATRIX-NEW_MATRIX))
		list2return.append(sub_list)
	return numpy.array(list2return)

def invert_matrix(MATRIX, MARKER_CODING):
	"""
		Return the inverted matrix (0 converted to 1 and 1 converted to 0) and taking in account missing data
	"""
	
	list2return = []
	for n in range(len(MARKER_CODING)):
		list2return.append(invert_code(MARKER_CODING[n], MATRIX[n,:]))
	return numpy.array(list2return)

def reshape_and_print_full_matrix(MARKER, MATRIX, PREFIX):
	"""
		Reformat matrix to square format and print to output file
	"""
	
	matrix_size = MATRIX.shape
	final_matrix = []
	
	# Creating the empty matrix
	for n in range(matrix_size[1]):
		sub_list = [0]*matrix_size[1]
		final_matrix.append(sub_list)
	
	# Creating the matrix
	for n in range(matrix_size[0]):
		for i in range(matrix_size[1]):
			j = (n + i)% matrix_size[1]
			final_matrix[i][j] = MATRIX[n,i]
			final_matrix[j][i] = MATRIX[n,i]
	final_matrix = numpy.array(final_matrix)
	
	outfile = open(PREFIX+'_REC.tab','w')
	# Printing marker header
	outfile.write('\t'.join(['ID']+MARKER))
	outfile.write('\n')
	# Pinting line by line matrix
	for i in range(len(MARKER)):
		mot = [MARKER[i]]
		for k in range(len(MARKER)):
			mot.append(str(final_matrix[i,k]))
		outfile.write('\t'.join(mot))
		outfile.write('\n')
	outfile.close()

def invert_code(MARKER_CODING, CODE):
	INV_CODE = []
	for i in range(int(CODE.shape[0]/2)):
		if MARKER_CODING == ['nn','np'] or MARKER_CODING == ['ll','lm'] or MARKER_CODING == ['hh','k-']:
			INV_CODE.append(numpy.nan)
			if CODE[(i*2)+1] == 1.0:
				INV_CODE.append(0)
			elif CODE[(i*2)+1] == 0.0:
				INV_CODE.append(1)
			else:
				INV_CODE.append(numpy.nan)
		elif MARKER_CODING == ['hh','hk', 'kk']:
			if CODE[(i*2)] == 0.0 and CODE[(i*2)+1] == 1.0:
				INV_CODE.append(0)
				INV_CODE.append(1)
			else:
				if CODE[(i*2)] == 1.0:
					INV_CODE.append(0)
				elif CODE[(i*2)] == 0.0:
					INV_CODE.append(1)
				else:
					INV_CODE.append(numpy.nan)
				if CODE[(i*2)+1] == 1.0:
					INV_CODE.append(0)
				elif CODE[(i*2)+1] == 0.0:
					INV_CODE.append(1)
				else:
					INV_CODE.append(numpy.nan)
		else:
			sys.exit('The program existed without finishing: Unmanaged marker coding: '+','.join(MARKER_CODING))
	INV_CODE = numpy.array(INV_CODE)
	return INV_CODE

def printSegDist(MATRIX, OUT):
	
	outfile = open(OUT+'_SegDist.tab','w')
	
	file = open(MATRIX)
	header = file.readline().split()
	
	# Checking
	if not ('Marker' in header):
		sys.exit('No Marker column found in the matrix file')
	if not ('coding' in header):
		sys.exit('No coding column found in the matrix file')
	if not ('ratio' in header):
		sys.exit('No ratio column found in the matrix file')
	# Checking done
	
	for line in file:
		data = line.split()
		if data:
			mName = data[header.index('Marker')]
			coding = data[header.index('coding')].split(',')
			ratio = list(map(float, data[header.index('ratio')].split(',')))
			liste = []
			
			if len(coding) != len(ratio):
				sys.exit('Expected ration does not match marker coding!!!\n'+line)
			
			freq_obs = []
			for i in range(len(coding)):
				freq_obs.append(data[4:].count(coding[i]))
			total = numpy.sum(freq_obs)
			
			freq_theor = []
			for i in range(len(coding)):
				freq_theor.append(ratio[i]*total)
			
			chi2, p = stats.chisquare(freq_obs, f_exp=freq_theor)
			outfile.write('\t'.join([mName, str(-1*math.log(p, 10))]))
			outfile.write('\n')
	outfile.close()

def printSegDistProp(MATRIX, OUT):
	
	outfile = open(OUT+'_SegDist.tab','w')
	
	file = open(MATRIX)
	header = file.readline().split()
	
	# Checking
	if not ('Marker' in header) and  not ('marker' in header):
		sys.exit('No Marker column found in the matrix file')
	if not ('coding' in header):
		sys.exit('No coding column found in the matrix file')
	if not ('ratio' in header):
		sys.exit('No ratio column found in the matrix file')
	# Checking done
	
	for line in file:
		data = line.split()
		if data:
			if 'Marker' in header:
				mName = data[header.index('Marker')]
			else:
				mName = data[header.index('marker')]
			coding = data[header.index('coding')].split(',')
			ratio = list(map(float, data[header.index('ratio')].split(',')))
			liste = []
			
			if len(coding) != len(ratio):
				sys.exit('Expected ration does not match marker coding!!!\n'+line)
			
			freq_obs = []
			for i in range(len(coding)):
				freq_obs.append(data[4:].count(coding[i]))
			total = numpy.sum(freq_obs)
			
			freq_theor = []
			for i in range(len(coding)):
				freq_theor.append(ratio[i]*total)
			
			chi2, p = stats.chisquare(freq_obs, f_exp=freq_theor)
			outfile.write('\t'.join([mName, str((numpy.mean(abs(numpy.array(freq_obs) - numpy.array(freq_theor)))/sum(freq_obs)))]))
			outfile.write('\n')
	outfile.close()

def marker_decode(MARKER_CODING, RATIO, MARKER, CODE, REVERSE, INDI_NUMBER):
	
	if REVERSE:
		liste = [MARKER,','.join(MARKER_CODING), RATIO,'1']
		MARKER_CODING.reverse()
		for n in range(INDI_NUMBER):
			if str(CODE[n*2]) == 'nan' and str(CODE[n*2+1]) == 'nan':
				liste.append('--')
			elif str(CODE[n*2]) == 'nan':
				liste.append(MARKER_CODING[int(CODE[n*2+1])])
			else:
				liste.append(MARKER_CODING[int(CODE[n*2])+int(CODE[n*2+1])])
				
	else:
		liste = [MARKER, ','.join(MARKER_CODING), '0']
		for n in range(INDI_NUMBER):
			if str(CODE[n*2]) == 'nan' and str(CODE[n*2+1]) == 'nan':
				liste.append('--')
			elif str(CODE[n*2]) == 'nan':
				liste.append(MARKER_CODING[int(CODE[n*2+1])])
			else:
				liste.append(MARKER_CODING[int(CODE[n*2])+int(CODE[n*2+1])])
	return '\t'.join(liste)

def PhaseMarker(MATRIX, OUT):
	
	outfile = open(OUT+'_phased.tab','w')
	
	file = open(MATRIX)
	header = file.readline().split()
	if 'ratio' in header:
		indi_number = len(header[4:])
	else:
		indi_number = len(header[3:])
	outfile.write('\t'.join(header))
	outfile.write('\n')
	prec_coded_marker = ""
	for line in file:
		data = line.split()
		if data:
			if prec_coded_marker == "":
				coded = recode_matrix(line,'ratio' in header)
				prec_coded_marker = numpy.array(coded[1])
				prec_phase = 0
				outfile.write(line)
			else:
				coded = recode_matrix(line,'ratio' in header)
				coded_marker = numpy.array(coded[1])
				invert_coded_marker = invert_code(coded[2], coded_marker)
				ValuePhase0 = numpy.nanmean(abs(prec_coded_marker-coded_marker))
				ValuePhase1 = numpy.nanmean(abs(prec_coded_marker-invert_coded_marker))
				if ValuePhase0 <= ValuePhase1:
					if prec_phase == 0:
						outfile.write(line)
						prec_phase = 0
					else:
						outfile.write(marker_decode(coded[2], data[header.index('ratio')], data[0], coded_marker, True, indi_number))
						outfile.write('\n')
						prec_phase = 1
				else:
					if prec_phase == 0:
						outfile.write(marker_decode(coded[2], data[header.index('ratio')], data[0], coded_marker, True, indi_number))
						outfile.write('\n')
						prec_phase = 1
					else:
						outfile.write(line)
						prec_phase = 0
				prec_coded_marker = coded_marker
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n This program calculates pairwise recombination from a marker matrix file.")
	# Wrapper options. 
	parser.add_option( '-m', '--matrix', dest='matrix', default=None, help='The matrix marker file')
	parser.add_option( '-o', '--output', dest='output', default=None, help='Output prefix name')
	parser.add_option( '-p', '--phased', dest='phased', default='n', help='The matrix has been phased (y or n), [default: %default]')
	parser.add_option( '-s', '--steps', dest='steps', default=None, help='Analysis to perform\t\t\t\t\t\t'
	'R: Calculate recombination rate\t\t\t\t\t'
	'S: Calculate segregation distortions (-log10(Chi-square test))\t\t\t\t'
	's: Calculate segregation distortions (as proportion of individuals deviating from expectation)\t\t\t\t'
	'P: Phase markers based on their order\t\t\t\t')
	
	(options, args) = parser.parse_args()
	
	if options.matrix == None:
		sys.exit('Please provide a reference multifasta file to --matrix argument')
	if options.output == None:
		sys.exit('Please provide a list of  file to --output argument')
	if options.steps == None:
		sys.exit('Please provide a list of  file to --steps argument')
	V = os.getpid()
	
	t0 = datetime.datetime.now()
	
	if 'R' in options.steps:
		# converting data into array
		datas = convert(options.matrix)
		marker = datas[0]
		matrix = datas[1]
		header = datas[2]
		marker_coding = datas[3]
		datas = None
		sys.stdout.write("Total time : "+str(datetime.datetime.now()-t0)+"\n")
		sys.stdout.flush()
		toto = os.system("pmap %s | grep total" % V)
	
	if 'R' in options.steps:
	
		# calculate pairwise recombination events
		pairwise_recomb = calculate_pairwise_recomb_rate(matrix, matrix)
		sys.stdout.write("Pairwise calculation done\nTotal time until now: "+str(datetime.datetime.now()-t0)+"\n")
		sys.stdout.flush()
		toto = os.system("pmap %s | grep total" % V)
		
		# rephasing data if necessary
		if options.phased == 'n':
			# inverting matrix
			inverted_matrix = invert_matrix(matrix, marker_coding)
			# calculating recombination rate on the inverted matrix
			pairwise_recomb_inv = calculate_pairwise_recomb_rate(matrix, inverted_matrix)
			# returning minimal recombination rate value
			pairwise_recomb = numpy.minimum(pairwise_recomb, pairwise_recomb_inv)
			sys.stdout.write("Rephasing calculation done\nTotal time until now: "+str(datetime.datetime.now()-t0)+"\n")
			sys.stdout.flush()
			toto = os.system("pmap %s | grep total" % V)
		
		# printing recombination matrix
		reshape_and_print_full_matrix(marker, pairwise_recomb, options.output)
		sys.stdout.write("Printing recombination done \nTotal time until now: "+str(datetime.datetime.now()-t0)+"\n")
		sys.stdout.flush()
		toto = os.system("pmap %s | grep total" % V)
	
	if 'S' in options.steps:
		printSegDist(options.matrix, options.output)
	
	if 's' in options.steps:
		printSegDistProp(options.matrix, options.output)
	
	if 'P' in options.steps:
		PhaseMarker(options.matrix, options.output)
		
if __name__ == "__main__": __main__()