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
import optparse
import sys
import time
import math
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def filter_on_read_cov(DATA, FORMAT, MINCOV, MAXCOV, MINALCOV, WINFREQ, GLOB_allele_ratio, GLOB_coverage, TO_EXCLUDE, VCF_HEADER, DICO_ACC_TO_TREAT):
	"""
		Convert VCF calling to a list of recoded variant

		param DATA: list of coded variant calling
		type DATA: list
		param FORMAT: list of vcf tags
		type FORMAT: list
		param MINCOV: minimal site coverage per accessions
		type MINCOV: int
		param MAXCOV: maximal site coverage per accessions
		type MAXCOV: int
		param MINALCOV: Minimum allele coverage threshold
		type MINALCOV: int
		param WINFREQ: ]0.01:0.1[
		type WINFREQ: str
		param GLOB_allele_ratio: Global variables directely filled
		type GLOB_allele_ratio: list
		param GLOB_coverage: Global variables directely filled
		type GLOB_coverage: list
		param TO_EXCLUDE: Accessions to exclude
		type TO_EXCLUDE: set
		param VCF_HEADER: VCF header
		type VCF_HEADER: list
		param DICO_ACC_TO_TREAT: list of treated accession (excluding not_used and to_exclude)
		type DICO_ACC_TO_TREAT: list
		return: list
	"""
	
	no_cov, no_freq, no_read, no_cov_plus, troncated_format, too_cov, more_than_di_allele = 0, 0, 0, 0, 0, 0, 0
	all_missing, monomorphous, dimorphous_sites, more_than_dimorphous = 0, 0, 0, 0
	
	if WINFREQ:
		winfreq = list(map(float, WINFREQ.split(':')))
	
	genotype_set = set()
	
	global_cov = []
	global_allele_ratio = []
	
	# Cheking if Flags are present in the FORMAT tag
	if 'AD' in FORMAT and 'DP' in FORMAT and 'GT' in FORMAT:
		liste2return = []
		liste2return_intermediate = []
		setnot2return = set()
		for ipos in range(len(DATA)):
			acc = DATA[ipos]
			data = acc.split(':')
			acc_id = VCF_HEADER[ipos]
			##### Looking for "real variant allele"
			if not(acc_id in TO_EXCLUDE):
				# Cheking if the accession is in the same format of the FORMAT tag
				if len(FORMAT) == len(data):
					# Cheking if there is information in the DP and AD tag
					if '.' in data[FORMAT.index('DP')] or '.' in data[FORMAT.index('AD')]:
						no_read += 1
					elif data[FORMAT.index('DP')] == '0':
						no_read += 1
					else:
						# additional step to identify between false and true multiallelic variant sites
						allele_cov_info = list(map(int, data[FORMAT.index('AD')].split(',')))
						# Checking for cases where the coverage is due to non A,T,G,C bases
						if sum(allele_cov_info) > 0:
							for n in range(len(allele_cov_info)):
								if allele_cov_info[n]/float(sum(allele_cov_info)) >= winfreq[0] and allele_cov_info[n] >= MINALCOV:
									setnot2return.add(n)
			##### End of looking for "real variant allele"
			
		# print (setnot2return)
			
		for ipos in range(len(DATA)):
			acc = DATA[ipos]
			data = acc.split(':')
			code = './.'
			acc_id = VCF_HEADER[ipos]
			# Cheking if accession needs to be treated
			if not(acc_id in TO_EXCLUDE):
				# Cheking if the accession is in the same format of the FORMAT tag
				if len(FORMAT) == len(data):
					# Cheking if there is information in the DP and AD tag
					if '.' in data[FORMAT.index('DP')] or '.' in data[FORMAT.index('AD')]:
						no_read += 1
					elif data[FORMAT.index('DP')] == '0':
						no_read += 1
					else:
						# additional step to identify between false and true multiallelic variant sites
						allele_cov_info = list(map(int, data[FORMAT.index('AD')].split(',')))
						liste_allele = []
						# Checking for cases where the coverage is due to non A,T,G,C bases
						if sum(allele_cov_info) > 0:
							for n in range(len(allele_cov_info)):
								if n in setnot2return:
									if allele_cov_info[n] > 0:
										liste_allele.append(str(n))
						# The genotype is mono allelique
						if len(liste_allele) == 1:
							# Cheking if the variant is too much covered
							coverage_value = [allele_cov_info[int(liste_allele[0])]]
							DP_tag = sum(coverage_value)
							if DP_tag > 0:
								if DP_tag <= MAXCOV:
									if acc_id in DICO_ACC_TO_TREAT:
										global_cov.append(DP_tag)
							if DP_tag > MAXCOV:
								too_cov += 1
							
							# Cheking if the variant is sufficiantly covered
							elif DP_tag >= MINCOV:
								if acc_id in DICO_ACC_TO_TREAT:
									global_allele_ratio.append(0)
								code = '/'.join(liste_allele*2)
							# The variant site is not sufficiently covered
							elif DP_tag != 0:
								no_cov += 1
							else:
								no_read += 1
								
						# The genotype is di-allelique
						elif len(liste_allele) == 2:
							# Cheking if the variant is too much covered
							coverage_value = []
							for n in liste_allele:
								coverage_value.append(allele_cov_info[int(n)])
							DP_tag = sum(coverage_value)
							if DP_tag > 0:
								if DP_tag <= MAXCOV:
									if acc_id in DICO_ACC_TO_TREAT:
										global_cov.append(DP_tag)
							
							if DP_tag > MAXCOV:
								too_cov += 1
							
							# Cheking if the variant is sufficiantly covered
							elif DP_tag >= MINCOV:
								min_freq = min(coverage_value)/float(sum(coverage_value))
								if acc_id in DICO_ACC_TO_TREAT:
									global_allele_ratio.append(min_freq)
								# Coding heterozygous based on winfreq
								if min_freq >= winfreq[1] and min(coverage_value) >= MINALCOV:
									code = '/'.join(liste_allele)
								elif min_freq < winfreq[0]:
									code = '/'.join([liste_allele[coverage_value.index(max(coverage_value))]]*2)
								# Coding ambiguous based on winfreq and or MINALCOV
								else:
									no_freq += 1
									
							# The variant site is not sufficiently covered
							elif DP_tag != 0:
								no_cov += 1
							else:
								no_read += 1
						else:
							more_than_di_allele += 1
				# Bad format for the accession calling (equivalent to no reads)
				else:
					no_read += 1
			# print acc_id, allele_cov_info, liste_allele, code, no_read, no_cov, too_cov, no_freq, no_cov_plus, troncated_format, more_than_di_allele, all_missing, monomorphous, dimorphous_sites, more_than_dimorphous
			liste2return_intermediate.append(code)
		
		# True di-allelique identification
		if len(setnot2return) == 0:
			all_missing = 1
			# print 'all_missing'
			liste2return = ['./.']*len(DATA)
		elif len(setnot2return) == 1:
			monomorphous = 1
			# print 'monomorphous'
			liste2return = ['./.']*len(DATA)
		elif len(setnot2return) == 2:
			dimorphous_sites = 1
			# print 'dimorphous_sites'
			# recoding variants
			listenot2return = list(map(str, list(setnot2return)))
			listenot2return.sort()
			for n in liste2return_intermediate:
				code = []
				for allele in n.split('/'):
					if allele == '.':
						code.append(allele)
					else:
						code.append(str(listenot2return.index(allele)))
				liste2return.append('/'.join(code))
			for n in global_cov:
				GLOB_coverage.append(n)
			for n in global_allele_ratio:
				if n != 0:
					GLOB_allele_ratio.append(n)
		else:
			more_than_dimorphous = 1
			# print 'more_than_dimorphous'
			liste2return = ['./.']*len(DATA)
	# Bad format for the FORMAT tag
	else:
		liste2return = ['./.']*len(DATA)
		troncated_format = 1
	
	# print liste2return_intermediate
	# print liste2return
	# print liste2return_intermediate == liste2return
	
	return [liste2return, no_read, no_cov, too_cov, no_freq, no_cov_plus, troncated_format, more_than_di_allele, all_missing, monomorphous, dimorphous_sites, more_than_dimorphous, liste2return_intermediate]

def get_missing(LISTE, ACCESSION_HEADER, DICO_ACC):
	"""
		Calculate the missing data proportion for the filter

		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param DICO_ACC: set of accession to work with
		type DICO_ACC: set()
		return: float
	"""
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	# calculating missing
	missing = 0
	for n in range(len(LISTE)):
		if ACCESSION_HEADER[n] in DICO_ACC:
			if LISTE[n] == './.':
				missing += 1
	return missing/float(len(DICO_ACC))

def KhiSquare(LISTE, ACCESSION_HEADER, DICO_ACC, POP_TYPE):
	"""
		Calculate KhiSquare P-Value

		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param DICO_ACC: set of accession to work with
		type DICO_ACC: set()
		param POP_TYPE: Population type
		type POP_TYPE: str
		return: list
	"""
	
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	
	# Counting genotype frequencies
	homo1, hetero, homo2, total = 0, 0, 0, 0
	for n in range(len(LISTE)):
		if ACCESSION_HEADER[n] in DICO_ACC:
			if LISTE[n] == '0/0':
				homo1 += 1
				total += 1
			elif LISTE[n] == '0/1':
				hetero += 1
				total += 1
			elif LISTE[n] == '1/1':
				homo2 += 1
				total += 1
			elif LISTE[n] != './.':
				sys.exit('Oups, there is a bug in the program2! KhiSquare '+LISTE[n])
	
	
	freq_obs_25_50_25 = np.array([homo1, hetero, homo2])
	freq_obs1_50_50 = np.array([homo1, hetero])
	freq_obs2_50_50 = np.array([homo2, hetero])
	freq_obs1_25_75 = np.array([homo1, hetero])
	freq_obs2_25_75 = np.array([homo2, hetero])
	
	
	freq_theor_25_50_25 = np.array([0.25*total, 0.5*total, 0.25*total])
	freq_theor_50_50 = np.array([0.5*total, 0.5*total])
	freq_theor_25_75 = np.array([0.25*total, 0.75*total])
	
	
	
	# Checking the population type and calculating p-value
	liste2return = []
	max_pvalue = 0
	if POP_TYPE == 'SELF':
		chi2, p = stats.chisquare(freq_obs_25_50_25, f_exp=freq_theor_25_50_25)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['hh,hk,kk']
	
	elif POP_TYPE == 'SELFPOL':
		chi2, p = stats.chisquare(freq_obs1_25_75, f_exp=freq_theor_25_75)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['hh,k-']
		chi2, p = stats.chisquare(freq_obs2_25_75, f_exp=freq_theor_25_75)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['hh,k-']
			
	elif POP_TYPE == 'BiP':
		chi2, p = stats.chisquare(freq_obs_25_50_25, f_exp=freq_theor_25_50_25)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['hh,hk,kk']
			
		chi2, p = stats.chisquare(freq_obs1_50_50, f_exp=freq_theor_50_50)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['nn,np']
			
		chi2, p = stats.chisquare(freq_obs2_50_50, f_exp=freq_theor_50_50)
		if p >= max_pvalue:
			max_pvalue = p
			min_chi2 = chi2
			coding = ['nn,np']
	
	return [max_pvalue, min_chi2, coding]

def decode(ACCESSION_HEADER, LISTE, ALLELES, TO_EXCLUDE, RAW_DATA, FORMAT, ADDCOV):
	"""
		Recode to give genotype
		
		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param TO_EXCLUDE: set of accession to work with
		type TO_EXCLUDE: set()
		param ALLELES: Liste of allele
		type ALLELES: list
		param RAW_DATA: List of raw data
		type RAW_DATA: list
		param FORMAT: list of vcf tags
		type FORMAT: list
		return: list
	"""
	
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	
	if ADDCOV == 'y':
		liste2return = []
		for n in range(len(LISTE)):
			if not(ACCESSION_HEADER[n] in TO_EXCLUDE):
				mot = []
				for geno in LISTE[n].split('/'):
					if geno == '.':
						mot.append(geno)
					else:
						mot.append(ALLELES[int(geno)])
				data = RAW_DATA[n].split(':')
				if len(data) == len(FORMAT):
					liste2return.append('/'.join(mot)+':'+data[FORMAT.index('AD')])
				else:
					liste2return.append('/'.join(mot)+':0,0')
		return liste2return
	
	else:
		liste2return = []
		for n in range(len(LISTE)):
			if not(ACCESSION_HEADER[n] in TO_EXCLUDE):
				mot = []
				for geno in LISTE[n].split('/'):
					if geno == '.':
						mot.append(geno)
					else:
						mot.append(ALLELES[int(geno)])
				data = RAW_DATA[n].split(':')
				liste2return.append('/'.join(mot))
		return liste2return

def ident_parent(ACCESSION_HEADER, LISTE, PARENT, KHISTAT, POP_TYPE):
	"""
		Identify parent marker
		
		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param PARENT: List of parent names
		type PARENT: list
		param KHISTAT: Allele coding
		type KHISTAT: list
		param POP_TYPE: Population type
		type POP_TYPE: str
		return: list with [parent marker, first or second]
	"""
	
	
	dico_parent = {}
	dico_parent['M'] = []
	dico_parent['O'] = []
	dico_parent['E'] = []
	for n in PARENT:
		genotype = LISTE[ACCESSION_HEADER.index(n)]
		geno = set(genotype.split('/'))
		if len(geno) == 1:
			if '.' in geno:
				dico_parent['M'].append(n)
			else:
				dico_parent['O'].append(n)
		else:
			dico_parent['E'].append(n)
	
	if POP_TYPE == "BiP":
		if len(PARENT) == 2:
			# observed segregation is 50/50
			if KHISTAT == ['nn,np']:
				# only one parent heterozygous (the other is missing or homozygous)
				if len(dico_parent['E']) == 1:
					return [dico_parent['E'][0], PARENT.index(dico_parent['E'][0])]
				# only one parent homozygous (the other is missing)
				elif len(dico_parent['O']) == 1:
					return [dico_parent['M'][0], PARENT.index(dico_parent['M'][0])]
				#both parent are heterozygous OR both parent are homozygous OR both parents are missing
				else:
					return ['unknown', 'x']
			
			# observed segregation is 25/50/25
			elif KHISTAT == ['hh,hk,kk']:
				return ['Bridge', 'x']
		elif len(PARENT) == 0:
			# observed segregation is 50/50
			if KHISTAT == ['nn,np']:
				return ['unknown', 'x']
			
			# observed segregation is 25/50/25
			elif KHISTAT == ['hh,hk,kk']:
				return ['Bridge', 'x']
		else:
			sys.exit('Unmanaged parent number')
	elif POP_TYPE == "SELF":
		if len(PARENT) == 0:
			return ['unknown', 'x']
		elif len(PARENT) == 1:
			# Provided parent is heterozygous
			if len(dico_parent['E']) == 1:
				return [dico_parent['E'][0], PARENT.index(dico_parent['E'][0])]
			# Provided parent is missing or homozygous
			elif len(dico_parent['M']) == 1:
				return ['unknown', 'x']
			else:
				return ['unknown', 'x']	
		elif len(PARENT) == 2:
			# One parent is homozygous and the other is heterozygous
			if len(dico_parent['E']) == 1 and len(dico_parent['O']) == 1:
				return [dico_parent['E'][0], PARENT.index(dico_parent['E'][0])]
			else:
				return ['unknown', 'x']	
		else:
			sys.exit('Unmanaged parent number')
	elif POP_TYPE == "SELFPOL":
		if len(PARENT) == 0:
			return ['unknown', 'x']
		elif len(PARENT) == 1:
			# Provided parent is heterozygous
			if len(dico_parent['E']) == 1:
				return [dico_parent['E'][0], PARENT.index(dico_parent['E'][0])]
			# Provided parent is missing or homozygous
			elif len(dico_parent['M']) == 1:
				return ['unknown', 'x']
			else:
				return ['unknown', 'x']	
		elif len(PARENT) == 2:
			# One parent is homozygous and the other is heterozygous
			if len(dico_parent['E']) == 1 and len(dico_parent['O']) == 1:
				return [dico_parent['E'][0], PARENT.index(dico_parent['E'][0])]
			else:
				return ['unknown', 'x']	
		else:
			sys.exit('Unmanaged parent number')
	else:
		sys.exit('Unmanaged population type')

def recode2onemap(ACCESSION_HEADER, LISTE, KHISTAT, POP_TYPE, DICO_ACC, PARENT_STATUS, MARKER_NAME):
	"""
		recode marker for onemap
		
		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param DICO_ACC: set of accession to work with
		type DICO_ACC: set()
		param KHISTAT: Allele coding
		type KHISTAT: str
		param POP_TYPE: Population type
		type POP_TYPE: str
		param KHISTAT: Allele coding
		type KHISTAT: list
		param PARENT_STATUS: Parent status (0 for parent 1, 1 for parent 2, 'x' for unknown parent)
		type PARENT_STATUS: value
		param MARKER_NAME: Matker name
		type MARKER_NAME: str
		return: str
	"""
	
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	
	
	# Counting genotype frequencies
	homo1, hetero, homo2, total = 0, 0, 0, 0
	for n in range(len(LISTE)):
		if ACCESSION_HEADER[n] in DICO_ACC:
			if LISTE[n] == '0/0':
				homo1 += 1
				total += 1
			elif LISTE[n] == '0/1':
				hetero += 1
				total += 1
			elif LISTE[n] == '1/1':
				homo2 += 1
				total += 1
			elif LISTE[n] != './.':
				sys.exit('Oups, there is a bug in the program2! recode2onemap '+LISTE[n])
	
	if homo2 < homo1:
		maxi = 'homo1'
	else:
		maxi = 'homo2'
	
	# Recoding
	liste2return = []
	
	if POP_TYPE == "BiP" or POP_TYPE == "SELF" or POP_TYPE == "SELFPOL":
		if KHISTAT == ['nn,np']:		
			if PARENT_STATUS == 0:
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('a')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('a')
						elif LISTE[n] == '0/1':
							liste2return.append('ab')
						elif LISTE[n] == '0/0':
							liste2return.append('-')
						elif LISTE[n] == '1/1':
							liste2return.append('-')
						elif LISTE[n] == './.':
							liste2return.append('-')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
				return '*'+MARKER_NAME + ' ' + 'D1.10' + '\t' + ','.join(liste2return)
			elif PARENT_STATUS == 1:
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('a')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('a')
						elif LISTE[n] == '0/1':
							liste2return.append('ab')
						elif LISTE[n] == '0/0':
							liste2return.append('-')
						elif LISTE[n] == '1/1':
							liste2return.append('-')
						elif LISTE[n] == './.':
							liste2return.append('-')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
				return '*'+MARKER_NAME + ' ' + 'D2.15' + '\t' + ','.join(liste2return)
			elif PARENT_STATUS == 'x':
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('a')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('a')
						elif LISTE[n] == '0/1':
							liste2return.append('ab')
						elif LISTE[n] == '0/0':
							liste2return.append('-')
						elif LISTE[n] == '1/1':
							liste2return.append('-')
						elif LISTE[n] == './.':
							liste2return.append('-')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
				return '*'+MARKER_NAME + ' ' + 'D1.10' + '\t' + ','.join(liste2return)
			else:
				sys.exit('Bug in Parent status'+PARENT_STATUS)
		elif KHISTAT == ['hh,hk,kk']:
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0':
						liste2return.append('a')
					elif LISTE[n] == '1/1':
						liste2return.append('b')
					elif LISTE[n] == '0/1':
						liste2return.append('ab')
					elif LISTE[n] == './.':
						liste2return.append('-')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
			return '*'+MARKER_NAME + ' ' + 'B3.7' + '\t' + ','.join(liste2return)
		elif KHISTAT == ['hh,k-']:
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0' and maxi == 'homo1':
						liste2return.append('o')
					elif LISTE[n] == '1/1' and maxi == 'homo2':
						liste2return.append('o')
					elif LISTE[n] == '0/1':
						liste2return.append('a')
					elif LISTE[n] == '0/0':
						liste2return.append('-')
					elif LISTE[n] == '1/1':
						liste2return.append('-')
					elif LISTE[n] == './.':
						liste2return.append('-')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
			return '*'+MARKER_NAME + ' ' + 'C.8' + '\t' + ','.join(liste2return)
		else:
			sys.exit('bug')

def recode2JoinMap(ACCESSION_HEADER, LISTE, KHISTAT, POP_TYPE, DICO_ACC, PARENT_STATUS, MARKER_NAME, COUNT):
	"""
		recode marker for Joinmap
		
		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param DICO_ACC: set of accession to work with
		type DICO_ACC: set()
		param KHISTAT: Allele coding
		type KHISTAT: str
		param POP_TYPE: Population type
		type POP_TYPE: str
		param KHISTAT: Allele coding
		type KHISTAT: list
		param PARENT_STATUS: Parent status (0 for parent 1, 1 for parent 2, 'x' for unknown parent)
		type PARENT_STATUS: value
		param MARKER_NAME: Marker name
		type MARKER_NAME: str
		param COUNT: Marker count
		type COUNT: str
		return: str
	"""
	
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	
	
	# Counting genotype frequencies
	homo1, hetero, homo2, total = 0, 0, 0, 0
	for n in range(len(LISTE)):
		if ACCESSION_HEADER[n] in DICO_ACC:
			if LISTE[n] == '0/0':
				homo1 += 1
				total += 1
			elif LISTE[n] == '0/1':
				hetero += 1
				total += 1
			elif LISTE[n] == '1/1':
				homo2 += 1
				total += 1
			elif LISTE[n] != './.':
				sys.exit('Oups, there is a bug in the program2! recode2JoinMap '+LISTE[n])
	
	if homo2 < homo1:
		maxi = 'homo1'
	else:
		maxi = 'homo2'
	
	# Recoding
	liste2return = []
	
	if POP_TYPE == "BiP" or POP_TYPE == "SELF" or POP_TYPE == "SELFPOL":
		if KHISTAT == ['nn,np']:		
			if PARENT_STATUS == 0:
				liste2return.append(MARKER_NAME+'\t<nnxnp>\t\t(nn,np) ; '+COUNT+'\n')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('nn')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('nn')
						elif LISTE[n] == '0/1':
							liste2return.append('np')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			elif PARENT_STATUS == 1:
				liste2return.append(MARKER_NAME+'\t<lmxll>\t\t(ll,lm) ; '+COUNT+'\n')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('ll')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('ll')
						elif LISTE[n] == '0/1':
							liste2return.append('lm')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			elif PARENT_STATUS == 'x':
				liste2return.append(MARKER_NAME+'\t<nnxnp>\t\t(nn,np) ; '+COUNT+'\n')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('nn')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('nn')
						elif LISTE[n] == '0/1':
							liste2return.append('np')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			else:
				sys.exit('Bug in Parent status'+PARENT_STATUS)
		elif KHISTAT == ['hh,hk,kk']:
			liste2return.append(MARKER_NAME+'\t<hkxhk>\t\t(hh,hk,kk) ; '+COUNT+'\n')
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0':
						liste2return.append('hh')
					elif LISTE[n] == '1/1':
						liste2return.append('kk')
					elif LISTE[n] == '0/1':
						liste2return.append('hk')
					elif LISTE[n] == './.':
						liste2return.append('--')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
		elif KHISTAT == ['hh,k-']:
			liste2return.append(MARKER_NAME+'\t<hkxhk>\t\t(hh,k-) ; '+COUNT+'\n')
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0' and maxi == 'homo1':
						liste2return.append('hh')
					elif LISTE[n] == '1/1' and maxi == 'homo2':
						liste2return.append('hh')
					elif LISTE[n] == '0/1':
						liste2return.append('k-')
					elif LISTE[n] == '0/0':
						liste2return.append('--')
					elif LISTE[n] == '1/1':
						liste2return.append('--')
					elif LISTE[n] == './.':
						liste2return.append('--')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
		else:
			sys.exit('bug')
		return ' '.join(liste2return)+'\n'

def recode2tab(ACCESSION_HEADER, LISTE, KHISTAT, POP_TYPE, DICO_ACC, PARENT_STATUS, MARKER_NAME):
	"""
		recode marker for Joinmap
		
		param LISTE: list of recoded variant calling
		type LISTE: list
		param ACCESSION_HEADER: list corresponding to accession header
		type ACCESSION_HEADER: list
		param DICO_ACC: set of accession to work with
		type DICO_ACC: set()
		param KHISTAT: Allele coding
		type KHISTAT: str
		param POP_TYPE: Population type
		type POP_TYPE: str
		param KHISTAT: Allele coding
		type KHISTAT: list
		param PARENT_STATUS: Parent status (0 for parent 1, 1 for parent 2, 'x' for unknown parent)
		type PARENT_STATUS: value
		param MARKER_NAME: Marker name
		type MARKER_NAME: str
		return: str
	"""
	
	# just in case
	if len(LISTE) != len(ACCESSION_HEADER):
		sys.exit('Oups, there is a bug in the program1!')
	
	
	# Counting genotype frequencies
	homo1, hetero, homo2, total = 0, 0, 0, 0
	for n in range(len(LISTE)):
		if ACCESSION_HEADER[n] in DICO_ACC:
			if LISTE[n] == '0/0':
				homo1 += 1
				total += 1
			elif LISTE[n] == '0/1':
				hetero += 1
				total += 1
			elif LISTE[n] == '1/1':
				homo2 += 1
				total += 1
			elif LISTE[n] != './.':
				sys.exit('Oups, there is a bug in the program2! recode2tab '+LISTE[n])
	
	if homo2 < homo1:
		maxi = 'homo1'
	else:
		maxi = 'homo2'
	
	# Recoding
	liste2return = []
	
	if POP_TYPE == "BiP" or POP_TYPE == "SELF" or POP_TYPE == "SELFPOL":
		if KHISTAT == ['nn,np']:		
			if PARENT_STATUS == 0:
				liste2return.append(MARKER_NAME)
				liste2return.append('nn,np')
				liste2return.append('0')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('nn')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('nn')
						elif LISTE[n] == '0/1':
							liste2return.append('np')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			elif PARENT_STATUS == 1:
				liste2return.append(MARKER_NAME)
				liste2return.append('ll,lm')
				liste2return.append('0')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('ll')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('ll')
						elif LISTE[n] == '0/1':
							liste2return.append('lm')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			elif PARENT_STATUS == 'x':
				liste2return.append(MARKER_NAME)
				liste2return.append('nn,np')
				liste2return.append('0')
				for n in range(len(LISTE)):
					if ACCESSION_HEADER[n] in DICO_ACC:
						if LISTE[n] == '0/0' and maxi == 'homo1':
							liste2return.append('nn')
						elif LISTE[n] == '1/1' and maxi == 'homo2':
							liste2return.append('nn')
						elif LISTE[n] == '0/1':
							liste2return.append('np')
						elif LISTE[n] == '0/0':
							liste2return.append('--')
						elif LISTE[n] == '1/1':
							liste2return.append('--')
						elif LISTE[n] == './.':
							liste2return.append('--')
						else:
							sys.exit('Probleme in allele coding'+LISTE[n])
			else:
				sys.exit('Bug in Parent status'+PARENT_STATUS)
		elif KHISTAT == ['hh,hk,kk']:
			liste2return.append(MARKER_NAME)
			liste2return.append('hh,hk,kk')
			liste2return.append('0')
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0':
						liste2return.append('hh')
					elif LISTE[n] == '1/1':
						liste2return.append('kk')
					elif LISTE[n] == '0/1':
						liste2return.append('hk')
					elif LISTE[n] == './.':
						liste2return.append('--')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
		elif KHISTAT == ['hh,k-']:
			liste2return.append(MARKER_NAME)
			liste2return.append('hh,k-')
			liste2return.append('0')
			for n in range(len(LISTE)):
				if ACCESSION_HEADER[n] in DICO_ACC:
					if LISTE[n] == '0/0' and maxi == 'homo1':
						liste2return.append('hh')
					elif LISTE[n] == '1/1' and maxi == 'homo2':
						liste2return.append('hh')
					elif LISTE[n] == '0/1':
						liste2return.append('k-')
					elif LISTE[n] == '0/0':
						liste2return.append('--')
					elif LISTE[n] == '1/1':
						liste2return.append('--')
					elif LISTE[n] == './.':
						liste2return.append('--')
					else:
						sys.exit('Probleme in allele coding'+LISTE[n])
		else:
			sys.exit('bug')
		return '\t'.join(liste2return)

def draw_statistics_plot(OUT, ALLELE_RATIO, COVERAGE_RATIO, MISSING, GLOB_KHI2, WINFREQ, PVALUE, MAX_MISSING, MINCOV):
	"""
		Draw plot
		
		param OUT: output file name
		type OUT: list
		param ALLELE_RATIO: numpy array
		type ALLELE_RATIO: numpy array
		param COVERAGE_RATIO: numpy array
		type COVERAGE_RATIO: numpy array
		param MISSING: numpy array
		type MISSING: numpy array
		param GLOB_KHI2: numpy array
		type GLOB_KHI2: numpy array
		param WINFREQ: list
		type WINFREQ: list
		param PVALUE: float
		type PVALUE: float
		param MAX_MISSING: float
		type MAX_MISSING: float
		param MINCOV: int
		type MINCOV: int
		return: void
	"""
	
	fig = plt.figure(1)
	ax = plt.subplot(221)
	n, bins, patches = ax.hist(COVERAGE_RATIO, 100, normed=1, facecolor='b', linewidth=0.001)
	ax.axvline(x=MINCOV, color='r')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_title('Site coverage (excluding 0)', fontweight='bold', fontsize=10)
	
	ax = plt.subplot(222)
	n, bins, patches = ax.hist(ALLELE_RATIO, 50, normed=1, facecolor='b', linewidth=0.001)
	if WINFREQ:
		winfreq = list(map(float, WINFREQ.split(':')))
		ax.axvline(x=winfreq[0], color='r')
		ax.axvline(x=winfreq[1], color='r')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_title('Minimal allele ratio (excluding 0)', fontweight='bold', fontsize=10)
	
	ax = plt.subplot(223)
	n, bins, patches = ax.hist(MISSING, 100, normed=1, facecolor='b', linewidth=0.001)
	ax.axvline(x=MAX_MISSING, color='r')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_title('Missing data proportion', fontweight='bold', fontsize=10)
	
	ax = plt.subplot(224)
	n, bins, patches = ax.hist(GLOB_KHI2, 100, normed=1, facecolor='b', linewidth=0.001)
	ax.axvline(x=math.log(PVALUE, 10), color='r')
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)
	ax.set_title('log(P-value) distribution', fontweight='bold', fontsize=10)
			
	fig.savefig(OUT)
	plt.close(fig)

def get_tags(SEQ_DIC, CHR, POS, OUT, MNAME):
	
	start = POS - 126
	end = POS + 125
	
	if start < 0:
		start = 0
	
	outmark = open(OUT, 'a')
	outmark.write('>'+MNAME+'\n')
	outmark.write(str(SEQ_DIC[CHR].seq)[start:end]+'\n')
	outmark.close()
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(
		description="This script filter VCF file to create matrix which can be used for genetical mapping studies",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr), Franc-Christophe BAURENS (franc-christophe.baurens@cirad.fr) "
		"and Olivier GARSMEUR (olivier.garsmeur@cirad.fr)")
	# Wrapper options.
	parser.add_option( '-v', '--vcf',		dest='vcf',			default=None,		help='The VCF file')
	parser.add_option( '-R', '--ref',		dest='ref',			default=None,		help='The Reference fasta file. [Default: %default]')
	parser.add_option( '-P', '--parent',	dest='parent',		default=None,		help='Parent of the population names (separated by ":").')
	parser.add_option( '-e', '--exclude',	dest='exclude',		default=None,		help='Accession to exclude from the final file.')
	parser.add_option( '-n', '--NoUsed',	dest='NoUsed',		default=None,		help='Accession to exclude from the filtering.')
	parser.add_option( '-m', '--MinCov',	dest='MinCov',		default='10',		help='Minimal read coverage for site. [Default: %default]')
	parser.add_option( '-M', '--MaxCov',	dest='MaxCov',		default='1000',		help='Maximal read coverage for site. [Default: %default]')
	parser.add_option( '-f', '--WinFreq',	dest='WinFreq',		default='0.05:0.1',	help='Window for minority allele coverage frequency to be insufficient to call a heterozygous but to high to call an homozygous (example: "0.05:0.1"). With the example if minority allele is in ]0.05:0.1] calling will become "./."')
	parser.add_option( '-c', '--MinAlCov',	dest='MinAlCov',	default='1',		help='Minimal read number of minor allele to call variant heterozygous (between 1 and infinity). [Default: %default]')
	parser.add_option( '-s', '--miss',		dest='miss',		default='0.2',		help='Maximal missing data proportion in the progeny (Excluding parents) (between 0 and 1). [Default: %default]')
	parser.add_option( '-p', '--pValue',	dest='pValue',		default='0.0001',	help='P-value threshold to keep marker (between 0 and 1). [Default: %default]')
	parser.add_option( '-T', '--pop',		dest='pop',			default='BiP',		help='Population type (Possible values: SELFPOL, SELF, BiP). [Default: %default]')
	parser.add_option( '-o', '--prefix',	dest='prefix',		default='Pop',		help='Prefix for output files. [Default: %default]')
	parser.add_option( '-a', '--addcov',	dest='addcov',		default='n',		help='Add coverage information to each genotype determined (y or n). [Default: %default]')
	parser.add_option( '-d', '--drawplot',	dest='drawplot',	default='n',		help='Draw statistic plot (y or n). [Default: %default]')
	parser.add_option( '-r', '--remove',	dest='remove',		default=None,		help='String to remove from marker name (may be too long for JoinMap format). Bay default marker name is "chromosome name"+"M"+"site position"')
	(options, args) = parser.parse_args()
	
	if options.vcf == None:
		sys.exit('Please provide a vcf file to -v options')
	
	if options.remove == None:
		to_replace = ""
	else:
		to_replace = options.remove
	
	if not(options.ref == None):
		sequence_dict = SeqIO.index(options.ref, "fasta")
		outTag = open(options.prefix+'_tags.fasta', 'w')
		outTag.close()
	
	# GLOBAL VARIABLES
	GLOB_allele_ratio = []
	GLOB_coverage = []
	GLOB_missing_per_marker = []
	GLOB_Khi2 = []
	DUPLICATED_markers = set()
	
	# identification of duplicated position
	file = open(options.vcf)
	all_positions = set()
	for line in file:
		data = line.split()
		if data:
			if data[0][0] != "#":
				what_the_fuck = ('M'.join(data[0:2])).replace(to_replace, "")
				if what_the_fuck in all_positions:
					DUPLICATED_markers.add(what_the_fuck)
				all_positions.add(what_the_fuck)
	
	# Getting parent name
	PARENT = []
	if not(options.parent == None):
		PARENT = options.parent.split(':')
	
	# To get final statistics
	DICO_FINAL_STAT = {}
	DICO_FINAL_STAT['unknown'] = 0
	DICO_FINAL_STAT['Bridge'] = 0
	for n in PARENT:
		DICO_FINAL_STAT[n] = 0
	
	# Getting accession to exclude from the output
	TO_EXCLUDE = set()
	if not(options.exclude == None):
		file = open(options.exclude)
		for line in file:
			data = line.split()
			TO_EXCLUDE.add(data[0])
		file.close()
	
	# Getting accession not to use in test
	NO_USED = set()
	if not(options.NoUsed == None):
		file = open(options.NoUsed)
		for line in file:
			data = line.split()
			NO_USED.add(data[0])
		file.close()
	
	# Getting options
	MAX_MISSING = float(options.miss)
	PVALUE = float(options.pValue)
	POP_TYPE = options.pop
	
	NO_READ, NO_COV, TOO_COV, NO_FREQ, NO_COV_PLUS, TRONCATED_FORMAT = 0, 0, 0, 0, 0, 0
	FLTR_MISS = 0
	FLTR_PVAL = 0
	PASSED = 0
	TO_MUCH_ALLELE = 0
	
	MORE_THAN_DI_ALLELE, ALL_MISSING, MONOMORPHOUS, DIMORPHOUS_SITES, MORE_THAN_DIMORPHOUS = 0, 0, 0, 0, 0
	
	# Creating filtered output file
	outvcf = open(options.prefix+'_sub.vcf','w')
	
	# Working line by line
	file = open(options.vcf)
	for line in file:
		data = line.split()
		if data[0] == "#CHROM":
			outvcf.write(line)
			header = data
			Accession_start = header.index('FORMAT')+1
			Accession_header = header[Accession_start:]
			
			# Getting accession to perform filter on
			dico_acc = set()
			dico_acc_to_genotype = set()
			dico_acc_to_treat = set()
			for n in Accession_header:
				if not(n in TO_EXCLUDE) and not (n in NO_USED) and not (n in PARENT):
					dico_acc.add(n)
				if not(n in TO_EXCLUDE) and not (n in PARENT):
					dico_acc_to_genotype.add(n)
				if not(n in TO_EXCLUDE) and not (n in NO_USED):
					dico_acc_to_treat.add(n)
			
			# Creating report output file
			outfile1 = open(options.prefix+'.tab', 'w')
			mot = ['marker', 'chomosome', 'position']
			for n in Accession_header:
				if not(n in TO_EXCLUDE):
					mot.append(n)
			mot.append('P-value')
			mot.append('ChiSquare\n')
			outfile1.write('\t'.join(mot))
			
			# Creating onemap output files
			outfile_onemapUn = open(options.prefix+'_onemap_unknown.tab', 'w')
			outfile_onemapBridge = open(options.prefix+'_onemap_Bridge.tab', 'w')
			for PID in PARENT:
				outfile_onemap = open(options.prefix+'_onemap_'+PID+'.tab', 'w')
				outfile_onemap.close()
			
			# Creating Joinmap output files
			outfile_JMUn = open(options.prefix+'_JM_unknown.loc', 'w')
			outfile_JMBridge = open(options.prefix+'_JM_Bridge.loc', 'w')
			for PID in PARENT:
				outfile_JM = open(options.prefix+'_JM_'+PID+'.loc', 'w')
				outfile_JM.close()
			
			# Creating tabulated output files
			entete = ['Marker', 'coding', 'rephased']
			for acc in Accession_header:
				if acc in dico_acc_to_genotype:
					entete.append(acc)
			outfile_Un = open(options.prefix+'_tab_unknown.tab', 'w')
			outfile_Un.write('\t'.join(entete)+'\n')
			outfile_Bridge = open(options.prefix+'_tab_Bridge.tab', 'w')
			outfile_Bridge.write('\t'.join(entete)+'\n')
			for PID in PARENT:
				outfile = open(options.prefix+'_tab_'+PID+'.tab', 'w')
				outfile.write('\t'.join(entete)+'\n')
				outfile.close()
			
		elif data[0][0] != "#":
			alleles = [data[header.index('REF')]]
			alleles = alleles + data[header.index('ALT')].split(',')
			chrom = data[header.index('#CHROM')]
			pos = data[header.index('POS')]
			format = data[header.index('FORMAT')].split(':')
			# Filtering on coverage and allele frequency
			liste = filter_on_read_cov(data[Accession_start:], format, int(options.MinCov), int(options.MaxCov), int(options.MinAlCov), options.WinFreq, GLOB_allele_ratio, GLOB_coverage, TO_EXCLUDE, Accession_header, dico_acc_to_treat)
			NO_READ += liste[1]
			NO_COV += liste[2]
			TOO_COV += liste[3]
			NO_FREQ += liste[4]
			NO_COV_PLUS += liste[5]
			TRONCATED_FORMAT += liste[6]
			MORE_THAN_DI_ALLELE += liste[7]
			ALL_MISSING += liste[8]
			MONOMORPHOUS += liste[9]
			DIMORPHOUS_SITES += liste[10]
			MORE_THAN_DIMORPHOUS += liste[11]
			TO_DECODE = liste[12]
			# Calculating missing data proportion
			MISSING = get_missing(liste[0], Accession_header, dico_acc)
			if MISSING < 1:
				GLOB_missing_per_marker.append(MISSING)
			# Filtering on missing data
			if MISSING <= MAX_MISSING:
				Khi_stat = KhiSquare(liste[0], Accession_header, dico_acc, POP_TYPE)
				GLOB_Khi2.append(math.log(Khi_stat[0], 10))
				if Khi_stat[0] >= PVALUE:
					PASSED += 1
					id_parent = ident_parent(Accession_header, liste[0], PARENT, Khi_stat[2], POP_TYPE)
					marker_name = chrom+'M'+pos
					marker_name = marker_name.replace(to_replace, "")
					if not (marker_name in DUPLICATED_markers):
						if not(options.ref == None):
							get_tags(sequence_dict, chrom, int(pos), options.prefix+'_tags.fasta', marker_name)
						DICO_FINAL_STAT[id_parent[0]] += 1
						if id_parent[0] == 'unknown':
							outfile_onemapUn.write(recode2onemap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
							outfile_JMUn.write(recode2JoinMap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name, str(DICO_FINAL_STAT[id_parent[0]]))+'\n')
							outfile_Un.write(recode2tab(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
						elif id_parent[0] == 'Bridge':
							outfile_onemapBridge.write(recode2onemap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
							outfile_JMBridge.write(recode2JoinMap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name, str(DICO_FINAL_STAT[id_parent[0]]))+'\n')
							outfile_Bridge.write(recode2tab(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
						else:
							outfile_onemap = open(options.prefix+'_onemap_'+id_parent[0]+'.tab', 'a')
							outfile_onemap.write(recode2onemap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
							outfile_onemap.close()
							
							outfile_JM = open(options.prefix+'_JM_'+id_parent[0]+'.loc', 'a')
							outfile_JM.write(recode2JoinMap(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name, str(DICO_FINAL_STAT[id_parent[0]]))+'\n')
							outfile_JM.close()
							
							outfile = open(options.prefix+'_tab_'+id_parent[0]+'.tab', 'a')
							outfile.write(recode2tab(Accession_header, liste[0], Khi_stat[2], POP_TYPE, dico_acc_to_genotype, id_parent[1], marker_name)+'\n')
							outfile.close()
						outvcf.write(line)
						genotype = decode(Accession_header, TO_DECODE, alleles, TO_EXCLUDE, data[Accession_start:], format, options.addcov)
						outfile1.write('\t'.join([chrom+'M'+pos, chrom, pos]+genotype+list(map(str, Khi_stat[0:2])))+'\n')
				else:
					FLTR_PVAL += 1
			else:
				FLTR_MISS += 1
		else:
			outvcf.write(line)
	file.close()
	outfile_onemapUn.close()
	outfile_JMUn.close()
	outfile_Un.close()
	outfile_onemapBridge.close()
	outfile_JMBridge.close()
	outfile_Bridge.close()
	
	# Adding end lines to Joimap files
	for n in DICO_FINAL_STAT:
		outfile = open(options.prefix+'_JM_'+n+'.loc','a')
		outfile.write('\nindividual names:\n')
		for acc in Accession_header:
			if acc in dico_acc_to_genotype:
				outfile.write(acc+'\n')
		outfile.close()
	
	# Adding header in onemap files
	for n in DICO_FINAL_STAT:
		to_add = str(len(dico_acc_to_genotype)) + ' ' + str(DICO_FINAL_STAT[n]) + '\n'
		outfile = open(options.prefix+'_onemap_'+n+'.tab','r')
		total = to_add + outfile.read()
		outfile.close()
		
		outfile = open(options.prefix+'_onemap_'+n+'.tab','w')
		outfile.write(total)
		outfile.close()
	
	# Adding header to Joimap files
	for n in DICO_FINAL_STAT:
		to_add = time.strftime(';Any %m %d %H:%M:%S CET %y', time.localtime())+'\n\nname= '+options.prefix+'\npopt= CP\nnloc= '+str(DICO_FINAL_STAT[n])+'\nnind= '+str(len(dico_acc_to_genotype))+'\n\n'
		outfile = open(options.prefix+'_JM_'+n+'.loc','r')
		total = to_add + outfile.read()
		outfile.close()
		
		outfile = open(options.prefix+'_JM_'+n+'.loc','w')
		outfile.write(total)
		outfile.close()
		
	
	# Writting global statistics
	outfile_stat = open(options.prefix+'_report.tab', 'w')
	outfile_stat.write("PASSED OPTIONS\n")
	outfile_stat.write('--vcf: '+options.vcf+'\n')
	if not(options.parent == None):
		outfile_stat.write('--parent: '+options.parent+'\n')
	if not(options.exclude == None):
		outfile_stat.write('--exclude: '+options.exclude+'\n')
	if not(options.NoUsed == None):
		outfile_stat.write('--NoUsed: '+options.NoUsed+'\n')
	outfile_stat.write('--MinCov: '+options.MinCov+'\n')
	outfile_stat.write('--MaxCov: '+options.MaxCov+'\n')
	outfile_stat.write('--WinFreq: '+options.WinFreq+'\n')
	outfile_stat.write('--MinAlCov: '+options.MinAlCov+'\n')
	outfile_stat.write('--miss: '+options.miss+'\n')
	outfile_stat.write('--pValue: '+options.pValue+'\n')
	outfile_stat.write('--pop: '+options.pop+'\n')
	outfile_stat.write('--prefix: '+options.prefix+'\n')
	outfile_stat.write("\n=============================================\nREPORT ON DATA POINTS\n")
	outfile_stat.write('Removed line due to duplicated position in vcf file: '+str(len(DUPLICATED_markers)*2)+'\n')
	outfile_stat.write('Converted to missing due to no reads: '+str(NO_READ)+'('+str(NO_READ/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write('Converted to missing due to insufficient reads coverage: '+str(NO_COV)+' ('+str(NO_COV/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write('Converted to missing due to overcoverage: '+str(TOO_COV)+'('+str(TOO_COV/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write('Converted to missing due to too many variant: '+str(MORE_THAN_DI_ALLELE)+'('+str(MORE_THAN_DI_ALLELE/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write('Converted to missing due to ambiguous minority allele frequency (--WinFreq argument): '+str(NO_FREQ)+' ('+str(NO_FREQ/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write('Converted to missing due to ambiguous minority allele coverage (--MinAlCov argument): '+str(NO_COV_PLUS)+' ('+str(NO_COV_PLUS/float((len(Accession_header)-len(TO_EXCLUDE))*(FLTR_MISS+FLTR_PVAL+PASSED))*100)+'%)\n')
	outfile_stat.write("\n=============================================\nREPORT MARKER POINT\n")
	outfile_stat.write('Total marker treated: '+str(FLTR_MISS+FLTR_PVAL+PASSED)+'\n')
	
	outfile_stat.write('Complete line coverted to missing due to unexpected vcf format (no AD or GT or DP tags): '+str(TRONCATED_FORMAT)+' ('+str(TRONCATED_FORMAT/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tMarker with missing data on all individuals: '+str(ALL_MISSING)+' ('+str(ALL_MISSING/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tMonomorphous makers (converted to missing for convenience): '+str(MONOMORPHOUS)+' ('+str(MONOMORPHOUS/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tDi-allelic markers: '+str(DIMORPHOUS_SITES)+' ('+str(DIMORPHOUS_SITES/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tMore than di-allelic markers (converted to missing for convenience): '+str(MORE_THAN_DIMORPHOUS)+' ('+str(MORE_THAN_DIMORPHOUS/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	
	outfile_stat.write('\tMarker removed based on missing data cutoff (--miss argument): '+str(FLTR_MISS)+' ('+str(FLTR_MISS/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tMarker removed based on ChiSquare cutoff (--pValue argument): '+str(FLTR_PVAL)+' ('+str(FLTR_PVAL/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	outfile_stat.write('\tSelected marker: '+str(PASSED)+' ('+str(PASSED/float(FLTR_MISS+FLTR_PVAL+PASSED)*100)+'%)\n')
	for n in DICO_FINAL_STAT:
		outfile_stat.write('Marker parsed in '+n+' file: '+str(DICO_FINAL_STAT[n])+'\n')
	
	# To draw plot
	
	if options.drawplot == 'y':
		GLOB_allele_ratio = np.array(GLOB_allele_ratio)
		GLOB_coverage = np.array(GLOB_coverage)
		GLOB_missing_per_marker = np.array(GLOB_missing_per_marker)
		GLOB_Khi2 = np.array(GLOB_Khi2)
		draw_statistics_plot(options.prefix+'.pdf', GLOB_allele_ratio, GLOB_coverage, GLOB_missing_per_marker, GLOB_Khi2, options.WinFreq, PVALUE, MAX_MISSING, int(options.MinCov))

		
		
	
if __name__ == "__main__": __main__()
