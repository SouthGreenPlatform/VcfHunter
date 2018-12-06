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
import configparser
import datetime
import optparse
import tempfile
import sys
import os
import traceback
import shutil
import time
import multiprocessing as mp
import subprocess
import threading
import random
import math
import gzip
from inspect import currentframe, getframeinfo
from operator import itemgetter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import utilsSR.utilsSR as utils
sys.stdout.write("module loaded\n")
sys.stdout.flush()

def run_analysis(LIB_DIC, ACC_ID, PLOIDY, OPTIONS, LOCA_PROGRAMS, CONFIG, DICO_CHR, QUEUE, PATHNAME, PARSEUNMAPPED):

	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	
	REF = CONFIG.get('Reference', 'genome')
	SAMTOOLS = LOCA_PROGRAMS.get('Programs','samtools')
	BWA = LOCA_PROGRAMS.get('Programs','bwa')
	JAVA = LOCA_PROGRAMS.get('Programs','java')
	PICARD = LOCA_PROGRAMS.get('Programs','picard')
	GATK = LOCA_PROGRAMS.get('Programs','gatk')
	UseUnifiedGenotyperForBaseRecal =  'no'
	PYTHON = LOCA_PROGRAMS.get('Programs','python')
	PREFIX = OPTIONS.prefix
	if CONFIG.has_section('Variant'):
		if CONFIG.has_option('Variant', 'UseUnifiedGenotyperForBaseRecal'):
			UseUnifiedGenotyperForBaseRecal =  CONFIG.get('Variant', 'UseUnifiedGenotyperForBaseRecal')
	
	if 'a' in OPTIONS.steps:
		#1 Mapping
		#2 Merging
		to_return = utils.run_step_A(ACC_ID, LIB_DIC, BWA, REF, TMP, JAVA, PICARD, SAMTOOLS, PREFIX, QUEUE, PARSEUNMAPPED)
	
	if 'b' in OPTIONS.steps:
		#3 removing duplicates
		to_return = utils.run_step_B(JAVA, PICARD, ACC_ID, TMP, PREFIX, QUEUE)
	
	if 'c' in OPTIONS.steps:
		#4 indel realignment
		to_return = utils.run_step_C(ACC_ID, JAVA, GATK, REF, CONFIG, PREFIX, QUEUE)
	
	if 'd' in OPTIONS.steps:
		#5 Base recalibration
		to_return = utils.run_step_D(CONFIG, ACC_ID, UseUnifiedGenotyperForBaseRecal, JAVA, GATK, REF, PLOIDY, PREFIX, QUEUE)
	
	if 'e' in OPTIONS.steps:
		#6 GVCF generation
		to_return = utils.run_step_E(ACC_ID, PYTHON, REF, DICO_CHR, PREFIX, QUEUE, PATHNAME)
	if to_return == 0:
		return 0
	else:
		return to_return

def main_run_analysis(job):

	try:
		rslt = run_analysis(job[0],job[1],job[2],job[3],job[4],job[5], job[6], job[7], job[8], job[9])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def main_combine(job):

	try:
		rslt = utils.create_pseudo_VCF(job[0],job[1],job[2],job[3],job[4],job[5],job[6],job[7])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def main_prepa_combine(job):

	try:
		rslt = utils.createGVCF(job[0],job[1],job[2])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def main_combine_Large(job):

	try:
		rslt = utils.create_pseudo_VCF_Large(job[0],job[1],job[2],job[3],job[4],job[5],job[6],job[7])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def main_merging_sub_vcf(job):

	try:
		rslt = utils.merge_sub_vcf(job[0],job[1],job[2],job[3])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\nThis program go through all steps needed to call SNP from filtered fastq file (From mapping to SNP calling).")
	# Wrapper options.
	parser.add_option( '-c', '--conf',		dest='conf',		default=None,		help='A configuration file containing path to references fastq files.'
	'The conf file should contain 2 sections ([Libraries] and [Reference]) and 2 additional ones ([Mapping], [Variant]).\t\t\t\t\t[Libraries] '
	'section should look like as follows :\t[Libraries]\t\t\t\t\t\tlib1 = genome_name path_to_mate1 path_to_mate2 ploidy\t\tlib2 = genome_name '
	'path_to_single ploidy\t\t\t...\t\t\t\t\t\t[Reference] section should look like as follows :\t[Reference]\t\t\t\t\t\tgenome = path_to_the_reference_sequence\t\t'
	'[Variant] section may contain 4 options and should look like :\t\t\t\t\t\t[Variant]\t\t\t\t\t\tindel = path to vcf_of_known_indels\t\t\t\t\tsnp = '
	'path to vcf_of_known_SNPs\t\t\t\tHCopt = additional options to pass to HaplotypeCaller\t\tUseUnifiedGenotyperForBaseRecal = yes or no (if not filled default = no)')
	parser.add_option( '-t', '--thread',	dest='thread',		default='1',		help='Max number of accessions treated at the same time (integer), [default: %default]')
	parser.add_option( '-q', '--queue',		dest='queue',		default=None,		help='Queue to use if SGE is installed on your machine. Do not fill otherwise, [default: %default]')
	parser.add_option( '-p', '--prefix',	dest='prefix',		default='All_lib',	help='Prefix for output bam and vcf containing all libraries. [default: %default]')
	parser.add_option( '-g', '--outgzip',	dest='outgzip',		default='n',		help='Output files in gzip format. [Default: %default]')
	parser.add_option( '-k', '--keepUn',	dest='keepUn',		default='n',		help='Keep unmapped reads in a file. [Default: %default]')
	parser.add_option( '-C', '--chrom',	dest='chrom',		default='all',		help='Chromosomes to work with (only for step f). If "all", all chromosomes will be used for calling. Either: a list of chromosome names separated by ":" [Default: %default]')
	parser.add_option( '-s', '--steps',		dest='steps',		default=None,		help='A string containing steps to perform:\t\t\t\t'
	'a: Aligning libraries\t\t\t\t\t'
	'b: Removing duplicates\t\t\t\t\t\t\t'
	'c: Indel realignment\t\t\t\t\t\t'
	'd: Bases recalibration\t\t\t\t\t\t'
	'e: Allele counting\t\t\t\t\t\t'
	'f: Genotype calling\t\t\t\t\t'
	'g: Merging genotype calling\t\t\t\t\t'
	'h: Mapping statistics calculation\t\t\t\t\t')
	(options, args) = parser.parse_args()
	
	if options.conf == None:
		sys.exit('--conf argument is missing')
	if options.steps == None:
		sys.exit('--steps argument is missing')
	
	nbProcs = int(options.thread)
	
	#Loading the file locating programs
	pathname = os.path.dirname(sys.argv[0])
	loca_programs = configparser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	#Loading the configuration file
	config = configparser.RawConfigParser()
	config.read(options.conf)
	
	ref = config.get('Reference', 'genome')
	samtools = loca_programs.get('Programs','samtools')
	bwa = loca_programs.get('Programs','bwa')
	java = loca_programs.get('Programs','java')
	picard = loca_programs.get('Programs','picard')
	gatk = loca_programs.get('Programs','gatk')
	python = loca_programs.get('Programs','python')
	vcfConcat = loca_programs.get('Programs','vcfconcat')
	
	# Checking options
	if options.keepUn != 'y' and options.keepUn != 'n':
		sys.exit('Please enter either "y" or "n" to --keepUn options\n')
	
	# Getting accessions to work with
	dico_lib = {}
	dico_ploidy = {}
	for n in config.options('Libraries'):
		data =  config.get('Libraries', n).split()
		if not(data[0] in dico_lib):
			dico_lib[data[0]] = {}
		dico_lib[data[0]][n] = data[1:-1]
		if data[0] in dico_ploidy:
			if dico_ploidy[data[0]] != data[-1]:
				sys.exit('This is embarrassing... The accession '+data[0]+' has two different ploidy level??')
		else:
			dico_ploidy[data[0]] = data[-1]
	
	# Checking files
	if 'a' in options.steps:
		BugInFile = 0
		for n in dico_lib:
			for lib in dico_lib[n]:
				for file in dico_lib[n][lib]:
					if not(os.path.isfile(file)):
						BugInFile = 1
						sys.stdout.write('This is embarrassing: the file '+file+' for library '+lib+' has not been found\n')
		if BugInFile:
			sys.exit('The program exited because some files were not found...')
	
	# Reference indexation
	sys.stdout.write('Checking for reference sequence index for picard ...\n')
	if not(os.path.isfile(ref+'.fai')):
		create_index = '%s faidx %s' % (samtools, ref)
		utils.run_job(getframeinfo(currentframe()), create_index, 'Error in creating .fai:\n')
		sys.stdout.write('Done\n')
	else:
		sys.stdout.write('Index found for the reference sequence.\n')
	
	sys.stdout.write('Checking for reference sequence index ...\n')
	if not(os.path.isfile(ref+'.amb')):
		sys.stdout.write('No index found for the reference sequence. Indexing ...\n')
		index = '%s index -a bwtsw %s 2>/dev/null' % (bwa, ref)
		utils.run_job(getframeinfo(currentframe()), index, 'Error in creating index for bwa:\n')
		sys.stdout.write('Done\n')
	else:
		sys.stdout.write('Index found for the reference sequence.\n')
	
	sys.stdout.write('Checking for reference sequence .dict for picard ...\n')
	if not(os.path.isfile('.'.join(ref.split('.')[:-1])+'.dict')):
		create_dict = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s' % (java, picard, ref, '.'.join(ref.split('.')[:-1])+'.dict')
		utils.run_job(getframeinfo(currentframe()), create_dict, 'Error in creating .dict:\n')
		sys.stdout.write('Done\n')
	else:
		sys.stdout.write('.dict found for the reference sequence.\n')
	sys.stdout.flush()
	
	# Obtaining chromosome informations
	dico_chr = {}
	if 'e' in options.steps or 'f' in options.steps or 'g' in options.steps or 'E' in options.steps or 'F' in options.steps:
		# calculating sequence length
		sequence_dict = SeqIO.index(ref, "fasta")
		for n in sequence_dict:
			dico_chr[n] = len(str(sequence_dict[n].seq))
		del sequence_dict
	
	# Running analysis
	if 'a' in options.steps or 'b' in options.steps or 'c' in options.steps or 'd' in options.steps or 'e' in options.steps:
		# threads = []
		listJobs = []
		for n in dico_lib:
			sys.stdout.flush()
			listJobs.append([dico_lib[n], n, dico_ploidy[n], options, loca_programs, config, dico_chr, options.queue, pathname, options.keepUn])
		
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_run_analysis, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
	
	# Creating pseudo VCF
	if 'f' in options.steps:
		
		# getting accession list
		liste_accessions = sorted(list(dico_lib.keys()))
		
		# Selecting chromosomes
		if options.chrom == 'all':
			ChrList = list(dico_chr.keys())
		else:
			ChrList = options.chrom.split(':')
		
		# Optimizing window size
		totalSize = 0
		for chr in ChrList:
			if not chr in dico_chr:
				sys.exit('This is embarrassing... chromosome: '+chr+' has not been found in the reference fasta file...\n')
			totalSize += dico_chr[chr]
		
		# window = int(totalSize/float(nbProcs))
		window = 1000000
		dicoWindow = {}
		for chr in ChrList:
			dicoWindow[chr] = []
			start = 0
			while start < dico_chr[chr]:
				if start + window <= dico_chr[chr]:
					dicoWindow[chr].append([start, start + window])
				else:
					dicoWindow[chr].append([start, dico_chr[chr]])
				start += window
		
		listJobs = []
		for chr in dicoWindow:
			for pos in dicoWindow[chr]:
				listJobs.append([liste_accessions, ref, options.prefix, dico_ploidy, dico_chr, chr, pos[0], pos[1]])
		
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_combine, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
		
		# Merging files by chromosomes
		listJobs = []
		for chr in dicoWindow:
			listJobs.append([options.prefix, chr, dicoWindow[chr], options.outgzip])
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_merging_sub_vcf, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
		
		sys.stdout.write('Setp f finished\n')
	
	# Merging pseudo VCF
	if 'g' in options.steps:
		utils.merge_vcf(options.prefix, dico_chr)
	
	# Mapping statistics calculation
	if 'h' in options.steps:
		utils.Calc_stats(options.prefix, dico_lib)
		
	# Merging GVCF into one large GVCF (For very large data set: greater than the limit allowed by unix)
	if 'E' in options.steps:
		
		# getting accession list
		liste_accessions = sorted(list(dico_lib.keys()))
		
		# Selecting chromosomes
		if options.chrom == 'all':
			ChrList = list(dico_chr.keys())
		else:
			ChrList = options.chrom.split(':')
		
		# Preparing chromosomal launching
		listJobs = []
		for chr in ChrList:
			listJobs.append([chr, liste_accessions, options.prefix])
		
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_prepa_combine, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
	
	# Creating the vcf
	if 'F' in options.steps:
		
		# getting accession list
		liste_accessions = sorted(list(dico_lib.keys()))
		
		# Selecting chromosomes
		if options.chrom == 'all':
			ChrList = list(dico_chr.keys())
		else:
			ChrList = options.chrom.split(':')
		
		# Optimizing window size
		totalSize = 0
		for chr in ChrList:
			if not chr in dico_chr:
				sys.exit('This is embarrassing... chromosome: '+chr+' has not been found in the reference fasta file...\n')
			totalSize += dico_chr[chr]
		
		# window = int(totalSize/float(nbProcs))
		window = 1000000
		dicoWindow = {}
		for chr in ChrList:
			dicoWindow[chr] = []
			start = 0
			while start < dico_chr[chr]:
				if start + window <= dico_chr[chr]:
					dicoWindow[chr].append([start, start + window])
				else:
					dicoWindow[chr].append([start, dico_chr[chr]])
				start += window
		
		listJobs = []
		for chr in dicoWindow:
			for pos in dicoWindow[chr]:
				listJobs.append([liste_accessions, ref, options.prefix, dico_ploidy, dico_chr, chr, pos[0], pos[1]])
		
		# for n in listJobs:
			# utils.create_pseudo_VCF_Large(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7])
		
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_combine_Large, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
		
		# Merging files by chromosomes
		listJobs = []
		for chr in dicoWindow:
			listJobs.append([options.prefix, chr, dicoWindow[chr], options.outgzip])
		pool = mp.Pool(processes=nbProcs)
		results = pool.map(main_merging_sub_vcf, listJobs)
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
		
		sys.stdout.write('Setp F finished\n')
	
		
if __name__ == "__main__": __main__()
