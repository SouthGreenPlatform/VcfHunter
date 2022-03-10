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
import optparse
import os
import shutil
import subprocess
import sys
import tempfile
import fileinput
import configparser
import operator
import time
import random
import datetime
import glob
import multiprocessing as mp
import subprocess
import threading
import gzip
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from inspect import currentframe, getframeinfo
import utilsSR.utilsSR as utils
sys.stdout.write("module loaded\n")
sys.stdout.flush()


def run_D_to_J(STEPS, LIB, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, JAVA, PICARD, REF, GATK, PLOIDY, UseUnifiedGenotyperForBaseRecal, BAMTOOLS, PYTHON, QUEUE, PATHNAME, DICO_CHR):
	if 'd' in STEPS:
		to_return = utils.run_step_D_RNAseq(LIB, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'e' in STEPS:
		to_return = utils.run_step_E_RNAseq(LIB, ACC, PREFIX, JAVA, PICARD, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'f' in STEPS:
		to_return = utils.run_step_F_RNAseq(ACC, JAVA, PICARD, PREFIX, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'g' in STEPS:
		to_return = utils.run_step_G_RNAseq(JAVA, PICARD, ACC, REF, PREFIX, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'h' in STEPS:
		to_return = utils.run_step_H_RNAseq(JAVA, GATK, ACC, REF, PREFIX, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'i' in STEPS:
		to_return = utils.run_step_I_RNAseq(ACC, JAVA, GATK, REF, PLOIDY, PREFIX, UseUnifiedGenotyperForBaseRecal, QUEUE)
		if to_return == 0:
			pass
		else:
			return to_return
	if 'j' in STEPS:
		to_return = utils.run_step_E(ACC, PYTHON, REF, DICO_CHR, PREFIX, QUEUE, PATHNAME)
		if to_return == 0:
			pass
		else:
			return to_return
	return to_return

def main_run_analysis(job):

	try:
		rslt = run_D_to_J(job[0],job[1],job[2],job[3],job[4],job[5],job[6],job[7],job[8],job[9],job[10],job[11],job[12],job[13],job[14],job[15], job[16], job[17], job[18])
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

def main_merging_sub_vcf(job):

	try:
		rslt = utils.merge_sub_vcf(job[0],job[1],job[2], job[3])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\nThis program go through all steps needed to call SNP from RNAseq dataset. From mapping to SNP calling.")
	# Wrapper options.
	parser.add_option( '-c', '--conf', dest='conf', default=None, help='A configuration file containing path to references fastq files.'
	'The conf file should contain 4 sections ([Libraries], [Reference], [star] and [General]).\t\t\t\t\t[Libraries] '
	'section should look like as follows :\t[Libraries]\t\t\t\t\t\tlib1 = genome_name path_to_mate1 path_to_mate2 ploidy\t\tlib2 = genome_name '
	'path_to_single ploidy\t\t\t...\t\t\t\t\t\t[Reference] section should look like as follows :\t[Reference]\t\t\t\t\t\tgenome = path_to_the_reference_sequence\t\t\t'
	'[star] section should contain star options and should look like :\t\t\t\t\t\t[star]\t\t\t\t\t\toptions = additional_options_to_pass\t\t\t'
	'[General] section may contain two options and should look like :\t\t\t\t\t[General]\t\t\t\t\t\tmax_size = max_read_length\t\t\t\t\tgff3 (optional) = '
	'path to a gff3 file used to calculate statistics on genes coverage\t\tUseUnifiedGenotyperForBaseRecal = yes or no (if not filled default = no)')
	parser.add_option( '-t', '--thread', dest='thread', default='1', help='Max number of accessions treated at the same time (integer), [default: %default]')
	parser.add_option( '-q', '--queue', dest='queue', default=None, help='Queue to use if SGE is installed on your machine. Do not fill otherwise, [default: %default]')
	parser.add_option( '-p', '--prefix', dest='prefix', default='All_lib', help='Prefix for output bam and vcf containing all libraries. [default: %default]')
	parser.add_option( '-g', '--outgzip',	dest='outgzip',		default='n',		help='Output files in gzip format. [Default: %default]')
	parser.add_option( '-s', '--steps', dest='steps', default=None, help='A string containing steps to perform:\t\t\t\t'
	'a: Indexing reference for first mapping step\t\t\t'
	'b: Merging fastq and first mapping step to generate splicing sites\t\t\t\t\t\t'
	'c: Indexing reference for second mapping step\t\t\t'
	'd: Aligning libraries\t\t\t\t\t\t'
	'e: Merging bam libraries with identical identifier\t\t'
	'f: Removing duplicates\t\t\t\t\t\t'
	'g: Reordering reads\t\t\t\t\t\t'
	'h: Splitting and trimming reads\t\t\t\t\t'
	'i: Indel realignment\t\t\t\t\t'
	'j: Allele counting\t\t\t\t\t\t'
	'k: Genotype calling\t\t\t\t\t'
	'l: Merging genotype calling\t\t\t\t\t'
	'm: Gene exon coverage statistics calculation\t')
	(options, args) = parser.parse_args()
	
	if options.conf == None:
		sys.exit('--conf argument is missing')
	if options.steps == None:
		sys.exit('--septs argument is missing')
	
	#Loading the file locating programs
	PATHNAME = os.path.dirname(sys.argv[0])
	LOCA_PROGRAMS = configparser.RawConfigParser()
	LOCA_PROGRAMS.read(PATHNAME+'/loca_programs.conf')
	
	#Loading the configuration file
	config = configparser.RawConfigParser()
	config.read(options.conf)
	
	STAR = LOCA_PROGRAMS.get('Programs','star')
	OUT_REF1 = options.prefix+'_ref_star_1'
	OUT_REF2 = options.prefix+'_ref_star_2'
	REF = config.get('Reference', 'genome')
	SJDBO = str(config.getint('General', 'max_size') - 1)
	PREFIX = options.prefix
	SJDBFCSE = PREFIX+'_JUNC_ESTIMATION_SJ.out.tab'
	STAR_OPT = config.get('star', 'options')
	PROC = options.thread
	JAVA = LOCA_PROGRAMS.get('Programs','java')
	PICARD = LOCA_PROGRAMS.get('Programs','picard')
	SAMTOOLS = LOCA_PROGRAMS.get('Programs','samtools')
	GATK = LOCA_PROGRAMS.get('Programs','gatk')
	BAMTOOLS = LOCA_PROGRAMS.get('Programs','bamtools')
	PYTHON = LOCA_PROGRAMS.get('Programs','python')
	QUEUE = options.queue
	nbProcs = int(options.thread)
	
	UseUnifiedGenotyperForBaseRecal =  'no'
	if config.has_section('General'):
		if config.has_option('General', 'UseUnifiedGenotyperForBaseRecal'):
			UseUnifiedGenotyperForBaseRecal =  config.get('General', 'UseUnifiedGenotyperForBaseRecal')
	
	if not(os.path.isfile('.'.join(REF.split('.')[:-1])+'.dict')):
		sys.stdout.write('No .dict file found for the reference. Running indexation...')
		create_dict = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s' % (JAVA, PICARD, REF, '.'.join(REF.split('.')[:-1])+'.dict')
		if QUEUE == None:
			utils.run_job(getframeinfo(currentframe()), create_dict, 'Error when creating the .dict file:\n')
		else:
			utils.run_qsub(QUEUE, [create_dict], 1, PREFIX+'-Dict', "12G", PREFIX)
		sys.stdout.write('Done\n')
	
	if not(os.path.isfile(REF+'.fai')):
		sys.stdout.write('No .fai file found for the reference. Running indexation...')
		create_index = '%s faidx %s' % (SAMTOOLS, REF)
		if QUEUE == None:
			utils.run_job(getframeinfo(currentframe()), create_index, 'Error when creating the .fai file:\n')
		else:
			utils.run_qsub(QUEUE, [create_index], 1, PREFIX+'-Index', "12G", PREFIX)
		sys.stdout.write('Done\n')
	
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
	
	# Obtaining chromosome informations
	dico_chr = {}
	if 'j' in options.steps or 'k' in options.steps or 'l' in options.steps:
		# calculating sequence length
		list_chr_order = []
		sequence_dict = list(SeqIO.parse(REF, "fasta"))
		for i in range(len(sequence_dict)):
			list_chr_order.append([sequence_dict[i].id, len(str(sequence_dict[i].seq))])
			dico_chr[sequence_dict[i].id] = len(str(sequence_dict[i].seq))
		del sequence_dict
	
	# Reference indexation
	if 'a' in options.steps:
		if os.path.isdir(OUT_REF1):
			shutil.rmtree(OUT_REF1)
		os.mkdir(OUT_REF1)
		utils.run_step_A_RNAseq(STAR, PROC, OUT_REF1, REF, SJDBO, PREFIX, QUEUE)
	
	# Aligning reads for junction identification
	if 'b' in options.steps:
		TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
		utils.run_step_B_RNAseq(TMP, dico_lib, STAR, PROC, PREFIX, STAR_OPT, OUT_REF1, QUEUE)
		
	# Running reference indexation with intron/exon junction
	if 'c' in options.steps:
		if os.path.isdir(OUT_REF2):
			shutil.rmtree(OUT_REF2)
		os.mkdir(OUT_REF2)
		utils.run_step_C_RNAseq(STAR, PROC, OUT_REF2, REF, SJDBO, PREFIX, SJDBFCSE, QUEUE)
	
	# Running step d to j
	if 'd' in options.steps or 'e' in options.steps or 'f' in options.steps or 'g' in options.steps or 'h' in options.steps or 'i' in options.steps or 'j' in options.steps:
		# threads = []
		listJobs = []
		for acc in dico_lib:
			sys.stdout.flush()
			listJobs.append([options.steps, dico_lib[acc], acc, STAR, OUT_REF2, 4, STAR_OPT, PREFIX, JAVA, PICARD, REF, GATK, dico_ploidy[acc],UseUnifiedGenotyperForBaseRecal, BAMTOOLS, PYTHON, QUEUE, PATHNAME, dico_chr])
		
		pool = mp.Pool(processes=int(PROC))
		results = pool.map(main_run_analysis, listJobs)
		
		for n in results:
			if n != 0:
				sys.stderr.write('*******************************\n\n')
				sys.stderr.write(str(n[1])+'\n')
				sys.stderr.write('*******************************\n')
		
		NoBugg = 1
		for n in results:
			if n != 0:
				if n[0] == "Not applicable":
					sys.stdout.write('A problem occured. See error log for more information.\n')
				else:
					sys.stdout.write('A problem occured in '+n[0]+' accession. See error log for more information.\n')
				NoBugg = 0
		
		if NoBugg:
			if os.path.isdir(PREFIX+'_ref_star_2'):
				shutil.rmtree(PREFIX+'_ref_star_2')
		
	# Creating pseudo VCF
	if 'k' in options.steps:
		
		# getting accession list
		liste_accessions = list(dico_lib.keys())
		
		# Optimizing window size
		totalSize = 0
		for chr in dico_chr:
			totalSize += dico_chr[chr]
		# window = int(totalSize/float(nbProcs))
		window = 1000000
		dicoWindow = {}
		for chr in dico_chr:
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
				listJobs.append([liste_accessions, REF, PREFIX, dico_ploidy, list_chr_order, chr, pos[0], pos[1]])
		
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
	if 'l' in options.steps:
		utils.merge_vcf(options.prefix, dico_chr)
	
	# Collecting mapping statistics by accession
	if 'd' in options.steps:
		if not(os.path.isdir(PREFIX)):
			os.mkdir(PREFIX)
		utils.run_stat_step_D_RNAseq(PREFIX, config)
	
	# Collecting duplicate statistics by accession
	if 'f' in options.steps:
		if not(os.path.isdir(PREFIX)):
			os.mkdir(PREFIX)
		utils.run_stat_step_F_RNAseq(PREFIX, dico_lib)
	
	# calculating exon coverage
	if 'm' in options.steps:
		GFF3 = config.get('General', 'gff3')
		utils.run_step_M_RNAseq(GFF3, REF, LOCA_PROGRAMS, PREFIX, dico_lib)

	
if __name__ == "__main__": __main__()
