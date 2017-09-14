
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

import ConfigParser
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


def stop_err( msg ):
	sys.stderr.write( "%s\n" % msg )
	sys.exit()

def run_job (frameinfo, cmd_line, ERROR):
	print cmd_line
	try:
		tmp = tempfile.NamedTemporaryFile().name
		error = open(tmp, 'w')
		proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
		returncode = proc.wait()
		error.close()
		error = open( tmp, 'rb' )
		stderr = ''
		buffsize = 1048576
		try:
			while True:
				stderr += error.read( buffsize )
				if not stderr or len( stderr ) % buffsize != 0:
					break
		except OverflowError:
			pass
		error.close()
		os.remove(tmp)
		if returncode != 0:
			raise Exception, stderr
	except Exception, e:
		stop_err( 'Line : '+str(frameinfo.lineno)+' - '+ERROR + str( e ) )

def hold_job(ID_liste):
	time.sleep(10)
	qs=os.popen("qstat")
	nb_run = 0
	for n in qs:
		l = n.split()
		if len(l) > 2:
			if l[0] in ID_liste:
				nb_run = nb_run + 1
	while nb_run > 0:
		time.sleep(10)
		qs=os.popen("qstat")
		nb_run = 0
		for n in qs:
			l = n.split()
			if len(l) > 2:
				if l[0] in ID_liste:
					nb_run = nb_run + 1

def job_ID2list(MOT):
	for n in MOT:
		return(n.split()[0])

def run_qsub(QUEUE, COMMANDE_LIST, PROC, JID, LMEM, PREFIX):
	
	list_job = []
	for n in COMMANDE_LIST:
		if LMEM == None:
			if PROC > 1:
				sys.stdout.write("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'\n")
				qs=os.popen("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
			else:
				sys.stdout.write("qsub -q "+QUEUE+" -terse -b yes -V -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'\n")
				qs=os.popen("qsub -q "+QUEUE+" -terse -b yes -V -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
		else:
			if PROC > 1:
				sys.stdout.write("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'\n")
				qs=os.popen("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
			else:
				sys.stdout.write("qsub -q "+QUEUE+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'\n")
				qs=os.popen("qsub -q "+QUEUE+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o "+PREFIX+".txt -e "+PREFIX+"_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
	sys.stdout.flush()
	hold_job(list_job)

def run_step_A(ACC_ID, LIB_DIC, BWA, REF, TMP, JAVA, PICARD, SAMTOOLS, PREFIX, QUEUE):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	if not(os.path.isdir(ACC_ID)):
		os.mkdir(ACC_ID)
	#############
	#1 Mapping
	sys.stdout.write('Initiating mapping step.\n')
	to_mapp = []
	for n in LIB_DIC:
		sys.stdout.write('Mapping '+n+'....\n')
		rg_tag = '@RG"\t"ID:'+ACC_ID+'"\t"SM:'+ACC_ID+'"\t"LB:'+ACC_ID+'"\t"PU:whatever"\t"PL:ILLUMINA'
		if QUEUE == None:
			if len(LIB_DIC[n]) == 1:
				map = '%s mem -t %s -M %s -R %s %s > %s' % (BWA, 1, REF, rg_tag, LIB_DIC[n][0], ACC_ID+'/'+TMP+n+'.sam')
			elif len(LIB_DIC[n]) == 2:
				map = '%s mem -t %s -M %s -R %s %s %s > %s' % (BWA, 1, REF, rg_tag, LIB_DIC[n][0], LIB_DIC[n][1], ACC_ID+'/'+TMP+n+'.sam')
			else:
				return 'Problem in the configuration file in libraries section on '+ACC_ID+'\n'
		else:
			if len(LIB_DIC[n]) == 1:
				map = '%s mem -t %s -M %s -R %s %s > %s' % (BWA, 4, REF, rg_tag, LIB_DIC[n][0], ACC_ID+'/'+TMP+n+'.sam')
			elif len(LIB_DIC[n]) == 2:
				map = '%s mem -t %s -M %s -R %s %s %s > %s' % (BWA, 4, REF, rg_tag, LIB_DIC[n][0], LIB_DIC[n][1], ACC_ID+'/'+TMP+n+'.sam')
			else:
				return 'Problem in the configuration file in libraries section on '+ACC_ID+'\n'
		to_mapp.append(map)
	
	if QUEUE == None:
		for comd in to_mapp:
			sys.stdout.write(comd)
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (bwa-mem):\n')
	else:
		run_qsub(QUEUE, to_mapp, 4, 'mapp-'+ACC_ID, "4G", PREFIX)
	sys.stdout.write('done\n')
		
	# checking step
	for n in LIB_DIC:
		if not(os.path.isfile(ACC_ID+'/'+TMP+n+'.sam')):
			return 'An error was encountered in step a during mapping step. The program exited without finishing on '+ACC_ID+'\n'
	sys.stdout.write('Mapping step completed.\n')
	
	
	#############
	#1bis Calculating statistics on mapping
	mapping_stat = []
	for n in LIB_DIC:
		stats = '%s stats -r %s %s > %s' % (SAMTOOLS, REF, ACC_ID+'/'+TMP+n+'.sam', ACC_ID+'/'+n+'.sam.stat')
		mapping_stat.append(stats)
	
	if QUEUE == None:
		for comd in mapping_stat:
			sys.stdout.write(comd)
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (bwa-mem):\n')
	else:
		run_qsub(QUEUE, mapping_stat, 1, 'stats-'+ACC_ID, "4G", PREFIX)
	
	# doing plots
	PLOTBAMSTAT = '/'.join(SAMTOOLS.split('/')[0:-1]+['plot-bamstats'])
	mapping_stat = []
	for n in LIB_DIC:
		stats = '%s -p %s %s' % (PLOTBAMSTAT, ACC_ID+'/STATS/', ACC_ID+'/'+n+'.sam.stat')
		mapping_stat.append(stats)
	
	if QUEUE == None:
		for comd in mapping_stat:
			sys.stdout.write(comd)
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (bwa-mem):\n')
	else:
		run_qsub(QUEUE, mapping_stat, 1, 'stats-'+ACC_ID, "4G", PREFIX)
		
	sys.stdout.write('done\n')
	
	
	#############
	#2 Removing unmapped reads
	sys.stdout.write('Initiating unmapped removal.\n')
	to_filter = []
	for n in LIB_DIC:
		unMapRm = '%s view -bS -uF 4 %s | %s view -b -uF 256 - > %s' % (SAMTOOLS, ACC_ID+'/'+TMP+n+'.sam', SAMTOOLS, ACC_ID+'/'+TMP+n+'.bam')
		to_filter.append(unMapRm)
	
	if QUEUE == None:
		for comd in to_filter:
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (samtools):\n')
	else:
		run_qsub(QUEUE, to_filter, 2, 'filter-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')
	
	#############
	#3 Merging
	sys.stdout.write('Initiating merging step.\n')
	if not(os.path.isdir(ACC_ID+'/'+TMP)):
		os.mkdir(ACC_ID+'/'+TMP)
	
	merge = ''
	to_merge = []
	for n in LIB_DIC:
		merge = merge + 'INPUT='+ACC_ID+'/'+TMP+n+'.bam '
	MERGE = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s MergeSamFiles %s OUTPUT=%s MERGE_SEQUENCE_DICTIONARIES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=%s' % (JAVA, PICARD, merge, ACC_ID+'/'+ACC_ID+'_merged.bam', ACC_ID+'/'+TMP)
	to_merge.append(MERGE)
	
	if QUEUE == None:
		for comd in to_merge:
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (picard tools):\n')
	else:
		run_qsub(QUEUE, to_merge, 1, 'merge-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')
	
	# checking step
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_merged.bam')):
		return 'An error was encountered in step a during merging step. The program exited without finishing on '+ACC_ID+'\n'
	
	# Removing intermediate files
	for n in LIB_DIC:
		os.remove(ACC_ID+'/'+TMP+n+'.sam')
		os.remove(ACC_ID+'/'+TMP+n+'.bam')
	
	shutil.rmtree(ACC_ID+'/'+TMP)
	sys.stdout.write('Merging step completed.\n')
	sys.stdout.write('Step a done for accession'+ACC_ID+'.\n')
	
	return 0

def run_step_B(JAVA, PICARD, ACC_ID, TMP, PREFIX, QUEUE):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	############################
	#3 removing duplicates
	sys.stdout.write('Initiating removing duplicates step.\n')
	to_rmdup = []
	rmdup = '%s -XX:ParallelGCThreads=1 -Xmx12G -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true QUIET=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VERBOSITY=WARNING TMP_DIR=%s CREATE_INDEX=true' % (JAVA, PICARD, ACC_ID+'/'+ACC_ID+'_merged.bam', ACC_ID+'/'+ACC_ID+'_rmdup.bam', ACC_ID+'/'+ACC_ID+'_duplicate', ACC_ID+'/'+TMP)
	to_rmdup.append(rmdup)
	
	if QUEUE == None:
		for comd in to_rmdup:
			run_job(getframeinfo(currentframe()), comd, 'Error in step B (picard tools):\n')
	else:
		run_qsub(QUEUE, to_rmdup, 1, 'rmdup-'+ACC_ID, "20G", PREFIX)
	
	# checking step
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_rmdup.bam')):
		return 'An error was encountered in step b during duplicate removal step. The program exited without finishing on '+ACC_ID+'\n'
			
	shutil.rmtree(ACC_ID+'/'+TMP)
	sys.stdout.write('Removing duplicates done\n')
	sys.stdout.write('Step b done for accession'+ACC_ID+'.\n')
	
	return 0

def run_step_C(ACC_ID, JAVA, GATK, REF, CONFIG, PREFIX, QUEUE):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	##########################
	#4 indel realignment
	sys.stdout.write('Initiating indel realignment step.\n')
	known_indels = False
	sys.stdout.write('Looking for a vcf file of known indels...\n')
	if CONFIG.has_section('Variant'):
		if CONFIG.has_option('Variant', 'indel'):
			known_indels = True
		else:
			sys.stdout.write('Not found: The program will do without.\n')
	else:
		sys.stdout.write('Not found: The program will do without.\n')
	
	
	# Looking for available bam files
	rmdup = False
	if os.path.isfile(ACC_ID+'/'+ACC_ID+'_rmdup.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_rmdup.bam'
		rmdup = True
	elif os.path.isfile(ACC_ID+'/'+ACC_ID+'_merged.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_merged.bam'
	else:
		return 'No bam files found in the good name format for '+ACC_ID+'... Please for all accession bam files should be named accession_merged.bam or accession_rmdup.bam or accession_realigned.bam\n'
	
	# Realignment
	if known_indels:
		realTC = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T RealignerTargetCreator -o %s -I %s -R %s -known %s' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'_RTC.intervals', bam, REF, config.get('Variant', 'indel'))
	else:
		realTC = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T RealignerTargetCreator -o %s -I %s -R %s' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'_RTC.intervals', bam, REF)
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), realTC, 'Error in step C (GATK realTC):\n')
	else:
		run_qsub(QUEUE, [realTC], 1, 'realTC-'+ACC_ID, "12G", PREFIX)
	
	# checking step
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_RTC.intervals')):
		return 'An error was encountered in step c during indel realignment step: File '+ACC_ID+'/'+ACC_ID+'_RTC.intervals'+' as not been found.\n'
	
	# IndelRealigner
	Ireal = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T IndelRealigner -o %s -I %s -targetIntervals %s -R %s' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'_realigned.bam', bam, ACC_ID+'/'+ACC_ID+'_RTC.intervals', REF)
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), Ireal, 'Error in step C (GATK Ireal):\n')
	else:
		run_qsub(QUEUE, [Ireal], 1, 'Ireal-'+ACC_ID, "12G", PREFIX)
	
	# checking step
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam')):
		return 'An error was encountered in step c during indel realignment step: File '+ACC_ID+'/'+ACC_ID+'_realigned.bam'+' as not been found.\n'
			
	# removing file:
	os.remove(ACC_ID+'/'+ACC_ID+'_RTC.intervals')
	sys.stdout.write('Indel realignment done.\n')
	sys.stdout.write('Step c done for accession'+ACC_ID+'.\n')
	
	return 0

def run_step_D(CONFIG, ACC_ID, UseUnifiedGenotyperForBaseRecal, JAVA, GATK, REF, PLOIDY, PREFIX, QUEUE):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	#########################
	#5 Base recalibration
	HCopt = ''
	known_indels = False
	known_snp = False
	sys.stdout.write('Looking for a vcf file of known variant...\n')
	if CONFIG.has_section('Variant'):
		if CONFIG.has_option('Variant', 'HCopt'):
			HCopt = CONFIG.get('Variant', 'HCopt')
		sys.stdout.write('Looking for a vcf file of known indels...\n')
		if CONFIG.has_option('Variant', 'indel'):
			known_indels = CONFIG.get('Variant', 'indel')
		else:
			sys.stdout.write('Not found: The program will do without but calculation time will be much longer...\n')
		
		sys.stdout.write('Looking for a vcf file of known snps...\n')
		if CONFIG.has_option('Variant', 'snp'):
			known_snp = CONFIG.get('Variant', 'snp')
		else:
			sys.stdout.write('Not found: The program will do without but calculation time will be much longer...\n')
	else:
		sys.stdout.write('Not found: The program will do without but calculation time will be much longer...\n')
	
	if os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_realigned.bam'
	else:
		return 'An error was encountered in step d during indel realignment step: File '+ACC_ID+'/'+ACC_ID+'_realigned.bam'+' as not been found.\n'
	
	if known_indels == False & known_snp == False:
		
		sys.stdout.write('Initiating first variant calling. Needed because no variant files are provided to do the base re-calibration...\n')
		# UnifiedGenotyper
		if UseUnifiedGenotyperForBaseRecal == "yes":
			hapcall = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T UnifiedGenotyper -I %s -o %s -R %s -stand_call_conf 100 -stand_emit_conf 100 --output_mode EMIT_VARIANTS_ONLY --genotyping_mode DISCOVERY --sample_ploidy %s' % (JAVA, GATK, bam, ACC_ID+'/'+ACC_ID+'.vcf', REF, PLOIDY)
		else:
			hapcall = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T HaplotypeCaller -I %s -o %s -R %s -dontUseSoftClippedBases -stand_call_conf 100 -stand_emit_conf 100 --output_mode EMIT_VARIANTS_ONLY --genotyping_mode DISCOVERY --sample_ploidy %s' % (JAVA, GATK, bam, ACC_ID+'/'+ACC_ID+'.vcf', REF, PLOIDY)
		
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), hapcall, 'Error in step D (GATK HCreal):\n')
		else:
			run_qsub(QUEUE, [hapcall], 1, 'HCreal-'+ACC_ID, "12G", PREFIX)
		sys.stdout.write('done\n')
		
		sys.stdout.write('Variant filtration...\n')
		# VariantFiltration
		varflt = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T VariantFiltration -V %s --out %s -R %s --clusterSize 3 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterExpression "QD < 1.5" --filterExpression "DP < 20" --maskExtension 0 --filterName HARD_TO_VALIDATE --filterName QD_FILTER --filterName DP_FILTER --maskName Mask --clusterWindowSize 10' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'.vcf', ACC_ID+'/'+ACC_ID+'_filtered.vcf', REF)
		
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), varflt, 'Error in step D (GATK varflt):\n')
		else:
			run_qsub(QUEUE, [varflt], 1, 'VFreal-'+ACC_ID, "12G", PREFIX)
		sys.stdout.write('done\n')
		
		sys.stdout.write('Variant selection...\n')
		file_known_snp = ""
		# SelectVariants
		selvar = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T SelectVariants -V %s -o %s -R %s -select "vc.isNotFiltered()" -selectType SNP' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'_filtered.vcf', ACC_ID+'/'+ACC_ID+'_selected.vcf', REF)
		file_known_snp = ACC_ID+'/'+ACC_ID+'_selected.vcf'
		
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), selvar, 'Error in step D (GATK selvar):\n')
		else:
			run_qsub(QUEUE, [selvar], 1, 'SVreal-'+ACC_ID, "12G", PREFIX)
		sys.stdout.write('done\n')
		
	sys.stdout.write('Initiating base recalibration...\n')
	# BaseRecalibrator
	if known_indels == False & known_snp == False:
		knownsites = file_known_snp
	elif known_indels == False:
		knownsites = known_snp
	elif known_snp == False:
		knownsites = known_indels
	else:
		knownsites = known_snp+'-knownSites '+known_indels
	baserecal = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T BaseRecalibrator -I %s -knownSites %s -o %s -R %s' % (JAVA, GATK, bam, knownsites, ACC_ID+'/'+ACC_ID+'_BR_recal.grp', REF)
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), baserecal, 'Error in step D (GATK baserecal 1):\n')
	else:
		run_qsub(QUEUE, [baserecal], 1, 'BR1real-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')
	
	sys.stdout.write('Initiating analysis of covariation remaining after recalibration...\n')
	# BaseRecalibrator
	if known_indels == False & known_snp == False:
		knownsites = file_known_snp
	elif known_indels == False:
		knownsites = known_snp
	elif known_snp == False:
		knownsites = known_indels
	else:
		knownsites = known_snp+'-knownSites '+known_indels
	baserecal = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T BaseRecalibrator -I %s -knownSites %s -BQSR %s -R %s -o %s' % (JAVA, GATK, bam, knownsites, ACC_ID+'/'+ACC_ID+'_BR_recal.grp', REF, ACC_ID+'/'+ACC_ID+'_BR_post_recal.grp')
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), baserecal, 'Error in step D (GATK baserecal 2):\n')
	else:
		run_qsub(QUEUE, [baserecal], 1, 'BR2real-'+ACC_ID, "12G", PREFIX)
	
	sys.stdout.write('Generating before/after plot...\n')
	# AnalyzeCovariates
	anacov = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T AnalyzeCovariates -R %s -before %s -after %s -plots %s' % (JAVA, GATK, REF, ACC_ID+'/'+ACC_ID+'_BR_recal.grp', ACC_ID+'/'+ACC_ID+'_BR_post_recal.grp', ACC_ID+'/'+ACC_ID+'_recalibration_plot.pdf')
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), anacov, 'Error in step D (GATK anacov):\n')
	else:
		run_qsub(QUEUE, [anacov], 1, 'ACreal-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')
	
	sys.stdout.write('Applying recalibration...\n')
	# PrintReads
	printread = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T PrintReads -BQSR %s -I %s -o %s -R %s' % (JAVA, GATK, ACC_ID+'/'+ACC_ID+'_BR_recal.grp', bam, ACC_ID+'/'+ACC_ID+'_real_recal.bam', REF)
	
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), printread, 'Error in step D (GATK printread):\n')
	else:
		run_qsub(QUEUE, [printread], 1, 'PRreal-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')

	# checking step
	sys.stdout.write('Step d final checking...\n')
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_real_recal.bam')):
		return 'An error was encountered in step d during recalibration step: File '+ACC_ID+'/'+ACC_ID+'_real_recal.bam'+' as not been found.\n'
	sys.stdout.write('OK\n')
	
	# removing file:
	os.remove(ACC_ID+'/'+ACC_ID+'_BR_recal.grp')
	os.remove(ACC_ID+'/'+ACC_ID+'_BR_post_recal.grp')
	if known_indels == False & known_snp == False:
		os.remove(ACC_ID+'/'+ACC_ID+'.vcf')
		os.remove(ACC_ID+'/'+ACC_ID+'_filtered.vcf')
		os.remove(ACC_ID+'/'+ACC_ID+'.vcf.idx')
		os.remove(ACC_ID+'/'+ACC_ID+'_filtered.vcf.idx')
	sys.stdout.write('Step d done for accession'+ACC_ID+'.\n')
	
	return 0

def run_step_E(ACC_ID, PYTHON, GATK, REF, PLOIDY, CONFIG, DICO_CHR, PREFIX, QUEUE, PATHNAME):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	#######################
	#6 pseudo GVCF generation
	
	if os.path.isfile(ACC_ID+'/'+ACC_ID+'_real_recal.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_real_recal.bam'
	elif os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_realigned.bam'
		sys.stdout.write('Warning no recalibrated bam file found for '+ACC_ID+'... Continue using realigned bam instead.\n')
		sys.stdout.flush()
	else:
		return 'An error was encountered in step e during recalibration step: File '+ACC_ID+'/'+ACC_ID+'_real_recal.bam'+' as not been found.\n'
	
	# Running allele count in bam
	sys.stdout.write('Initiating allele count\n')
	allele_count = '%s %s/count_allele_number.py -r %s -b %s -o %s' % (PYTHON, PATHNAME, REF, bam, ACC_ID+'/'+ACC_ID+'_allele_count')
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), allele_count, 'Error in step E:\n')
	else:
		run_qsub(QUEUE, [allele_count], 1, 'AlCount-'+ACC_ID, "12G", PREFIX)
	
	# checking step
	for n in DICO_CHR:
		if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_allele_count_'+n+'.gz')):
			sys.exit('ERROR : File '+ACC_ID+'/'+ACC_ID+'_allele_count_'+n+'.gz'+' as not been found.\n')
	sys.stdout.write('Step e done for accession '+ACC_ID+'.\n')
	
	return 0

def combin(n, k):
	"""Nombre de combinaisons de n objets pris k a k"""
	if k > n//2:
		k = n-k
	x = 1
	y = 1
	i = n-k+1
	while i <= n:
		x = (x*i)//y
		y += 1
		i += 1
	return x

def binom(k,n,p):
	x = combin(n,k)*(p**k)*((1-p)**(n-k))
	return x

def genotype_accession(COVERAGE, ALLELE, ERROR, PLOIDY):
	
	# WARNING The statistics will not be accurate on low coverage
	
	# To manage huge factorial
	nb_tirage = sum(COVERAGE)
	if nb_tirage >= 1000:
		coverage = list(map(int, [x/float(nb_tirage)*100 for x in COVERAGE]))
	else:
		coverage = list(COVERAGE)
	nb_tirage = sum(coverage)
	
	dico_proba = {}
	# calculating probabilities for homozygous state
	p = 1-ERROR
	for n in range(len(ALLELE)):
		dico_proba['/'.join([str(n)]*int(PLOIDY))] = binom(coverage[n],nb_tirage,p)/binom(int(round(nb_tirage*p)),nb_tirage,p)
	
		if PLOIDY == '2':
			# calculating probabilities for heterozygous state
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n),str(k)])] = (binom(coverage[n], nb_tirage, 0.5)*binom(coverage[k], nb_tirage, 0.5))/(binom(int(round(nb_tirage*0.5)),nb_tirage,0.5)*binom(int(round(nb_tirage*0.5)),nb_tirage,0.5))
							comb_done.append(couple)
	
		if PLOIDY == '3':
			# calculating probabilities for heterozygous state for triploids
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 2.0/3.0)*binom(coverage[n], nb_tirage, 1.0/3.0))/(binom(int(round(nb_tirage*(2.0/3.0))),nb_tirage,(2.0/3.0))*binom(int(round(nb_tirage*(1.0/3.0))),nb_tirage,(1.0/3.0)))
							comb_done.append(couple)
						couple = sorted([n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 2.0/3.0)*binom(coverage[k], nb_tirage, 1.0/3.0))/(binom(int(round(nb_tirage*(2.0/3.0))),nb_tirage,(2.0/3.0))*binom(int(round(nb_tirage*(1.0/3.0))),nb_tirage,(1.0/3.0)))
							comb_done.append(couple)
		
		if PLOIDY == '4':
			# calculating probabilities for heterozygous state for tetraploid
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 3.0/4.0)*binom(coverage[n], nb_tirage, 1.0/4.0))/(binom(int(round(nb_tirage*(3.0/4.0))),nb_tirage,(3.0/4.0))*binom(int(round(nb_tirage*(1.0/4.0))),nb_tirage,(1.0/4.0)))
							comb_done.append(couple)
						couple = sorted([n,n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 2.0/4.0)*binom(coverage[k], nb_tirage, 2.0/4.0))/(binom(int(round(nb_tirage*(2.0/4.0))),nb_tirage,(2.0/4.0))*binom(int(round(nb_tirage*(2.0/4.0))),nb_tirage,(2.0/4.0)))
							comb_done.append(couple)
						couple = sorted([n,n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 3.0/4.0)*binom(coverage[k], nb_tirage, 1.0/4.0))/(binom(int(round(nb_tirage*(3.0/4.0))),nb_tirage,(3.0/4.0))*binom(int(round(nb_tirage*(1.0/4.0))),nb_tirage,(1.0/4.0)))
							comb_done.append(couple)
		
		if PLOIDY == '12':
			# calculating probabilities for heterozygous state for tetraploid
			set_done = set()
			comb_done = []
			for n in range(len(ALLELE)):
				set_done.add(n)
				for k in range(len(ALLELE)):
					if not(k in set_done):
						couple = sorted([n,k,k,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 11.0/12.0)*binom(coverage[n], nb_tirage, 1.0/12.0))/(binom(int(round(nb_tirage*(11.0/12.0))),nb_tirage,(11.0/12.0))*binom(int(round(nb_tirage*(1.0/12.0))),nb_tirage,(1.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,k,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 10.0/12.0)*binom(coverage[n], nb_tirage, 2.0/12.0))/(binom(int(round(nb_tirage*(10.0/12.0))),nb_tirage,(10.0/12.0))*binom(int(round(nb_tirage*(2.0/12.0))),nb_tirage,(2.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,k,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 9.0/12.0)*binom(coverage[n], nb_tirage, 3.0/12.0))/(binom(int(round(nb_tirage*(9.0/12.0))),nb_tirage,(9.0/12.0))*binom(int(round(nb_tirage*(3.0/12.0))),nb_tirage,(3.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,k,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 8.0/12.0)*binom(coverage[n], nb_tirage, 4.0/12.0))/(binom(int(round(nb_tirage*(8.0/12.0))),nb_tirage,(8.0/12.0))*binom(int(round(nb_tirage*(4.0/12.0))),nb_tirage,(4.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,k,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 7.0/12.0)*binom(coverage[n], nb_tirage, 5.0/12.0))/(binom(int(round(nb_tirage*(7.0/12.0))),nb_tirage,(7.0/12.0))*binom(int(round(nb_tirage*(5.0/12.0))),nb_tirage,(5.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,k,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[k], nb_tirage, 6.0/12.0)*binom(coverage[n], nb_tirage, 6.0/12.0))/(binom(int(round(nb_tirage*(6.0/12.0))),nb_tirage,(6.0/12.0))*binom(int(round(nb_tirage*(6.0/12.0))),nb_tirage,(6.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,k,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 7.0/12.0)*binom(coverage[k], nb_tirage, 5.0/12.0))/(binom(int(round(nb_tirage*(7.0/12.0))),nb_tirage,(7.0/12.0))*binom(int(round(nb_tirage*(5.0/12.0))),nb_tirage,(5.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,k,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 8.0/12.0)*binom(coverage[k], nb_tirage, 4.0/12.0))/(binom(int(round(nb_tirage*(8.0/12.0))),nb_tirage,(8.0/12.0))*binom(int(round(nb_tirage*(4.0/12.0))),nb_tirage,(4.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,k,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 9.0/12.0)*binom(coverage[k], nb_tirage, 3.0/12.0))/(binom(int(round(nb_tirage*(9.0/12.0))),nb_tirage,(9.0/12.0))*binom(int(round(nb_tirage*(3.0/12.0))),nb_tirage,(3.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,n,k,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k), str(k)])] = (binom(coverage[n], nb_tirage, 10.0/12.0)*binom(coverage[k], nb_tirage, 2.0/12.0))/(binom(int(round(nb_tirage*(10.0/12.0))),nb_tirage,(10.0/12.0))*binom(int(round(nb_tirage*(2.0/12.0))),nb_tirage,(2.0/12.0)))
							comb_done.append(couple)
							
						couple = sorted([n,n,n,n,n,n,n,n,n,n,n,k])
						if not(couple in comb_done):
							dico_proba['/'.join([str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(n), str(k)])] = (binom(coverage[n], nb_tirage, 11.0/12.0)*binom(coverage[k], nb_tirage, 1.0/12.0))/(binom(int(round(nb_tirage*(11.0/12.0))),nb_tirage,(11.0/12.0))*binom(int(round(nb_tirage*(1.0/12.0))),nb_tirage,(1.0/12.0)))
							comb_done.append(couple)
		
	# getting best probability
	best_genotype = '/'.join(['.']*int(PLOIDY))
	best_value = 0
	for genotype in dico_proba:
		if best_value < dico_proba[genotype]:
			best_value = dico_proba[genotype]
			best_genotype = genotype
		elif best_value == dico_proba[genotype]:
			best_genotype = '/'.join(['.']*int(PLOIDY))
	return best_genotype
	
def get_min_acc_pos(DICO):
	
	min_pos = False
	for n in DICO:
		if min_pos:
			min_pos = min(min_pos, DICO[n]) 
		else:
			min_pos = DICO[n]
	return min_pos

def get_allele_coverage (ALLELE, FORMAT, LISTE):
	
	if not(ALLELE in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '*']):
		return 0
	else:
		return int(LISTE[FORMAT.index(ALLELE)])

def create_pseudo_VCF(LIST_ACC, REF, PREFIX, DICO_PLOIDY, DICO_CHR, CHR, START, END):
	
	#1- preparing output
	outfile = open(PREFIX+'_'+CHR+'_'+str(START)+'_'+str(END)+'_allele_count.vcf', 'w')
	outfile.write("##fileformat=VCFv4.2\n")
	outfile.write("##reference=file:///"+REF+"\n")
	outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	outfile.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
	for n in DICO_CHR:
		outfile.write("##contig=<ID="+n+",length="+str(DICO_CHR[n])+">\n")
	liste2print = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
	for acc in LIST_ACC:
		liste2print.append(acc)
	outfile.write('\t'.join(liste2print))
	outfile.write('\n')
	
	#2- Initiating reading files and variables
	dico_accession_infile = {}
		
	dico_accession_infile[CHR] = {}
	# sys.stdout.write(CHR+'\n')
	for acc in LIST_ACC:
		dico_accession_infile[CHR][acc] = gzip.open(acc+'/'+acc+'_allele_count_'+CHR+'.gz', 'rb')

	#3- Initiating variables
	dico_accession_positions = {}
	dico_accession_line = {}
	for acc in LIST_ACC:
		dico_accession_positions[acc] = 0
		dico_accession_line[acc] = 0
	
	liste_allele = ['A', 'C', 'G', 'T', '*']
		
	#4a- Initiating the loop
	
	#4aa- phasing all files to the expected window
	old_position = 0
	for acc in LIST_ACC:
		# initiating the line
		dico_accession_line[acc] = dico_accession_infile[CHR][acc].readline().split()
		if dico_accession_line[acc] != []:
			while int(dico_accession_line[acc][1]) <= START:
				dico_accession_line[acc] = dico_accession_infile[CHR][acc].readline().split()
				if dico_accession_line[acc] == []:
					break
	
	#4ab- recording position informations
	for acc in LIST_ACC:
		if dico_accession_line[acc] == []:
			dico_accession_positions[acc] = 1000000000000000
		else:
			dico_accession_positions[acc] = int(dico_accession_line[acc][1])
	min_position = get_min_acc_pos(dico_accession_positions)
	
	# All file are not empty
	if min_position != 1000000000000000 and min_position <= END:
	
		#4ac- recording total allele coverage information and reference allele
		reference = False
		dico_allele_cov = {}
		for allele in liste_allele:
			dico_allele_cov[allele] = 0
			
		for acc in LIST_ACC:
			if dico_accession_positions[acc] == min_position:
				reference = dico_accession_line[acc][2] # line to get reference allele
				format = dico_accession_line[acc][4].split(':')
				allele_count = dico_accession_line[acc][5].split(':')
				for allele in dico_allele_cov:
					dico_allele_cov[allele] += get_allele_coverage (allele, format, allele_count)
		
		#4ad- recording expected allele
		allele_to_keep = [reference]
		for allele in dico_allele_cov:
			if allele != reference:
				if dico_allele_cov[allele]:
					allele_to_keep.append(allele)
		# print allele_to_keep
		
		#4ae- printing results to output if necessary
		if len(allele_to_keep) > 1:
			liste2print = [CHR, str(min_position), '.', allele_to_keep[0], ','.join(allele_to_keep[1:]), '.', '.', '.', 'GT:AD:DP']
			for acc in LIST_ACC:
				liste_cov = []
				liste_cov_int = []
				if dico_accession_positions[acc] == min_position:
					format = dico_accession_line[acc][4].split(':')
					allele_count = dico_accession_line[acc][5].split(':')
					for allele in allele_to_keep:
						allele_coverage = get_allele_coverage (allele, format, allele_count)
						liste_cov.append(str(allele_coverage))
						liste_cov_int.append(allele_coverage)
					depth = dico_accession_line[acc][3]
					genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.005, DICO_PLOIDY[acc])
				else:
					for allele in allele_to_keep:
						liste_cov.append('0')
					depth = '0'
					genotype = '/'.join(['.']*int(DICO_PLOIDY[acc]))
				liste2print.append(':'.join([genotype,','.join(liste_cov),depth]))
			# print liste2print
			outfile.write('\t'.join(liste2print))
			outfile.write('\n')
			
	
		#4b- we are in the loop
		while old_position != min_position:
			
			#4ba- Reading new lines
			for acc in LIST_ACC:
				if dico_accession_positions[acc] == min_position:
					dico_accession_line[acc] = dico_accession_infile[CHR][acc].readline().split()
			
			#4bb- recording position informations
			for acc in LIST_ACC:
				if dico_accession_line[acc] == []:
					dico_accession_positions[acc] = 1000000000000000
				else:
					dico_accession_positions[acc] = int(dico_accession_line[acc][1])
				# print acc, dico_accession_line[acc]
			min_position = get_min_acc_pos(dico_accession_positions)
			# print min_position
			
			if min_position == 1000000000000000 or min_position > END:
				break
			
			#4bc- recording total allele coverage information and reference allele
			reference = False
			dico_allele_cov = {}
			for allele in liste_allele:
				dico_allele_cov[allele] = 0
				
			for acc in LIST_ACC:
				if dico_accession_positions[acc] == min_position:
					reference = dico_accession_line[acc][2] # line to get reference allele
					format = dico_accession_line[acc][4].split(':')
					allele_count = dico_accession_line[acc][5].split(':')
					for allele in dico_allele_cov:
						dico_allele_cov[allele] += get_allele_coverage (allele, format, allele_count)
			# print dico_allele_cov
			
			#4bd- recording expected allele
			allele_to_keep = [reference]
			for allele in dico_allele_cov:
				if allele != reference:
					if dico_allele_cov[allele]:
						allele_to_keep.append(allele)
			# print allele_to_keep
			
			#4be- printing results to output if necessary
			if len(allele_to_keep) > 1:
				liste2print = [CHR, str(min_position), '.', allele_to_keep[0], ','.join(allele_to_keep[1:]), '.', '.', '.', 'GT:AD:DP']
				for acc in LIST_ACC:
					liste_cov = []
					liste_cov_int = []
					if dico_accession_positions[acc] == min_position:
						format = dico_accession_line[acc][4].split(':')
						allele_count = dico_accession_line[acc][5].split(':')
						for allele in allele_to_keep:
							allele_coverage = get_allele_coverage (allele, format, allele_count)
							liste_cov.append(str(allele_coverage))
							liste_cov_int.append(allele_coverage)
						depth = dico_accession_line[acc][3]
						genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.005, DICO_PLOIDY[acc])
					else:
						for allele in allele_to_keep:
							liste_cov.append('0')
						depth = '0'
						genotype = '/'.join(['.']*int(DICO_PLOIDY[acc]))
					liste2print.append(':'.join([genotype,','.join(liste_cov),depth]))
				# print liste2print
				outfile.write('\t'.join(liste2print))
				outfile.write('\n')
	
	#5- Closing files
	for acc in LIST_ACC:
		dico_accession_infile[CHR][acc].close()
	return 0

def merge_vcf(PREFIX, DICO_CHR):
	
	# recording VCF headers
	dico_header = {}
	accessions = set()
	for chr in DICO_CHR:
		file = open(PREFIX+'_'+chr+'_all_allele_count.vcf','r')
		if chr in dico_header:
			return 1
		if not ('header' in dico_header):
			dico_header['header'] = []
			for line in file:
				data = line.split()
				if line[0:2] == "##":
					dico_header['header'].append(line)
				elif data[0] == "#CHROM":
					dico_header[chr] = data
					for acc in data[9:]:
						accessions.add(acc)
					break
		else:
			for line in file:
				data = line.split()
				if data[0] == "#CHROM":
					dico_header[chr] = data
					for acc in data[9:]:
						accessions.add(acc)
					break
		file.close()
	
	# Sorting accessions
	accessions = sorted(list(accessions))
	
	# recording chromosome order
	chromosome_order = []
	for line in dico_header['header']:
		data = line.split()
		if "##contig" in data[0]:
			chromosome_order.append(data[0].split(',')[0].split('=')[2])
	
	# Printing the final VCF
	outfile = open(PREFIX+'_all_allele_count.vcf','w')
	outfile.write(''.join(dico_header['header']))
	outfile.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+accessions)+'\n')
	for chr in chromosome_order:
		file = open(PREFIX+'_'+chr+'_all_allele_count.vcf','r')
		for line in file:
			data = line.split()
			if data[0][0] != "#":
				liste2print = data[0:9]
				for acc in accessions:
					liste2print.append(data[dico_header[chr].index(acc)])
				outfile.write('\t'.join(liste2print)+'\n')
		file.close()
	outfile.close()
	return 0

def merge_sub_vcf(PREFIX, CHR, LIST):
	
	# recording VCF headers
	dico_header = {}
	accessions = []
	for pos in LIST:
		POS = '-'.join(map(str,pos))
		file = open(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf','r')
		if POS in dico_header:
			return 1
		if not ('header' in dico_header):
			dico_header['header'] = []
			for line in file:
				data = line.split()
				if line[0:2] == "##":
					dico_header['header'].append(line)
				elif data[0] == "#CHROM":
					dico_header[POS] = data
					for acc in data[9:]:
						accessions.append(acc)
					break
		else:
			for line in file:
				data = line.split()
				if data[0] == "#CHROM":
					dico_header[POS] = data
					break
		file.close()
	
	# recording chromosome order
	chromosome_order = []
	for line in dico_header['header']:
		data = line.split()
		if "##contig" in data[0]:
			chromosome_order.append(data[0].split(',')[0].split('=')[2])
	
	# Printing the final VCF
	outfile = open(PREFIX+'_'+CHR+'_all_allele_count.vcf','w')
	outfile.write(''.join(dico_header['header']))
	outfile.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+accessions)+'\n')
	for pos in LIST:
		POS = '-'.join(map(str,pos))
		file = open(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf','r')
		for line in file:
			data = line.split()
			if data[0][0] != "#":
				liste2print = data[0:9]
				for acc in accessions:
					liste2print.append(data[dico_header[POS].index(acc)])
				outfile.write('\t'.join(liste2print)+'\n')
		file.close()
	outfile.close()
	
	for pos in LIST:
		os.remove(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf')
	return 0

def perform_count(BAMREADCOUNT, REF, BAM, OUT):
	
	bam_count = '%s -b 10 -q 1 -d 1000000 -w 0 -f %s %s > %s' % (BAMREADCOUNT, REF, BAM, OUT)
	os.system(bam_count)

def generate_pseudo_vcf_oldies(TAB, REF, OUT):
	
	# calculating sequence length
	dico_chr = {}
	sequence_dict = SeqIO.index(REF, "fasta")
	for n in sequence_dict:
		dico_chr[n] = [len(str(sequence_dict[n].seq)), str(sequence_dict[n].seq)]
	del sequence_dict
	
	# Verification
	file = open(TAB)
	i = 0
	chr = ""
	dico_allele = {}
	liste_allele = ['A', 'C', 'G', 'T']
	mot_allele = ':'.join(liste_allele)
	for line in file:
		data = line.split()
		i += 1
		if data:
			if chr != data[0]:
				if not(chr == ""):
					outfile.close()
				chr = data[0]
				outfile = gzip.open(OUT+'_'+chr+'.gz', 'wb')
			# recording information on the current line 
			for n in data[4:]:
				allele_stat = n.split(':')
				dico_allele[allele_stat[0]] = allele_stat[1]
			# preparing and printing result
			liste_mot = []
			Reference_allele = data[2].replace('a', 'A').replace('c', 'C').replace('g', 'G').replace('t', 'T').replace('n', 'N')
			for n in liste_allele:
				liste_mot.append(dico_allele[n])
			outfile.write('\t' . join([data[0], data[1], Reference_allele, data[3], mot_allele, ':'.join(liste_mot)]))
			outfile.write('\n')
	outfile.close()
	os.remove(TAB)
	
	for n in dico_chr:
		if not(os.path.isfile(OUT+'_'+n+'.gz')):
			outfile = gzip.open(OUT+'_'+n+'.gz', 'wb')
			outfile.close()

def generate_pseudo_vcf(TAB, REF, OUT):
	
	# calculating sequence length
	sys.stdout.write("Loading reference...")
	sys.stdout.flush()
	dico_chr = {}
	sequence_dict = SeqIO.index(REF, "fasta")
	for n in sequence_dict:
		dico_chr[n] = [len(str(sequence_dict[n].seq)), str(sequence_dict[n].seq)]
	del sequence_dict
	sys.stdout.write("done\n")
	sys.stdout.flush()
	
	# Verification
	file = open(TAB)
	i = 0
	chr = ""
	dico_allele = {}
	dico_deletions = {}
	for line in file:
		liste_allele = []
		data = line.split()
		pos = int(data[1])
		Final_cov = int(data[3])
		sum_deletion = 0
		i += 1
		if data:
			if chr != data[0]:
				if not(chr == ""):
					outfile.close()
				chr = data[0]
				outfile = gzip.open(OUT+'_'+chr+'.gz', 'wb')
			# Recording information on the current line 
			for n in data[4:]:
				if n[0] != "=": # To remove the column begining with "=" which I don't know what it means...
					allele_stat = n.split(':')
					dico_allele[allele_stat[0]] = allele_stat[1]
					liste_allele.append(allele_stat[0])
					# To manage deletions
					if allele_stat[0][0] == '-':
						cov = int(allele_stat[1]) # Getting number of reads with the deletion
						for i in range(len(allele_stat[0][1:])): # Recording this number for adjacent sites with this same deletion (for final coverage calculation)
							if not(pos+i in dico_deletions):
								dico_deletions[pos+i] = 0
							dico_deletions[pos+i] += cov
						Final_cov -= cov
						sum_deletion += cov
			# Preparing and printing result
			liste_mot = []
			for n in liste_allele:
				liste_mot.append(dico_allele[n])
			
			Reference_allele = data[2].replace('a', 'A').replace('c', 'C').replace('g', 'G').replace('t', 'T').replace('n', 'N')
			if pos in dico_deletions:
				if dico_deletions[pos] != sum_deletion:
					liste_allele.append("*")
					liste_mot.append(str(dico_deletions[pos]-sum_deletion))
				else:
					liste_allele.append("*")
					liste_mot.append(str(dico_deletions[pos]))
				outfile.write('\t' . join([data[0], data[1], Reference_allele, str(Final_cov+dico_deletions[pos]), ':'.join(liste_allele), ':'.join(liste_mot)]))
				del dico_deletions[pos]
			else:
				liste_allele.append("*")
				liste_mot.append(str(0))
				outfile.write('\t' . join([data[0], data[1], Reference_allele, data[3], ':'.join(liste_allele), ':'.join(liste_mot)]))
			outfile.write('\n')
	outfile.close()
	os.remove(TAB)
	
	for n in dico_chr:
		if not(os.path.isfile(OUT+'_'+n+'.gz')):
			outfile = gzip.open(OUT+'_'+n+'.gz', 'wb')
			outfile.close()

def Calc_stats(PREFIX, DICO_LIB):
	
	outfile1 = open(PREFIX+'_lib.stats','w')
	outfile2 = open(PREFIX+'_acc.stats','w')
	
	outfile1.write('\t'.join(['Library', 'Accession', 'RawReadNumber', 'ReadsPaired', 'ReadsSingle', 'ReadsMapped', 'ReadsMapped%', 'ReadsUnmapped', 'ReadsUnmapped%', 'ReadMappedPaired', 'ReadMappedPaired(%Pairs)', 'ReadsMappedSingle', 'ReadsMappedPairedProperPair', 'ReadsMappedPairedProperPair(%MappedPairs)', 'MultipleHitsReads', 'MultipleHitsReads(%Mapped)']))
	outfile1.write('\n')
	
	
	outfile2.write('\t'.join(['Accession', 'RawReadNumber', 'ReadsPaired', 'ReadsSingle', 'ReadsMapped', 'ReadsMapped%', 'ReadsUnmapped', 'ReadsUnmapped%', 'ReadMappedPaired', 'ReadMappedPaired(%Pairs)', 'ReadsMappedSingle', 'ReadsMappedPairedProperPair', 'ReadsMappedPairedProperPair(%MappedPairs)', 'MultipleHitsReads', 'MultipleHitsReads(%Mapped)']))
	outfile2.write('\n')
	
	dicoInfoGlob = {}
	dicoInfoLibs = {}
	for n in DICO_LIB:
		dicoInfoLibs[n] = {}
		dicoInfoGlob[n] = {}
		dicoInfoGlob[n]['total'] = 0
		dicoInfoGlob[n]['paired'] = 0
		dicoInfoGlob[n]['single'] = 0
		dicoInfoGlob[n]['mapped'] = 0
		dicoInfoGlob[n]['unmapped'] = 0
		dicoInfoGlob[n]['PairedMapped'] = 0
		dicoInfoGlob[n]['MultipleHits'] = 0
		dicoInfoGlob[n]['PoperlyPaired'] = 0
		listLibs = list(DICO_LIB[n].keys())
		for lib in listLibs:
			dicoInfoLibs[n][lib] = {}
			if os.path.isfile('/'.join([n,lib+'.sam.stat'])):
				file = open('/'.join([n,lib+'.sam.stat']))
				for line in file:
					data = line.split()
					if data:
						if data[0] == 'FFQ':
							break
						if ' '.join(data[1:4]) == 'raw total sequences:' and len(data) == 5:
							dicoInfoLibs[n][lib]['total'] = int(data[4])
							dicoInfoGlob[n]['total'] += int(data[4])
						elif ' '.join(data[1:3]) == 'reads paired:' and len(data) == 9:
							dicoInfoLibs[n][lib]['paired'] = int(data[3])
							dicoInfoLibs[n][lib]['single'] = dicoInfoLibs[n][lib]['total'] - int(data[3])
							dicoInfoGlob[n]['paired'] += dicoInfoLibs[n][lib]['paired']
							dicoInfoGlob[n]['single'] += dicoInfoLibs[n][lib]['single']
						elif ' '.join(data[1:3]) == 'reads mapped:' and len(data) == 4:
							dicoInfoLibs[n][lib]['mapped'] = int(data[3])
							dicoInfoGlob[n]['mapped'] += dicoInfoLibs[n][lib]['mapped']
						elif ' '.join(data[1:5]) == 'reads mapped and paired:' and len(data) == 15:
							dicoInfoLibs[n][lib]['PairedMapped'] = int(data[5])
							dicoInfoGlob[n]['PairedMapped'] += dicoInfoLibs[n][lib]['PairedMapped']
						elif ' '.join(data[1:3]) == 'reads MQ0:' and len(data) == 8:
							dicoInfoLibs[n][lib]['MultipleHits'] = int(data[3])
							dicoInfoGlob[n]['MultipleHits'] += dicoInfoLibs[n][lib]['MultipleHits']
						elif ' '.join(data[1:4]) == 'reads properly paired:' and len(data) == 9:
							dicoInfoLibs[n][lib]['PoperlyPaired'] = int(data[4])
							dicoInfoGlob[n]['PoperlyPaired'] += dicoInfoLibs[n][lib]['PoperlyPaired']
						elif ' '.join(data[1:3]) == 'reads unmapped:' and len(data) == 4:
							dicoInfoLibs[n][lib]['unmapped'] = int(data[3])
							dicoInfoGlob[n]['unmapped'] += dicoInfoLibs[n][lib]['unmapped']
				file.close()
				if dicoInfoLibs[n][lib]['total'] == 0:
					outfile1.write('\t'.join(list(map(str,[lib, n, dicoInfoLibs[n][lib]['total'], dicoInfoLibs[n][lib]['paired'], dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['mapped'], 0, dicoInfoLibs[n][lib]['unmapped'], 0, dicoInfoLibs[n][lib]['PairedMapped'], 0, dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['PoperlyPaired'], 0, dicoInfoLibs[n][lib]['MultipleHits'], 0]))))
				elif dicoInfoLibs[n][lib]['mapped'] == 0:
					outfile1.write('\t'.join(list(map(str,[lib, n, dicoInfoLibs[n][lib]['total'], dicoInfoLibs[n][lib]['paired'], dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['mapped'], dicoInfoLibs[n][lib]['mapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['unmapped'], dicoInfoLibs[n][lib]['unmapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['PairedMapped'], dicoInfoLibs[n][lib]['PairedMapped']/float(dicoInfoLibs[n][lib]['paired'])*100, dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['PoperlyPaired'], 0, dicoInfoLibs[n][lib]['MultipleHits'], 0]))))
				elif dicoInfoLibs[n][lib]['paired'] == 0:
					outfile1.write('\t'.join(list(map(str,[lib, n, dicoInfoLibs[n][lib]['total'], dicoInfoLibs[n][lib]['paired'], dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['mapped'], dicoInfoLibs[n][lib]['mapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['unmapped'], dicoInfoLibs[n][lib]['unmapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['PairedMapped'], 0, dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['PoperlyPaired'], 0, dicoInfoLibs[n][lib]['MultipleHits'], dicoInfoLibs[n][lib]['MultipleHits']/float(dicoInfoLibs[n][lib]['mapped'])*100]))))
				elif dicoInfoLibs[n][lib]['PairedMapped'] == 0:
					outfile1.write('\t'.join(list(map(str,[lib, n, dicoInfoLibs[n][lib]['total'], dicoInfoLibs[n][lib]['paired'], dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['mapped'], dicoInfoLibs[n][lib]['mapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['unmapped'], dicoInfoLibs[n][lib]['unmapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['PairedMapped'], dicoInfoLibs[n][lib]['PairedMapped']/float(dicoInfoLibs[n][lib]['paired'])*100, dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['PoperlyPaired'], dicoInfoLibs[n][lib]['PoperlyPaired']/0, dicoInfoLibs[n][lib]['MultipleHits'], dicoInfoLibs[n][lib]['MultipleHits']/float(dicoInfoLibs[n][lib]['mapped'])*100]))))
				else:
					outfile1.write('\t'.join(list(map(str,[lib, n, dicoInfoLibs[n][lib]['total'], dicoInfoLibs[n][lib]['paired'], dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['mapped'], dicoInfoLibs[n][lib]['mapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['unmapped'], dicoInfoLibs[n][lib]['unmapped']/float(dicoInfoLibs[n][lib]['total'])*100, dicoInfoLibs[n][lib]['PairedMapped'], dicoInfoLibs[n][lib]['PairedMapped']/float(dicoInfoLibs[n][lib]['paired'])*100, dicoInfoLibs[n][lib]['single'], dicoInfoLibs[n][lib]['PoperlyPaired'], dicoInfoLibs[n][lib]['PoperlyPaired']/float(dicoInfoLibs[n][lib]['PairedMapped'])*100, dicoInfoLibs[n][lib]['MultipleHits'], dicoInfoLibs[n][lib]['MultipleHits']/float(dicoInfoLibs[n][lib]['mapped'])*100]))))
				outfile1.write('\n')
			else:
				sys.stdout.write('This is embarrassing: the file '+'/'.join([n,lib+'.sam.stat'])+' has not been found\n')
	
		if dicoInfoGlob[n]['total'] == 0:
			outfile2.write('\t'.join(list(map(str,[n, dicoInfoGlob[n]['total'], dicoInfoGlob[n]['paired'], dicoInfoGlob[n]['single'], dicoInfoGlob[n]['mapped'], 0, dicoInfoGlob[n]['unmapped'], 0, dicoInfoGlob[n]['PairedMapped'], 0, dicoInfoGlob[n]['single'], dicoInfoGlob[n]['PoperlyPaired'], 0, dicoInfoGlob[n]['MultipleHits'], 0]))))
		elif dicoInfoGlob[n]['mapped'] == 0:
			outfile2.write('\t'.join(list(map(str,[n, dicoInfoGlob[n]['total'], dicoInfoGlob[n]['paired'], dicoInfoGlob[n]['single'], dicoInfoGlob[n]['mapped'], dicoInfoGlob[n]['mapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['unmapped'], dicoInfoGlob[n]['unmapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['PairedMapped'], dicoInfoGlob[n]['PairedMapped']/float(dicoInfoGlob[n]['paired'])*100, dicoInfoGlob[n]['single'], dicoInfoGlob[n]['PoperlyPaired'], 0, dicoInfoGlob[n]['MultipleHits'], 0]))))
		elif dicoInfoGlob[n]['paired'] == 0:
			outfile2.write('\t'.join(list(map(str,[n, dicoInfoGlob[n]['total'], dicoInfoGlob[n]['paired'], dicoInfoGlob[n]['single'], dicoInfoGlob[n]['mapped'], dicoInfoGlob[n]['mapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['unmapped'], dicoInfoGlob[n]['unmapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['PairedMapped'], 0, dicoInfoGlob[n]['single'], dicoInfoGlob[n]['PoperlyPaired'], 0, dicoInfoGlob[n]['MultipleHits'], dicoInfoGlob[n]['MultipleHits']/float(dicoInfoGlob[n]['mapped'])*100]))))
		elif dicoInfoGlob[n]['PairedMapped'] == 0:
			outfile2.write('\t'.join(list(map(str,[n, dicoInfoGlob[n]['total'], dicoInfoGlob[n]['paired'], dicoInfoGlob[n]['single'], dicoInfoGlob[n]['mapped'], dicoInfoGlob[n]['mapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['unmapped'], dicoInfoGlob[n]['unmapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['PairedMapped'], dicoInfoGlob[n]['PairedMapped']/float(dicoInfoGlob[n]['paired'])*100, dicoInfoGlob[n]['single'], dicoInfoGlob[n]['PoperlyPaired'], dicoInfoGlob[n]['PoperlyPaired']/0, dicoInfoGlob[n]['MultipleHits'], dicoInfoGlob[n]['MultipleHits']/float(dicoInfoGlob[n]['mapped'])*100]))))
		else:
			outfile2.write('\t'.join(list(map(str,[n, dicoInfoGlob[n]['total'], dicoInfoGlob[n]['paired'], dicoInfoGlob[n]['single'], dicoInfoGlob[n]['mapped'], dicoInfoGlob[n]['mapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['unmapped'], dicoInfoGlob[n]['unmapped']/float(dicoInfoGlob[n]['total'])*100, dicoInfoGlob[n]['PairedMapped'], dicoInfoGlob[n]['PairedMapped']/float(dicoInfoGlob[n]['paired'])*100, dicoInfoGlob[n]['single'], dicoInfoGlob[n]['PoperlyPaired'], dicoInfoGlob[n]['PoperlyPaired']/float(dicoInfoGlob[n]['PairedMapped'])*100, dicoInfoGlob[n]['MultipleHits'], dicoInfoGlob[n]['MultipleHits']/float(dicoInfoGlob[n]['mapped'])*100]))))
		outfile2.write('\n')
	
	outfile1.close()
	outfile2.close()