
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
import time
import gzip
import shutil
import tempfile
import itertools
import subprocess
from inspect import currentframe, getframeinfo
from Bio.Seq import Seq
from Bio import SeqIO

def stop_err( msg ):
	sys.stderr.write( "%s\n" % msg )
	sys.exit()

def run_job (frameinfo, cmd_line, ERROR):
	try:
		tmp = tempfile.NamedTemporaryFile().name
		error = open(tmp, 'w')
		proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
		returncode = proc.wait()
		error.close()
		error = open( tmp, 'r' )
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
			raise Exception
		# elif stderr:
			# raise Exception
	except Exception:
		stop_err( 'Line : '+str(frameinfo.lineno)+' - '+ERROR + str( stderr ) )

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

def Record_CDS_and_UTR(GFF3, dico_gene_gff3):
	
	"""
		Record in a dictionnary CDS and UTR positions.
		
		:param GFF3: A gff3 file.
		:type GFF3: gff3
		:return: Dictionary containing for each gene CDS, UTR and INTRON position.
		:rtype: void
	"""
	
	# initializing dictionaries
	dico_parent = {}
	# recording for each gene CDS and exon position
	file = open(GFF3)
	for line in file:
		data = line.replace('\n','').split('\t')
		if data:
			if data[0][0] != "#":
				# recording gene in dictionary
				if data[2] == 'gene':
					dico = gffline2dic(data[8])
					dico_gene_gff3[dico['ID']] = {}
					dico_gene_gff3[dico['ID']]['exon'] = set()
					dico_gene_gff3[dico['ID']]['chr'] = data[0]
				elif data[2] == 'mRNA':
					dico = gffline2dic(data[8])
					dico_parent[dico['ID']] = dico['Parent']
				elif data[2] == 'exon':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					debut = int(data[3])
					while debut <= int(data[4]):
						dico_gene_gff3[nom_gene]['exon'].add(debut)
						debut += 1
				elif data[2] == 'CDS':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					debut = int(data[3])
					while debut <= int(data[4]):
						dico_gene_gff3[nom_gene]['exon'].add(debut)
						debut += 1
	file.close()
	
	for n in dico_gene_gff3:
		if len(dico_gene_gff3[n]['exon']) == 0:
			sys.exit('This is embarrassing... The program exited without finishing because a gene without CDS and exon was found: '+n)

def gffline2dic(LINE):
	
	"""
		Return a dictionnary of the information contained in 9th column of a gff3 file
		
		:param LINE: The 9th column of a gff3 file.
		:type LINE: str
		:return: Dictionary containing each element contained in the line.
		:rtype: dictionary
	"""
	# initializing dictionary
	dico_line = {}
	# recording information in dictionary
	to_split = LINE.split(';')
	for n in to_split:
		data = n.split('=')
		dico_line[data[0]] = data[1]
	return dico_line

def calcul_cov(LOCA_PROGRAMS, SAM, TYPE, OUT):
	"""
		Calculate the coverage of a sam or bam file, site by site
		
		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam or bam file.
		:type SAM: str
		:param TYPE: The format of the input file
		:type TYPE: str ("sam" | "bam")
		:param OUT: The name of the output file
		:type OUT: str
		:return: void
	"""
	
	if TYPE == 'sam':
		#Convert the sam file to the bam format
		sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, SAM+'_BAM.bam')
		run_job(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam (calcul_cov):\n')
		SAM = SAM+'_BAM.bam'
	elif TYPE != 'bam':
		mot = TYPE+' argument passed in --type is not recognized'
		sys.exit(mot)
	
	cal_cov = ('%s depth %s > %s') % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	run_job(getframeinfo(currentframe()), cal_cov, 'Error in calculating coverage:\n')
	
	#remove the intermediate bam file.
	if TYPE == 'sam':
		os.remove(SAM)

def calculate_annotation_coverage(OUT, dico_GFF3, BAM, REF, LOCA_PROGRAMS):
	# calculating coverage
	calcul_cov(LOCA_PROGRAMS, BAM, 'bam', BAM+'.cov')
	
	# recording in dictionnary covered sites
	dico_cov = {}
	file = open(BAM+'.cov')
	for line in file:
		data = line.split()
		if data:
			if not(data[0] in dico_cov):
				dico_cov[data[0]] = set()
			dico_cov[data[0]].add(int(data[1])) 
	file.close()
	
	# calculating exon coverage proportion gene by gene
	outfile = open(OUT,'w')
	for n in dico_GFF3:
		total_len = len(dico_GFF3[n]['exon'])
		cov_len = 0
		if dico_GFF3[n]['chr'] in dico_cov:
			for k in dico_GFF3[n]['exon']:
				if k in dico_cov[dico_GFF3[n]['chr']]:
					cov_len += 1
		outfile.write('\t'.join([n, str(float(cov_len)/total_len*100)])+'\n')
	outfile.close()


	
##############################################
#          Process RNASeq Only
##############################################


def run_step_A_RNAseq(STAR, PROC, OUT_REF1, REF, SJDBO, JID, QUEUE):
	index = '%s --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s' %	(STAR, PROC, OUT_REF1, REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), index, 'Error in step A:\n')
	else:
		run_qsub(QUEUE, [index], PROC, JID+'-INDEX1', "12G", JID)
	sys.stdout.write("Step a: reference indexation done\n")
	sys.stdout.flush()

def run_step_B_RNAseq(TMP, DICO_LIB, STAR, PROC, PREFIX, STAR_OPT, OUT_REF1, QUEUE):
	# Concatenation of reads for junction identification
	os.mkdir(TMP)
	list_pair1 = []
	list_pair2 = []
	dico_single = set()
	for n in DICO_LIB:
		for j in DICO_LIB[n]:
			if len(DICO_LIB[n][j]) == 1:
				dico_single.add(DICO_LIB[n][j][0])
			elif len(DICO_LIB[n][j]) == 2:
				list_pair1.append(DICO_LIB[n][j][0])
				list_pair2.append(DICO_LIB[n][j][1])
			else:
				sys.exit('Problem in the configuration file in libraries section')
	
	# Mapping reads first time for junction identification
	if len(list_pair1) > 0:
	
		fichier_final = open(TMP+'/mate1.fq', "w")
		for i in list_pair1:
			shutil.copyfileobj(open(i, 'r'), fichier_final)
		fichier_final.close()
		
		fichier_final = open(TMP+'/mate2.fq', "w")
		for i in list_pair2:
			shutil.copyfileobj(open(i, 'r'), fichier_final)
		fichier_final.close()
		
		mappP = '%s --runThreadN %s --genomeDir %s --readFilesIn %s %s --outFileNamePrefix %s %s' % (STAR, PROC, OUT_REF1, TMP+'/mate1.fq', TMP+'/mate2.fq', PREFIX+'_JUNC_ESTIMATION_pair', STAR_OPT)
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), mappP, 'Error in step B:\n')
		else:
			run_qsub(QUEUE, [mappP], PROC, PREFIX+'-MapP1', "12G", PREFIX)
		# checking step
		if not(os.path.isfile(PREFIX+'_JUNC_ESTIMATION_pairAligned.out.sam')):
			sys.exit('An error was encountered in step b. The program exited without finishing during mapping step')
		
	if len(dico_single) > 0:
		fichier_final = open(TMP+'/single.fq', "w")
		for i in dico_single:
			shutil.copyfileobj(open(i, 'r'), fichier_final)
		fichier_final.close()
		mappS = '%s --runThreadN %s --genomeDir %s --readFilesIn %s --outFileNamePrefix %s %s' % (STAR, PROC, OUT_REF1, TMP+'/single.fq', PREFIX+'_JUNC_ESTIMATION_single', STAR_OPT)
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), mappS, 'Error in step B:\n')
		else:
			run_qsub(QUEUE, [mappS], PROC, PREFIX+'-MapS1', "12G", PREFIX)
		# checking step
		if not(os.path.isfile(PREFIX+'_JUNC_ESTIMATION_singleAligned.out.sam')):
			sys.exit('An error was encountered in step b. The program exited without finishing')

	# removing folders
	shutil.rmtree(TMP)
	shutil.rmtree(OUT_REF1)
	
	# Merging JUNCTION files if needed
	if len(list_pair1) > 0 and len(dico_single) > 0:
		dico_sites = {}
		dico_list_sites = {}
		file = open(PREFIX+'_JUNC_ESTIMATION_pairSJ.out.tab')
		for line in file:
			data = line.split()
			if data:
				line_id = ' '.join(data[0:3])
				intermediate_list = map(int, data[1:3])
				if not(data[0]) in dico_list_sites:
					dico_list_sites[data[0]] = []
				if intermediate_list in dico_list_sites[data[0]]:
					if int(data[3]) == dico_sites[line_id][0] and int(data[4]) == dico_sites[line_id][1] and int(data[5]) == dico_sites[line_id][2]:
						dico_sites[line_id][3] += int(data[6])
						dico_sites[line_id][4] += int(data[7])
						dico_sites[line_id][5] = max(int(data[7]), dico_sites[line_id][5])
					else:
						print (data[0:3], data[3], dico_sites[line_id][0], data[4], dico_sites[line_id][1], data[5], dico_sites[line_id][2])
						sys.exit('An error was encountered in step b. The program exited without finishing during mapping step')
				else:
					dico_list_sites[data[0]].append(intermediate_list)
					dico_sites[line_id] = map(int, data[3:])
		file.close()
		
		file = open(PREFIX+'_JUNC_ESTIMATION_singleSJ.out.tab')
		for line in file:
			data = line.split()
			if data:
				line_id = ' '.join(data[0:3])
				intermediate_list = map(int, data[1:3])
				if not(data[0]) in dico_list_sites:
					dico_list_sites[data[0]] = []
				if intermediate_list in dico_list_sites[data[0]]:
					if int(data[3]) == dico_sites[line_id][0] and int(data[4]) == dico_sites[line_id][1] and int(data[5]) == dico_sites[line_id][2]:
						dico_sites[line_id][3] += int(data[6])
						dico_sites[line_id][4] += int(data[7])
						dico_sites[line_id][5] = max(int(data[7]), dico_sites[line_id][5])
					else:
						print (data[0:3], data[3], dico_sites[line_id][0], data[4], dico_sites[line_id][1], data[5], dico_sites[line_id][2])
						sys.exit('An error was encountered in step b. The program exited without finishing during mapping step')
				else:
					dico_list_sites[data[0]].append(intermediate_list)
					dico_sites[line_id] = map(int, data[3:])
		file.close()
		
		outfile = open(PREFIX+'_JUNC_ESTIMATION_SJ.out.tab','w')
		for n in dico_list_sites:
			dico_list_sites[n].sort(cmp=lambda x,y: cmp(x[0],y[0]))
			for j in dico_list_sites[n]:
				line_id = ' '.join([n]+map(str, j))
				outfile.write('\t'.join(map(str,[n]+ j+dico_sites[line_id])))
				outfile.write('\n')
		outfile.close()
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleSJ.out.tab')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleAligned.out.sam')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.final.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.progress.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairSJ.out.tab')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairAligned.out.sam')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.final.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.progress.out')
	elif len(list_pair1) > 0:
		os.renames(PREFIX+'_JUNC_ESTIMATION_pairSJ.out.tab', PREFIX+'_JUNC_ESTIMATION_SJ.out.tab')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairAligned.out.sam')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.final.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_pairLog.progress.out')
	elif len(dico_single) > 0:
		os.renames(PREFIX+'_JUNC_ESTIMATION_singleSJ.out.tab', PREFIX+'_JUNC_ESTIMATION_SJ.out.tab')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleAligned.out.sam')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.final.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.out')
		os.remove(PREFIX+'_JUNC_ESTIMATION_singleLog.progress.out')
	else:
		sys.exit('An error was encountered in step b. The program exited without finishing. The problem is on fastq files...')
	sys.stdout.write("Step b: Spliced sites estimation done\n")
	sys.stdout.flush()

def run_step_C_RNAseq(STAR, PROC, OUT_REF2, REF, SJDBO, JID, SJDBFCSE, QUEUE):
	index = '%s --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbOverhang %s --sjdbFileChrStartEnd %s' % (STAR, PROC, OUT_REF2, REF, SJDBO, SJDBFCSE)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), index, 'Error in step C:\n')
	else:
		run_qsub(QUEUE, [index], PROC, JID+'-INDEX2', "12G", JID)
	sys.stdout.write("Step c: reference indexation done\n")
	sys.stdout.flush()

def run_step_D_RNAseq (LIBS, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, QUEUE):
	if not (os.path.isdir(ACC)):
		os.mkdir(ACC)
	for lib in LIBS:
		rg_tag = 'ID:'+ACC+'"\t"SM:'+ACC+'"\t"LB:'+ACC+'"\t"PU:whatever"\t"PL:ILLUMINA'
		if len(LIBS[lib]) == 1:
			mapp2 = '%s --runThreadN %s --genomeDir %s --readFilesIn %s --outFileNamePrefix %s --outSAMattrRGline %s %s' % (STAR, PROC, OUT_REF2, LIBS[lib][0], ACC+'/'+lib, rg_tag, STAR_OPT)
		elif len(LIBS[lib]) == 2:
			mapp2 = '%s --runThreadN %s --genomeDir %s --readFilesIn %s %s --outFileNamePrefix %s --outSAMattrRGline %s %s' % (STAR, PROC, OUT_REF2, LIBS[lib][0], LIBS[lib][1], ACC+'/'+lib, rg_tag, STAR_OPT)
		else:
			return 'Problem in the configuration file in libraries section for accession '+ACC+'\n'
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), mapp2, 'Error in step D:\n')
		else:
			run_qsub(QUEUE, [mapp2], PROC, ACC+'-Map2', "12G", PREFIX)
		# checking step
		if not(os.path.isfile(ACC+'/'+lib+'Aligned.out.sam')):
			return 'An error was encountered in step d. The program exited without finishing for accession '+ACC+'\n'
		os.remove(ACC+'/'+lib+'Log.progress.out')
		os.remove(ACC+'/'+lib+'Log.out')
	sys.stdout.write("Step d: Mapping done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_stat_step_D_RNAseq(PREFIX, CONFIG):
	outfile = open(PREFIX+'/'+PREFIX+'_mapping.tab','w')
	outfile.write('Accession\tPairedReads\tSingleReads\tUniquelyMappedPairs\tMultipleMappedPairs\tUniquelyMappedSingle\tMultipleMappedSingle\n')
	dico_lib_stat = {}
	for n in CONFIG.options('Libraries'):
		data =  CONFIG.get('Libraries', n).split()
		if not(data[0] in dico_lib_stat):
			dico_lib_stat[data[0]] = []
		if len(data[1:-1]) == 2:
			dico_lib_stat[data[0]].append([n,'pair'])
		elif len(data[1:-1]) == 1:
			dico_lib_stat[data[0]].append([n,'single'])
		else:
			sys.exit('An error was encountered in step d because an problem was encountered in the configuration file in libraries section. The program exited without finishing')
	for acc in dico_lib_stat:
		nb_single = 0
		nb_mapped_single = 0
		nb_multi_single = 0
		nb_pair = 0
		nb_mapped_pair = 0
		nb_multi_pair = 0
		for k in dico_lib_stat[acc]:
			file = open(acc+'/'+k[0]+'Log.final.out')
			for line in file:
				data = line.split('|')
				if data:
					if data[0] == '                          Number of input reads ':
						number = int(data[1].split()[0])
						if k[1] == 'pair':
							nb_pair += number
						else:
							nb_single += number
					elif data[0] == '                   Uniquely mapped reads number ':
						number = int(data[1].split()[0])
						if k[1] == 'pair':
							nb_mapped_pair += number
						else:
							nb_mapped_single += number
					elif data[0] == '        Number of reads mapped to multiple loci ':
						number = int(data[1].split()[0])
						if k[1] == 'pair':
							nb_multi_pair += number
						else:
							nb_multi_single += number
					elif data[0] == '        Number of reads mapped to too many loci ':
						number = int(data[1].split()[0])
						if k[1] == 'pair':
							nb_multi_pair += number
						else:
							nb_multi_single += number
		outfile.write('\t'.join([acc,str(nb_pair),str(nb_single),str(nb_mapped_pair),str(nb_multi_pair),str(nb_mapped_single),str(nb_multi_single)])+'\n')
	outfile.close()

def run_step_E_RNAseq(LIB, ACC, PREFIX, JAVA, PICARD, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	to_merge = ''
	for n in LIB:
		to_merge = to_merge + 'INPUT='+ACC+'/'+n+'Aligned.out.sam '
	merge = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s MergeSamFiles %s OUTPUT=%s MERGE_SEQUENCE_DICTIONARIES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=%s' % (JAVA, PICARD, to_merge, ACC+'/'+ACC+'_merged.bam', ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), merge, 'Error in step E:\n')
	else:
		run_qsub(QUEUE, [merge], 1, ACC+'-Merge', "12G", PREFIX)
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_merged.bam')):
		return 'An error was encountered in step e. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	for n in LIB:
		os.remove(ACC+'/'+n+'Aligned.out.sam')
	sys.stdout.write("Step e: Merging done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_F_RNAseq(ACC, JAVA, PICARD, PREFIX, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	rmdup = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true QUIET=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VERBOSITY=WARNING TMP_DIR=%s' % (JAVA, PICARD, ACC+'/'+ACC+'_merged.bam', ACC+'/'+ACC+'_rmdup.bam', ACC+'/'+ACC+'_duplicate', ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), rmdup, 'Error in step F:\n')
	else:
		run_qsub(QUEUE, [rmdup], 1, ACC+'-RmDup', "12G", PREFIX)
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_rmdup.bam')):
		return 'An error was encountered in step f. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	os.remove(ACC+'/'+ACC+'_merged.bam')
	os.remove(ACC+'/'+ACC+'_merged.bai')
	sys.stdout.write("Step f: duplicate removal done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_stat_step_F_RNAseq(PREFIX, DICO_LIB):
	outfile = open(PREFIX+'/'+PREFIX+'_rmdup_stat.tab','w')
	outfile.write('LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n')
	for acc in DICO_LIB:
		file = open(acc+'/'+acc+'_duplicate')
		i = 0
		while i < 8:
			i += 1
			line = file.readline()
		outfile.write(line)
	outfile.close()

def run_step_G_RNAseq(JAVA, PICARD, ACC, REF, PREFIX, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	reorder = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s CREATE_INDEX=true TMP_DIR=%s' % (JAVA, PICARD, ACC+'/'+ACC+'_rmdup.bam',  ACC+'/'+ACC+'_reorder.bam', REF, ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), reorder, 'Error in step G:\n')
	else:
		run_qsub(QUEUE, [reorder], 1, ACC+'-Reord', "12G", PREFIX)
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_reorder.bam')):
		return 'An error was encountered in step g. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	os.remove(ACC+'/'+ACC+'_rmdup.bam')
	sys.stdout.write("Step g: Bam reodering done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_H_RNAseq(JAVA, GATK, ACC, REF, PREFIX, QUEUE):
	trim = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T SplitNCigarReads -I %s -o %s -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -R %s -U ALLOW_N_CIGAR_READS' % (JAVA, GATK, ACC+'/'+ACC+'_reorder.bam', ACC+'/'+ACC+'_trim.bam', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), trim, 'Error in step H:\n')
	else:
		run_qsub(QUEUE, [trim], 1, ACC+'-Trim', "12G", PREFIX)
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_trim.bam')):
		return 'An error was encountered in step h. The program exited without finishing for accession '+ACC+'\n'
	os.remove(ACC+'/'+ACC+'_reorder.bam')
	os.remove(ACC+'/'+ACC+'_reorder.bai')
	sys.stdout.write("Step h: Read splitting done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_I_RNAseq(ACC, JAVA, GATK, REF, PLOIDY, PREFIX, UseUnifiedGenotyperForBaseRecal, QUEUE):
	
	realTC = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T RealignerTargetCreator -o %s -I %s -R %s' % (JAVA, GATK, ACC+'/'+ACC+'_RTC.intervals', ACC+'/'+ACC+'_trim.bam', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), realTC, 'Error in step I:\n')
	else:
		run_qsub(QUEUE, [realTC], 1, ACC+'-realTC', "12G", PREFIX)
	if not(os.path.isfile(ACC+'/'+ACC+'_RTC.intervals')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'

	Ireal = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T IndelRealigner -o %s -I %s -targetIntervals %s -R %s' % (JAVA, GATK, ACC+'/'+ACC+'_realigned.bam', ACC+'/'+ACC+'_trim.bam', ACC+'/'+ACC+'_RTC.intervals', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), Ireal, 'Error in step I:\n')
	else:
		run_qsub(QUEUE, [Ireal], 1, ACC+'-Ireal', "12G", PREFIX)
	if not(os.path.isfile(ACC+'/'+ACC+'_realigned.bam')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'
	os.remove(ACC+'/'+ACC+'_RTC.intervals')
	
	if not(os.path.isfile(ACC+'/'+ACC+'_realigned.bam')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'
	sys.stdout.write("Step i: Base recalibration done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_M_RNAseq(GFF3, REF, LOCA_PROGRAMS, PREFIX, DICO_LIB):
	# recording gff3 informations
	sys.stdout.write('recording gff3 informations\n')
	dico_gene_gff3 = {}
	Record_CDS_and_UTR(GFF3, dico_gene_gff3)
	
	# calculating exon coverage by accession
	sys.stdout.write('calculating exon coverage by accession\n')
	list_lib = []
	list_bam_cov = []
	for ACC_ID in DICO_LIB:
		sys.stdout.write(ACC_ID+'\n')
		sys.stdout.flush()
		if os.path.isfile(ACC_ID+'/'+ACC_ID+'_real_recal.bam'):
			bam = ACC_ID+'/'+ACC_ID+'_real_recal.bam'
			list_bam_cov.append(bam+'.cov')
		elif os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam'):
			bam = ACC_ID+'/'+ACC_ID+'_realigned.bam'
			list_bam_cov.append(bam+'.cov')
			sys.stdout.write('Warning no recalibrated bam file found for '+ACC_ID+'... Continue using realigned bam instead.\n')
			sys.stdout.flush()
		else:
			return 'ERROR : Neither '+ACC_ID+'/'+ACC_ID+'_real_recal.bam or '+ACC_ID+'/'+ACC_ID+'_realigned.bam were found.\n'
		list_lib.append(ACC_ID)
		calculate_annotation_coverage(ACC_ID+'/'+ACC_ID+'_exon_coverage.cov', dico_gene_gff3, bam, REF, LOCA_PROGRAMS)

	
	# Loading final dictionary
	final_dico = {}
	for acc in list_lib:
		file = open(acc+'/'+acc+'_exon_coverage.cov')
		for line in file:
			data = line.split()
			if data:
				if not(data[0] in final_dico):
					final_dico[data[0]] = {}
				final_dico[data[0]][acc] = data[1]
		file.close()
	
	# merging all informations in a single file
	outfile = open(PREFIX+'/'+PREFIX+'_exon_coverage.cov','w')
	outfile.write('gene_id\t'+'\t'.join(list_lib)+'\n')
	for acc in final_dico:
		outfile.write(acc)
		for k in list_lib:
			outfile.write('\t'+final_dico[acc][k])
		outfile.write('\n')
	outfile.close()
	
	for acc in list_bam_cov:
		os.remove(acc)

##############################################
#          Process ReSeq Only
##############################################

def run_step_A(ACC_ID, LIB_DIC, BWA, REF, TMP, JAVA, PICARD, SAMTOOLS, PREFIX, QUEUE, PARSEUNMAPPED):

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
			sys.stdout.write(comd+'\n')
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
			sys.stdout.write(comd+'\n')
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (Map Stat):\n')
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
			sys.stdout.write(comd+'\n')
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (Plot stat):\n')
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
			sys.stdout.write(comd+'\n')
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
			sys.stdout.write(comd+'\n')
			run_job(getframeinfo(currentframe()), comd, 'Error in step A (picard tools):\n')
	else:
		run_qsub(QUEUE, to_merge, 1, 'merge-'+ACC_ID, "12G", PREFIX)
	sys.stdout.write('done\n')
	
	# checking step
	if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_merged.bam')):
		return 'An error was encountered in step a during merging step. The program exited without finishing on '+ACC_ID+'\n'
	
	#############
	#3 Parsing unmapped read if requested
	if PARSEUNMAPPED == 'y':
		sys.stdout.write('Initiating unmapped removal.\n')
		to_filter = []
		for n in LIB_DIC:
			unMapkeep = '%s view -bS -uf 4 %s > %s' % (SAMTOOLS, ACC_ID+'/'+TMP+n+'.sam', ACC_ID+'/'+n+'Unmapped.bam')
			to_filter.append(unMapkeep)
		
		if QUEUE == None:
			for comd in to_filter:
				sys.stdout.write(comd+'\n')
				run_job(getframeinfo(currentframe()), comd, 'Error in step A (samtools):\n')
		else:
			run_qsub(QUEUE, to_filter, 2, 'filter-'+ACC_ID, "12G", PREFIX)
		sys.stdout.write('done\n')
		
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

##############################################
#          Common to Process
##############################################

def run_step_E(ACC_ID, PYTHON, REF, DICO_CHR, PREFIX, QUEUE, PATHNAME):

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

def genotype_accession(COVERAGE, ALLELE, ERROR, PLOIDY, LRT, ALLCOMB):
	
	# WARNING The statistics will not be accurate on low coverage
	
	
	# To manage huge factorial
	nb_tirage = sum(COVERAGE)
	if nb_tirage >= 1000:
		coverage = list(map(int, [x/float(nb_tirage)*100 for x in COVERAGE]))
	else:
		coverage = list(COVERAGE)
	nb_tirage = sum(coverage)
	
	# Obtaining combinations
	if ALLCOMB:
		FinalComb = list(itertools.combinations_with_replacement(range(len(ALLELE)), int(PLOIDY)))
	else:
		ListCombPair = list(itertools.combinations_with_replacement(range(len(ALLELE)), 2))
		FinalComb = set()
		for n in ListCombPair:
			N = set(n)
			if len(N) == 2:
				IntermList = list(itertools.combinations_with_replacement(n, int(PLOIDY)))
				for k in IntermList:
					FinalComb.add(k)
	
	dico_proba = {}
	# calculating probabilities for homozygous state
	p = 1-ERROR
	for n in FinalComb:
		N = set(n)
		Key = '/'.join(list(map(str, n)))
		if len(N) == 1:
			Value = N.pop()
			Proba = binom(coverage[Value],nb_tirage,p)/binom(int(round(nb_tirage*p)),nb_tirage,p)
		else:
			Proba = 1
			for k in N:
				Ratio = n.count(k)/float(PLOIDY)
				Proba = Proba*(binom(coverage[k], nb_tirage, Ratio)/binom(int(round(nb_tirage*Ratio)),nb_tirage,Ratio))
		dico_proba[Key] = Proba
	
	if LRT:
		# getting best and second best probability
		best_genotype = '/'.join(['.']*int(PLOIDY))
		best_value = 0
		second_best = 0
		for genotype in dico_proba:
			if best_value < dico_proba[genotype]:
				if best_value != 0:
					second_best = best_value
				best_value = dico_proba[genotype]
				best_genotype = genotype
			elif best_value == dico_proba[genotype]:
				best_genotype = '/'.join(['.']*int(PLOIDY))
			elif second_best < dico_proba[genotype]:
				second_best = dico_proba[genotype]
		if best_value == 0:
			return best_genotype, (1)
		elif second_best == 0:
			return best_genotype, (9999)
		else:
			return best_genotype, best_value/second_best
	
	else:
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
	outfile = gzip.open(PREFIX+'_'+CHR+'_'+str(START)+'_'+str(END)+'_allele_count.vcf.gz', 'wt')
	outfile.write("##fileformat=VCFv4.2\n")
	outfile.write("##reference=file:///"+REF+"\n")
	outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	outfile.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
	for n in DICO_CHR:
		outfile.write("##contig=<ID="+n+",length="+str(DICO_CHR[n])+">\n")
	liste2print = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
	for acc in LIST_ACC:
		ACC = acc.split('/')[-1]
		liste2print.append(ACC)
	outfile.write('\t'.join(liste2print))
	outfile.write('\n')
	
	#2- Initiating reading files and variables
	dico_accession_infile = {}
		
	dico_accession_infile[CHR] = {}
	# sys.stdout.write(CHR+'\n')
	for acc in LIST_ACC:
		ACC = acc.split('/')[-1]
		dico_accession_infile[CHR][acc] = gzip.open(acc+'/'+ACC+'_allele_count_'+CHR+'.gz', 'rt')

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
		for allele in sorted(dico_allele_cov.keys()):
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
					genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.005, DICO_PLOIDY[acc], False, False)
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
			for allele in sorted(dico_allele_cov.keys()):
				if allele != reference:
					if dico_allele_cov[allele]:
						allele_to_keep.append(allele)
			
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
						genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.005, DICO_PLOIDY[acc], False, False)
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

def create_pseudo_VCF_Large(LIST_ACC, REF, PREFIX, DICO_PLOIDY, DICO_CHR, CHR, START, END):
	
	#1- preparing output
	outfile = gzip.open(PREFIX+'_'+CHR+'_'+str(START)+'_'+str(END)+'_allele_count.vcf.gz', 'wt')
	outfile.write("##fileformat=VCFv4.2\n")
	outfile.write("##reference=file:///"+REF+"\n")
	outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	outfile.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
	for n in DICO_CHR:
		outfile.write("##contig=<ID="+n+",length="+str(DICO_CHR[n])+">\n")
	liste2print = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
	for acc in LIST_ACC:
		ACC = acc.split('/')[-1]
		liste2print.append(ACC)
	outfile.write('\t'.join(liste2print))
	outfile.write('\n')
	
	# Opening big file
	file = gzip.open(PREFIX+'_'+CHR+'_Back.gz', 'rt')
	
	# Recording accession index in the line
	AccPosIndex = {}
	header = file.readline().split()
	LineLen = len(header) # just for verification
	for n in LIST_ACC:
		if not(n in header):
			sys.stdout.write('The accession '+n+' was not found in the file '+PREFIX+'_'+CHR+'Back.gz. Please check your configuration file.\n')
			return 1
		else:
			AccPosIndex[n] = header.index(n)
	
	# Initiating studied alleles (only SNP are managed)
	liste_allele = ['A', 'C', 'G', 'T', '*']
	
	# Performing the calling
	for line in file:
		data = line.split()
		if data:
			POS = int(data[1])
			if POS > START and POS <= END:
				#1- recording total allele coverage information and reference allele
				reference = data[2]
				format = data[3].split(':')
				
				dico_allele_cov = {}
				for allele in liste_allele:
					dico_allele_cov[allele] = 0
				
				for acc in LIST_ACC:
					allele_count = data[AccPosIndex[acc]].split(':')
					for allele in dico_allele_cov:
						dico_allele_cov[allele] += get_allele_coverage (allele, format, allele_count)
			
				#2- recording expected allele
				allele_to_keep = [reference]
				for allele in sorted(dico_allele_cov.keys()):
					if allele != reference:
						if dico_allele_cov[allele]:
							allele_to_keep.append(allele)
				
				#3- Calculating genotype and printing results to output if necessary
				if len(allele_to_keep) > 1:
					liste2print = [CHR, str(POS), '.', allele_to_keep[0], ','.join(allele_to_keep[1:]), '.', '.', '.', 'GT:AD:DP']
					for acc in LIST_ACC:
						liste_cov = []
						liste_cov_int = []
						allele_count = data[AccPosIndex[acc]].split(':')
						depth = allele_count[0]
						for allele in allele_to_keep:
							allele_coverage = get_allele_coverage (allele, format, allele_count)
							liste_cov.append(str(allele_coverage))
							liste_cov_int.append(allele_coverage)
						genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.005, DICO_PLOIDY[acc], False, False)
						liste2print.append(':'.join([genotype,','.join(liste_cov),depth]))
					outfile.write('\t'.join(liste2print))
					outfile.write('\n')
	
	file.close()
	outfile.close()
	
	return 0

def merge_vcf(PREFIX, DICO_CHR):
	
	# recording VCF headers
	dico_header = {}
	accessions = set()
	for chr in DICO_CHR:
		if os.path.exists(PREFIX+'_'+chr+'_all_allele_count.vcf'):
			file = open(PREFIX+'_'+chr+'_all_allele_count.vcf','r')
		elif os.path.exists(PREFIX+'_'+chr+'_all_allele_count.vcf.gz'):
			file = gzip.open(PREFIX+'_'+chr+'_all_allele_count.vcf.gz','rb')
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

def merge_sub_vcf(PREFIX, CHR, LIST, GZIP):
	
	# recording VCF headers
	dico_header = {}
	accessions = []
	for pos in LIST:
		POS = '-'.join(map(str,pos))
		file = gzip.open(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf.gz','rt')
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
	if GZIP == 'n':
		outfile = open(PREFIX+'_'+CHR+'_all_allele_count.vcf','w')
	elif GZIP == 'y':
		outfile = gzip.open(PREFIX+'_'+CHR+'_all_allele_count.vcf'+'.gz','wt')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	outfile.write(''.join(dico_header['header']))
	outfile.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+accessions)+'\n')
	for pos in LIST:
		POS = '-'.join(map(str,pos))
		file = gzip.open(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf.gz','rt')
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
		os.remove(PREFIX+'_'+CHR+'_'+str(pos[0])+'_'+str(pos[1])+'_allele_count.vcf.gz')
	return 0

def perform_count(BAMREADCOUNT, REF, BAM, OUT):
	
	bam_count = '%s -b 10 -q 1 -d 1000000 -w 0 -f %s %s | gzip > %s' % (BAMREADCOUNT, REF, BAM, OUT+'.gz')
	os.system(bam_count)

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
	file = gzip.open(TAB+'.gz', 'rt')
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
			if data[3] != "0":
				if chr != data[0]:
					if not(chr == ""):
						outfile.close()
					chr = data[0]
					outfile = gzip.open(OUT+'_'+chr+'.gz', 'wt')
				# Recording information on the current line
				for n in data[4:]:
					if n[0] != "=": # To remove the column beginning with "=" which I don't know what it means...
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
			outfile = gzip.open(OUT+'_'+n+'.gz', 'wt')
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

def createGVCF(CHR, ACCS, PREFIX):
	
	LargeFile = PREFIX+'_'+CHR+'_Back.gz'
	IntermFile = PREFIX+'_'+CHR+'_New.gz'
	
	for i in range(len(ACCS)):
		acc = ACCS[i]
		AccName = acc.split('/')[-1]
		AccFile = acc+'/'+AccName+'_allele_count_'+CHR+'.gz'
		if os.path.exists(LargeFile):
			file1 = gzip.open(LargeFile, 'rt')
			file2 = gzip.open(AccFile, 'rt')
			header = file1.readline().split()
			AccNumber = len(header)-4
			if AccName in header:
				sys.stdout.write('The accession '+AccName+' was already found in the file. It will thus not be included once again\n')
				file1.close()
				file2.close()
			else:
				print(AccName, CHR)
				outfile = gzip.open(IntermFile,'wt')
				header.append(AccName)
				outfile.write('\t'.join(header))
				outfile.write('\n')
				
				# Initiating data aggregation
				line1 = file1.readline()
				line2 = file2.readline()
				while line1 and line2:
					data1 = line1.split()
					data2 = line2.split()
					pos1 = int(data1[1])
					pos2 = int(data2[1])
					# print(pos1, pos2, AccName)
					if pos1 == pos2:
						GENO1 = data1[3].split(':')
						GENO2 = data2[4].split(':')
						Geno1 = set(GENO1)
						Geno2 = set(GENO2)
						# Cov1 = data1[4].split(':')
						Cov2 = data2[5].split(':')
						Ref1 = data1[2]
						Ref2 = data2[2]
						if Ref1 != Ref2:
							sys.exit('Oups, there is a bug 1\n')
						# Identification of what is missing from Geno1
						NotInBigFile = Geno2.difference(Geno1)
						if NotInBigFile:
							
							data1[3] = ':'.join(data1[3].split(':')+list(NotInBigFile))
							
							for Var in NotInBigFile:
								GENO1.append(Var)
							for i in range(len(data1[4:])):
									data1[i+4] = ':'.join(data1[i+4].split(':')+['0']*len(NotInBigFile))
						
						# Preparing final Genotype information for current accession
						FinalGENO = data1[3].split(':')
						
						FinalAccCov = [data2[3]]
						for Var in FinalGENO[1:]:
							if Var in Geno2:
								FinalAccCov.append(Cov2[GENO2.index(Var)])
							else:
								FinalAccCov.append('0')
						
						data1 = data1 + [':'.join(FinalAccCov)]
						
						outfile.write('\t'.join(data1))
						outfile.write('\n')
						
						line1 = file1.readline()
						line2 = file2.readline()
					elif pos1 > pos2: # (1) Adding a new line in the global file with missing data for all accessions already present in this vcf (2) reading new line of file2
					
						# print('yes1')
						GENO2 = '0:'+data2[4]
						Cov2 = data2[3]+':'+data2[5]
						MissCov = [':'.join(['0']*len(Cov2.split(':')))]*AccNumber
						DataToPrint = data2[0:3]+[GENO2]+MissCov+[Cov2]
						
						outfile.write('\t'.join(DataToPrint))
						outfile.write('\n')
						
						line2 = file2.readline()
					
					elif pos1 < pos2: # (1) Putting missing data to accession of file 2 (2) reading new line of file1
					
						# print('yes2')
						GENO1 = data1[3].split(':')
						MissCov = ':'.join(['0']*len(GENO1))
						data1 = data1+[MissCov]
						
						outfile.write('\t'.join(data1))
						outfile.write('\n')
						
						line1 = file1.readline()
				
				if line1: # We must put missing to the new accession for all remaining sites
					# print('yes3')
					data1 = line1.split()
					GENO1 = data1[3].split(':')
					MissCov = ':'.join(['0']*len(GENO1))
					data1 = data1+[MissCov]
					outfile.write('\t'.join(data1))
					outfile.write('\n')
					
					for line1 in file1:
						data1 = line1.split()
						GENO1 = data1[3].split(':')
						MissCov = ':'.join(['0']*len(GENO1))
						data1 = data1+[MissCov]
						outfile.write('\t'.join(data1))
						outfile.write('\n')
					
				elif line2: # We must put missing to all "old" accessions for all remaining sites
					
					# print('yes4')
					data2 = line2.split()
					GENO2 = '0:'+data2[4]
					Cov2 = data2[3]+':'+data2[5]
					MissCov = [':'.join(['0']*len(Cov2.split(':')))]*AccNumber
					DataToPrint = data2[0:3]+[GENO2]+MissCov+[Cov2]
					outfile.write('\t'.join(DataToPrint))
					outfile.write('\n')
					
					for line2 in file2:
						data2 = line2.split()
						GENO2 = '0:'+data2[4]
						Cov2 = data2[3]+':'+data2[5]
						MissCov = [':'.join(['0']*len(Cov2.split(':')))]*AccNumber
						DataToPrint = data2[0:3]+[GENO2]+MissCov+[Cov2]
						outfile.write('\t'.join(DataToPrint))
						outfile.write('\n')
				
				file1.close()
				file2.close()
				outfile.close()
				
				# Renaming New file
				os.rename(IntermFile, LargeFile)
			
		else:
			file1 = gzip.open(AccFile, 'rt')
			outfile = gzip.open(LargeFile,'wt')
			outfile.write('\t'.join(['#CHROM', 'POS', 'REF', 'ALLELE', AccName]))
			outfile.write('\n')
			for line in file1:
				data = line.split()
				if data:
					outfile.write('\t'.join(data[0:3]+['0:'+data[4], data[3]+':'+data[5]]))
					outfile.write('\n')
			outfile.close()
	return 0
			
	
	
	