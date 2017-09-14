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
import ConfigParser
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
sys.stdout.write("module loaded\n")
sys.stdout.flush()

def get_min_acc_pos(DICO):
	
	min_pos = False
	for n in DICO:
		if min_pos:
			min_pos = min(min_pos, DICO[n]) 
		else:
			min_pos = DICO[n]
	return min_pos

def get_allele_coverage (ALLELE, FORMAT, LISTE):
	
	if ALLELE == 'N':
		return 0
	else:
		return int(LISTE[FORMAT.index(ALLELE)])

def create_pseudo_VCF(LIST_ACC, REF, PREFIX, DICO_PLOIDY):
	
	#0- getting sequence length
	dico_chr = {}
	sequence_dict = SeqIO.index(REF, "fasta")
	for n in sequence_dict:
		dico_chr[n] = len(str(sequence_dict[n].seq))
	del sequence_dict
	
	#1- preparing output
	outfile = open(PREFIX+'_allele_count.vcf', 'w')
	outfile.write("##fileformat=VCFv4.2\n")
	outfile.write("##reference=file:///"+REF+"\n")
	outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	outfile.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
	for n in dico_chr:
		outfile.write("##contig=<ID="+n+",length="+str(dico_chr[n])+">\n")
	liste2print = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
	for acc in LIST_ACC:
		liste2print.append(acc)
	outfile.write('\t'.join(liste2print))
	outfile.write('\n')
	
	#2- Initiating reading files and variables
	dico_accession_infile = {}
	
	list_chr = list(dico_chr.keys())
	while list_chr:
		chr = list_chr[0]
		del list_chr[0]
		
		dico_accession_infile[chr] = {}
		sys.stdout.write(chr+'\n')
		for acc in LIST_ACC:
			dico_accession_infile[chr][acc] = gzip.open(acc+'/'+acc+'_allele_count_'+chr+'.gz', 'rb')
			# print chr, acc
	
		#3- Initiating variables
		dico_accession_positions = {}
		dico_accession_line = {}
		for acc in LIST_ACC:
			dico_accession_positions[acc] = 0
			dico_accession_line[acc] = 0
		
		liste_allele = ['A', 'C', 'G', 'T']
			
		#4a- Initiating the loop
		
		#4aa- getting first line
		old_position = 0
		for acc in LIST_ACC:
			# initiating the line
			dico_accession_line[acc] = dico_accession_infile[chr][acc].readline().split()
		
		#4ab- recording position informations
		for acc in LIST_ACC:
			if dico_accession_line[acc] == []:
				dico_accession_positions[acc] = 1000000000000000
			else:
				dico_accession_positions[acc] = int(dico_accession_line[acc][1])
		min_position = get_min_acc_pos(dico_accession_positions)
		# print min_position
		
		# All file are not empty
		if min_position != 1000000000000000:
		
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
			# print dico_allele_cov
			
			#4ad- recording expected allele
			allele_to_keep = [reference]
			for allele in dico_allele_cov:
				if allele != reference:
					if dico_allele_cov[allele]:
						allele_to_keep.append(allele)
			# print allele_to_keep
			
			#4ae- printing results to output if necessary
			if len(allele_to_keep) > 1:
				liste2print = [chr, str(min_position), '.', allele_to_keep[0], ','.join(allele_to_keep[1:]), '.', '.', '.', 'GT:AD:DP']
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
						genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.01, DICO_PLOIDY[acc])
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
						dico_accession_line[acc] = dico_accession_infile[chr][acc].readline().split()
				
				#4bb- recording position informations
				for acc in LIST_ACC:
					if dico_accession_line[acc] == []:
						dico_accession_positions[acc] = 1000000000000000
					else:
						dico_accession_positions[acc] = int(dico_accession_line[acc][1])
					# print acc, dico_accession_line[acc]
				min_position = get_min_acc_pos(dico_accession_positions)
				# print min_position
				
				if min_position == 1000000000000000:
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
					liste2print = [chr, str(min_position), '.', allele_to_keep[0], ','.join(allele_to_keep[1:]), '.', '.', '.', 'GT:AD:DP']
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
							genotype = genotype_accession(liste_cov_int, allele_to_keep, 0.01, DICO_PLOIDY[acc])
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
			dico_accession_infile[chr][acc].close()
		
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
			sys.exit('This is embarrassing... The program exited without finishing because a gene without CDS and exon was found:')

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

def run_qsub(QUEUE, COMMANDE_LIST, PROC, JID, LMEM):
	
	list_job = []
	for n in COMMANDE_LIST:
		if LMEM == None:
			if PROC > 1:
				sys.stdout.write("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'\n")
				sys.stdout.flush()
				qs=os.popen("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
			else:
				sys.stdout.write("qsub -q "+QUEUE+" -terse -b yes -V -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'\n")
				sys.stdout.flush()
				qs=os.popen("qsub -q "+QUEUE+" -terse -b yes -V -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
		else:
			if PROC > 1:
				sys.stdout.write("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'\n")
				sys.stdout.flush()
				qs=os.popen("qsub -q "+QUEUE+" -pe parallel_smp "+str(PROC)+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
			else:
				sys.stdout.write("qsub -q "+QUEUE+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'\n")
				sys.stdout.flush()
				qs=os.popen("qsub -q "+QUEUE+" -terse -b yes -V -l mem_free="+LMEM+" -N "+JID+" -o RnaSeqlog.txt -e log_error.txt '"+n+"'")
				list_job.append(job_ID2list(qs))
			
	hold_job(list_job)

def run_step_A(STAR, PROC, OUT_REF1, REF, SJDBO, JID, QUEUE):
	index = '%s --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s' %	(STAR, PROC, OUT_REF1, REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), index, 'Error in step A:\n')
	else:
		run_qsub(QUEUE, [index], PROC, JID+'-INDEX1', "12G")
	sys.stdout.write("Step a: reference indexation done\n")
	sys.stdout.flush()

def run_step_B(TMP, DICO_LIB, STAR, PROC, PREFIX, STAR_OPT, OUT_REF1, QUEUE):
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
			run_qsub(QUEUE, [mappP], PROC, PREFIX+'-MapP1', "12G")
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
			run_qsub(QUEUE, [mappS], PROC, PREFIX+'-MapS1', "12G")
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
						print data[0:3], data[3], dico_sites[line_id][0], data[4], dico_sites[line_id][1], data[5], dico_sites[line_id][2]
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
						print data[0:3], data[3], dico_sites[line_id][0], data[4], dico_sites[line_id][1], data[5], dico_sites[line_id][2]
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

def run_step_C(STAR, PROC, OUT_REF2, REF, SJDBO, JID, SJDBFCSE, QUEUE):
	index = '%s --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbOverhang %s --sjdbFileChrStartEnd %s' % (STAR, PROC, OUT_REF2, REF, SJDBO, SJDBFCSE)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), index, 'Error in step C:\n')
	else:
		run_qsub(QUEUE, [index], PROC, JID+'-INDEX2', "12G")
	sys.stdout.write("Step c: reference indexation done\n")
	sys.stdout.flush()

def run_step_D (LIBS, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, QUEUE):
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
			run_qsub(QUEUE, [mapp2], PROC, ACC+'-Map2', "12G")
		# checking step
		if not(os.path.isfile(ACC+'/'+lib+'Aligned.out.sam')):
			return 'An error was encountered in step d. The program exited without finishing for accession '+ACC+'\n'
		os.remove(ACC+'/'+lib+'Log.progress.out')
		os.remove(ACC+'/'+lib+'Log.out')
	sys.stdout.write("Step d: Mapping done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_E(LIB, ACC, PREFIX, JAVA, PICARD, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	to_merge = ''
	for n in LIB:
		to_merge = to_merge + 'INPUT='+ACC+'/'+n+'Aligned.out.sam '
	merge = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s MergeSamFiles %s OUTPUT=%s MERGE_SEQUENCE_DICTIONARIES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=%s' % (JAVA, PICARD, to_merge, ACC+'/'+ACC+'_merged.bam', ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), merge, 'Error in step E:\n')
	else:
		run_qsub(QUEUE, [merge], 1, ACC+'-Merge', "12G")
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_merged.bam')):
		return 'An error was encountered in step e. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	for n in LIB:
		os.remove(ACC+'/'+n+'Aligned.out.sam')
	sys.stdout.write("Step e: Merging done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_F(ACC, JAVA, PICARD, PREFIX, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	rmdup = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true QUIET=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VERBOSITY=WARNING TMP_DIR=%s' % (JAVA, PICARD, ACC+'/'+ACC+'_merged.bam', ACC+'/'+ACC+'_rmdup.bam', ACC+'/'+ACC+'_duplicate', ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), rmdup, 'Error in step F:\n')
	else:
		run_qsub(QUEUE, [rmdup], 1, ACC+'-RmDup', "12G")
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_rmdup.bam')):
		return 'An error was encountered in step f. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	os.remove(ACC+'/'+ACC+'_merged.bam')
	os.remove(ACC+'/'+ACC+'_merged.bai')
	sys.stdout.write("Step f: duplicate removal done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_G(JAVA, PICARD, ACC, REF, PREFIX, QUEUE):
	TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
	reorder = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s CREATE_INDEX=true TMP_DIR=%s' % (JAVA, PICARD, ACC+'/'+ACC+'_rmdup.bam',  ACC+'/'+ACC+'_reorder.bam', REF, ACC+'/'+TMP)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), reorder, 'Error in step G:\n')
	else:
		run_qsub(QUEUE, [reorder], 1, ACC+'-Reord', "12G")
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_reorder.bam')):
		return 'An error was encountered in step g. The program exited without finishing for accession '+ACC+'\n'
	shutil.rmtree(ACC+'/'+TMP)
	os.remove(ACC+'/'+ACC+'_rmdup.bam')
	sys.stdout.write("Step g: Bam reodering done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_H(JAVA, GATK, ACC, REF, PREFIX, QUEUE):
	trim = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T SplitNCigarReads -I %s -o %s -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -R %s -U ALLOW_N_CIGAR_READS' % (JAVA, GATK, ACC+'/'+ACC+'_reorder.bam', ACC+'/'+ACC+'_trim.bam', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), trim, 'Error in step H:\n')
	else:
		run_qsub(QUEUE, [trim], 1, ACC+'-Trim', "12G")
	# checking step
	if not(os.path.isfile(ACC+'/'+ACC+'_trim.bam')):
		return 'An error was encountered in step h. The program exited without finishing for accession '+ACC+'\n'
	os.remove(ACC+'/'+ACC+'_reorder.bam')
	os.remove(ACC+'/'+ACC+'_reorder.bai')
	sys.stdout.write("Step h: Read splitting done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_I(ACC, JAVA, GATK, REF, PLOIDY, PREFIX, UseUnifiedGenotyperForBaseRecal, QUEUE):
	
	realTC = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T RealignerTargetCreator -o %s -I %s -R %s' % (JAVA, GATK, ACC+'/'+ACC+'_RTC.intervals', ACC+'/'+ACC+'_trim.bam', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), realTC, 'Error in step I:\n')
	else:
		run_qsub(QUEUE, [realTC], 1, ACC+'-realTC', "12G")
	if not(os.path.isfile(ACC+'/'+ACC+'_RTC.intervals')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'

	Ireal = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s -T IndelRealigner -o %s -I %s -targetIntervals %s -R %s' % (JAVA, GATK, ACC+'/'+ACC+'_realigned.bam', ACC+'/'+ACC+'_trim.bam', ACC+'/'+ACC+'_RTC.intervals', REF)
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), Ireal, 'Error in step I:\n')
	else:
		run_qsub(QUEUE, [Ireal], 1, ACC+'-Ireal', "12G")
	if not(os.path.isfile(ACC+'/'+ACC+'_realigned.bam')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'
	os.remove(ACC+'/'+ACC+'_RTC.intervals')
	
	if not(os.path.isfile(ACC+'/'+ACC+'_realigned.bam')):
		return 'An error was encountered in step i. The program exited without finishing for accession '+ACC+'\n'
	sys.stdout.write("Step i: Base recalibration done for accession "+ACC+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_J(ACC_ID, ALLELE_COUNT, PYTHON, REF, QUEUE):

	sys.stdout.write('Working on '+ACC_ID+' accession\n')
	
	# calculating sequence length
	dico_chr = {}
	sequence_dict = SeqIO.index(REF, "fasta")
	for n in sequence_dict:
		dico_chr[n] = [len(str(sequence_dict[n].seq)), str(sequence_dict[n].seq)]
	del sequence_dict
	
	#######################
	#6 pseudo GVCF generation
	
	if os.path.isfile(ACC_ID+'/'+ACC_ID+'_real_recal.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_real_recal.bam'
	elif os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam'):
		bam = ACC_ID+'/'+ACC_ID+'_realigned.bam'
		sys.stdout.write('Warning no recalibrated bam file found for '+ACC_ID+'... Continue using realigned bam instead.\n')
		sys.stdout.flush()
	else:
		return 'ERROR : Neither '+ACC_ID+'/'+ACC_ID+'_real_recal.bam or '+ACC_ID+'/'+ACC_ID+'_realigned.bam were found.\n'
	
	# Running allele count in bam
	sys.stdout.write('Initiating allele count\n')
	allele_count = '%s %s -r %s -b %s -o %s' % (PYTHON, ALLELE_COUNT, REF, bam, ACC_ID+'/'+ACC_ID+'_allele_count')
	if QUEUE == None:
		run_job(getframeinfo(currentframe()), allele_count, 'Error in step J:\n')
	else:
		run_qsub(QUEUE, [allele_count], 1, ACC_ID+'-AlCount', "12G")
	
	# checking step
	for n in dico_chr:
		if not(os.path.isfile(ACC_ID+'/'+ACC_ID+'_allele_count_'+n+'.gz')):
			return 'ERROR : File '+ACC_ID+'/'+ACC_ID+'_allele_count_'+n+'.gz'+' as not been found.\n'
	sys.stdout.write("Step j: Allele count done for accession "+ACC_ID+"\n")
	sys.stdout.flush()
	
	return 0

def run_step_L(GFF3, REF, LOCA_PROGRAMS, PREFIX, DICO_LIB):
	# recording gff3 informations
	sys.stdout.write('recording gff3 informations\n')
	dico_gene_gff3 = {}
	Record_CDS_and_UTR(GFF3, dico_gene_gff3)
	
	# calculating exon coverage by accession
	sys.stdout.write('calculating exon coverage by accession\n')
	list_lib = []
	for ACC_ID in DICO_LIB:
		sys.stdout.write(ACC_ID+'\n')
		sys.stdout.flush()
		if os.path.isfile(ACC_ID+'/'+ACC_ID+'_real_recal.bam'):
			bam = ACC_ID+'/'+ACC_ID+'_real_recal.bam'
		elif os.path.isfile(ACC_ID+'/'+ACC_ID+'_realigned.bam'):
			bam = ACC_ID+'/'+ACC_ID+'_realigned.bam'
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
	
	# for acc in list_lib:
		# os.remove(acc+'/'+acc+'_realigned_recalibrated.bam.cov')

def run_stat_step_D(PREFIX, CONFIG):
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

def run_stat_step_F(PREFIX, DICO_LIB):
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

def run_D_to_J(STEPS, LIB, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, JAVA, PICARD, REF, GATK, PLOIDY, UseUnifiedGenotyperForBaseRecal, BAMTOOLS, ALLELE_COUNT, PYTHON, QUEUE):
	if 'd' in STEPS:
		to_return = run_step_D(LIB, ACC, STAR, OUT_REF2, PROC, STAR_OPT, PREFIX, QUEUE)
	if 'e' in STEPS:
		to_return = run_step_E(LIB, ACC, PREFIX, JAVA, PICARD, QUEUE)
	if 'f' in STEPS:
		to_return = run_step_F(ACC, JAVA, PICARD, PREFIX, QUEUE)
	if 'g' in STEPS:
		to_return = run_step_G(JAVA, PICARD, ACC, REF, PREFIX, QUEUE)
	if 'h' in STEPS:
		to_return = run_step_H(JAVA, GATK, ACC, REF, PREFIX, QUEUE)
	if 'i' in STEPS:
		to_return = run_step_I(ACC, JAVA, GATK, REF, PLOIDY, PREFIX, UseUnifiedGenotyperForBaseRecal, QUEUE)
	if 'j' in STEPS:
		to_return = run_step_J(ACC, ALLELE_COUNT, PYTHON, REF, QUEUE)
	if to_return == 0:
		return 0
	else:
		return to_return

def main_run_analysis(job):

	try:
		rslt = run_D_to_J(job[0],job[1],job[2],job[3],job[4],job[5],job[6],job[7],job[8],job[9],job[10],job[11],job[12],job[13],job[14],job[15],job[16], job[17])
	except Exception as e:
		print e
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
	'l: Coverage statistics\t')
	(options, args) = parser.parse_args()
	
	if options.conf == None:
		sys.exit('--conf argument is missing')
	if options.steps == None:
		sys.exit('--septs argument is missing')
	
	#Loading the file locating programs
	PATHNAME = os.path.dirname(sys.argv[0])
	LOCA_PROGRAMS = ConfigParser.RawConfigParser()
	LOCA_PROGRAMS.read(PATHNAME+'/loca_programs.conf')
	
	#Loading the configuration file
	config = ConfigParser.RawConfigParser()
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
	ALLELE_COUNT = LOCA_PROGRAMS.get('Programs','allele_count')
	PYTHON = LOCA_PROGRAMS.get('Programs','python')
	QUEUE = options.queue
	
	UseUnifiedGenotyperForBaseRecal =  'no'
	if config.has_section('General'):
		if config.has_option('General', 'UseUnifiedGenotyperForBaseRecal'):
			UseUnifiedGenotyperForBaseRecal =  config.get('General', 'UseUnifiedGenotyperForBaseRecal')
	
	if not(os.path.isfile('.'.join(REF.split('.')[:-1])+'.dict')):
		sys.stdout.write('No .dict file found for the reference. Running indexation...')
		create_dict = '%s -XX:ParallelGCThreads=1 -Xmx8G -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s' % (JAVA, PICARD, REF, '.'.join(REF.split('.')[:-1])+'.dict')
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), create_dict, 'Error when creating the .dict file:\n')
		else:
			run_qsub(QUEUE, [create_dict], 1, PREFIX+'-Dict', "12G")
		sys.stdout.write('Done\n')
	
	if not(os.path.isfile(REF+'.fai')):
		sys.stdout.write('No .fai file found for the reference. Running indexation...')
		create_index = '%s faidx %s' % (SAMTOOLS, REF)
		if QUEUE == None:
			run_job(getframeinfo(currentframe()), create_index, 'Error when creating the .fai file:\n')
		else:
			run_qsub(QUEUE, [create_index], 1, PREFIX+'-Index', "12G")
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
	
	# Reference indexation
	if 'a' in options.steps:
		if os.path.isdir(OUT_REF1):
			shutil.rmtree(OUT_REF1)
		os.mkdir(OUT_REF1)
		run_step_A(STAR, PROC, OUT_REF1, REF, SJDBO, PREFIX, QUEUE)
	

	# Aligning reads for junction identification
	if 'b' in options.steps:
		TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
		run_step_B(TMP, dico_lib, STAR, PROC, PREFIX, STAR_OPT, OUT_REF1, QUEUE)
		
	# Running reference indexation with intron/exon junction
	if 'c' in options.steps:
		if os.path.isdir(OUT_REF2):
			shutil.rmtree(OUT_REF2)
		os.mkdir(OUT_REF2)
		run_step_C(STAR, PROC, OUT_REF2, REF, SJDBO, PREFIX, SJDBFCSE, QUEUE)
	
	# Running step d to j
	if 'd' in options.steps or 'e' in options.steps or 'f' in options.steps or 'g' in options.steps or 'h' in options.steps or 'i' in options.steps or 'j' in options.steps:
		# threads = []
		listJobs = []
		for acc in dico_lib:
			sys.stdout.flush()
			listJobs.append([options.steps, dico_lib[acc], acc, STAR, OUT_REF2, 4, STAR_OPT, PREFIX, JAVA, PICARD, REF, GATK, dico_ploidy[acc],UseUnifiedGenotyperForBaseRecal, BAMTOOLS, ALLELE_COUNT, PYTHON, QUEUE])
		
		pool = mp.Pool(processes=int(PROC))
		results = pool.map(main_run_analysis, listJobs)
		
		if os.path.isdir(PREFIX+'_ref_star_2'):
			shutil.rmtree(PREFIX+'_ref_star_2')
		
		for n in results:
			if n != 0:
				sys.stdout.write(str(n)+'\n')
		
	# GVCF accession merging
	if 'k' in options.steps:
		create_pseudo_VCF(dico_lib.keys(), REF, PREFIX, dico_ploidy)
	
	# calculating exon coverage
	if 'l' in options.steps:
		GFF3 = config.get('General', 'gff3')
		run_step_L(GFF3, REF, LOCA_PROGRAMS, PREFIX, dico_lib)
	
	# Collecting statistics by accession
	if 'd' in options.steps:
		if not(os.path.isdir(PREFIX)):
			os.mkdir(PREFIX)
		run_stat_step_D(PREFIX, config)
	
	if 'f' in options.steps:
		if not(os.path.isdir(PREFIX)):
			os.mkdir(PREFIX)
		run_stat_step_F(PREFIX, dico_lib)

	
if __name__ == "__main__": __main__()
