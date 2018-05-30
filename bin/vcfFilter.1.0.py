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
import os
import datetime
import gzip

sys.stdout.write('modules loaded\n')


def filter_vcf(VCF, NAMES, OUTGROUP, PREFIX, RMTYPE, MINCOV, MINAL, NMISS, RMALALT, MAXCOV, MINFREQ, OUTGZIP):
	
	"""
		Filter vcf file and output a filtered vcf
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param NAMES: Path to file containing accession to treat
		:type NAMES: str
		:param OUTGROUP: Path to file containing accession not to consider for filtering but to print in the output
		:type OUTGROUP: str
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
		:param RMTYPE: Type of variant to remove
		:type RMTYPE: str
		:param MINCOV: Minimal coverage to keep a genotype
		:type MINCOV: int
		:param MINAL: Minimal allele coverage to keep a genotype
		:type MINAL: int
		:param NMISS: Maximal missing genotype to keep a variant site
		:type NMISS: int
		:param RMALALT: Number of alternative variant to remove the variant site.
		:type RMALALT: str
		:return: A filtered vcf file
		:rtype: void
		:param MAXCOV: Maximal coverage to keep a genotype
		:type MAXCOV: int
		:param MINFREQ: Maximal coverage to keep a genotype
		:type MINFREQ: float
		:param MINFREQ: Deciding if output should be gziped
		:type MINFREQ: str
	"""
	nb_remove = 0
	nb_kept = 0
	nb_autapo_remove = 0
	nb_allele_remove = 0
	nb_missing_remove = 0
	nb_tag_remove = 0
	nb_INDEL_remove = 0
	nb_SNP_remove = 0
	nb_BadFormat = 0
	
	# Recording accession to treate
	file = open(NAMES)
	DICO_NAME = set()
	for line in file:
		data = line.split()
		if data:
			DICO_NAME.add(data[0])
	file.close()
	
	
	# Recording outgroup
	if OUTGROUP == None:
		DICO_OUTGROUP = set()
	else:
		file = open(OUTGROUP)
		DICO_OUTGROUP = set()
		for line in file:
			data = line.split()
			if data:
				DICO_OUTGROUP.add(data[0])
		file.close()
	
	# recording variant status to exclude
	if RMTYPE == None:
		exclude = []
	else:
		exclude = RMTYPE.split(',')
	
	# recording variant alternatives to exclude
	if  RMALALT == None:
		exclude_al = []
	else:
		exclude_al = [int(i) for i in RMALALT.split(':')]
	
	# Creating output file
	if OUTGZIP == 'n':
		outfile = open(PREFIX+'_filt.vcf','w')
	elif OUTGZIP == 'y':
		outfile = gzip.open(PREFIX+'_filt.vcf.gz','wt')
	else:
		sys.exit('Wrong argument passed to --outgzip options. Argument accepted: y or n\n')
	
	# Reading vcf file
	PrintFilter = 1
	
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
	for line in file:
		data = line.split()
		if data:
			if data[0][0:8] == "##contig" and PrintFilter:
				if exclude:
					outfile.write('##Additionnal.filter=<ID=TAGRemoval,Date='+str(datetime.datetime.now())+',Description="Variant'+' '.join(exclude)+' are removed">\n')
				outfile.write('##Additionnal.filter=<ID=CoverageFiltration,Date='+str(datetime.datetime.now())+',Description="Genotype having less than '+str(MINCOV)+' x coverage and less than '+str(MINAL)+' x coverage for each allele are converted to missing">\n')
				outfile.write('##Additionnal.filter=<ID=MissingDataFiltration,Date='+str(datetime.datetime.now())+',Description="SNP with more than '+str(NMISS)+' genotype missing are removed">\n')
				PrintFilter = 0
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				list_to_print_in_output = []
				for n in header:
					list_to_print_in_output.append(n)
					if n == 'FORMAT':
						break
				for n in header:
					if (n in DICO_OUTGROUP) or (n in DICO_NAME):
						list_to_print_in_output.append(n)
				outfile.write('\t'.join(list_to_print_in_output)+'\n')
				FILTERpos = header.index('FILTER')
				REFpos = header.index('REF')
				ALTpos = header.index('ALT')
				FLAGformat = header.index('FORMAT')
				Accession_start = header.index('FORMAT')+1
			# Printing headers
			elif data[0][0] == '#':
				outfile.write(line)
			# Working on variant
			else:
				nb_alt = set()
				remove = 0
				# Identification we should remove the variant based on the FILTER tag
				filter_tag = data[FILTERpos].split(';')
				for n in filter_tag:
					if n in exclude:
						remove = 1
				if remove == 1:
					nb_tag_remove += 1
				# Initialisation of the corrected line
				list_to_print_in_output = data[0:Accession_start]
				# Identification genotype to convert to missing due to MINCOV and MINAL parameters
				nb_missing = 0
				ref_allele = data[REFpos]
				alt_allele = data[ALTpos].split(',')
				flag_format = data[FLAGformat].split(':')
				if 'AD' in flag_format and 'DP' in flag_format and 'GT' in flag_format:
					GTpos = flag_format.index('GT')
					DPpos = flag_format.index('DP')
					ADpos = flag_format.index('AD')
					dico_autapo = {}
					for accession in header:
						# working on accession to treat
						if accession in DICO_NAME:
							filtered_accession = filter_accession(accession, data, header, flag_format, GTpos, DPpos, ADpos, MINCOV, MINAL, MAXCOV, MINFREQ)
							nb_alt.update(filtered_accession[2])
							list_to_print_in_output.append(filtered_accession[0])
							if filtered_accession[1]:
								nb_missing += 1
							for z in filtered_accession[2]:
								if z in dico_autapo:
									dico_autapo[z] += 1
								else:
									dico_autapo[z] = 1
						# working on outgroup (not taken in acount for variant filtering)
						elif accession in DICO_OUTGROUP:
							list_to_print_in_output.append(filter_accession(accession, data, header, flag_format, GTpos, DPpos, ADpos, MINCOV, MINAL, MAXCOV, MINFREQ)[0])
					# Validating there is not to much missing data
					if nb_missing > NMISS:
						remove = 1
						nb_missing_remove += 1
					if 'AUTAPO' in exclude:
						# Validating if it is an autapomorphy
						if cherche_autapo(dico_autapo):
							remove = 1
							nb_autapo_remove += 1
					if 'SNP' in exclude:
						# Validating if it is a SNP alone
						dic_snp = set(alt_allele)
						dic_snp.add(ref_allele)
						if IsSNP(dic_snp):
							remove = 1
							nb_SNP_remove += 1
					if 'INDELS' in exclude:
						# Validating if it is a SNP alone
						dic_snp = set(alt_allele)
						dic_snp.add(ref_allele)
						if IsIndel(dic_snp):
							remove = 1
							nb_INDEL_remove += 1
					nb_alt.discard('.')
					if len(nb_alt) in exclude_al:
						remove = 1
						nb_allele_remove += 1
					if remove:
						nb_remove += 1
					else:
						nb_kept += 1
						outfile.write('\t'.join(list_to_print_in_output)+'\n')
				else:
					nb_BadFormat += 1
					nb_remove += 1
	sys.stdout.write('Removed variant: '+str(nb_remove)+'\n')
	sys.stdout.write('\tRemoved variant (Bad format): '+str(nb_BadFormat)+'\n')
	sys.stdout.write('\tRemoved variant (missing): '+str(nb_missing_remove)+'\n')
	sys.stdout.write('\tRemoved variant (tag): '+str(nb_tag_remove)+'\n')
	sys.stdout.write('\tRemoved variant (autapomorphy): '+str(nb_autapo_remove)+'\n')
	sys.stdout.write('\tRemoved variant (SNP): '+str(nb_SNP_remove)+'\n')
	sys.stdout.write('\tRemoved variant (INDEL): '+str(nb_INDEL_remove)+'\n')
	sys.stdout.write('\tRemoved variant (bad allele number): '+str(nb_allele_remove)+'\n')
	sys.stdout.write('Kept variant: '+str(nb_kept)+'\n')

def filter_accession(ACCESSION, DATA, HEADER, FLAG_FORMAT, GTPOS, DPPOS, ADPOS, MINCOV, MINAL, MAXCOV, MINFREQ):
	
	"""
		Filter accession based on MINCOV and MINAL parameters
		
		:param ACCESSION: Accession to treat.
		:type ACCESSION: str
		:param DATA: A list corresponding to a vcf line.
		:type DATA: list
		:param HEADER: A list corresponding to the vcf header line.
		:type HEADER: list
		:param MINCOV: Minimal coverage to keep a genotype.
		:type MINCOV: int
		:param MINAL: Minimal allele coverage to keep a genotype.
		:type MINAL: int
		:param FLAG_FORMAT: FORMAT of the information of the accession.
		:type FLAG_FORMAT: int
		:param GTPOS: Position of the GT tag in the data information of the accession.
		:type GTPOS: int
		:param DPPOS: Position of the DP tag in the data information of the accession.
		:type DPPOS: int
		:param ADPOS: Position of the AD tag in the data information of the accession.
		:type ADPOS: int
		:return: A list with [0] --> gentype calling recalculated, [1] --> missing genotype (bolean), [2] --> set listing alleles found in the genotype
		:rtype: list
		:param MAXCOV: Maximal coverage to keep a genotype
		:type MAXCOV: int
		:param MINFREQ: Maximal coverage to keep a genotype
		:type MINFREQ: float
	"""
	
	ACCinfo = DATA[HEADER.index(ACCESSION)]
	
	missing = False
	to_print = ''
	accession_var = ACCinfo.split(':')
	# Looking if a genotype has been called
	geno = accession_var[GTPOS].replace('/','|').split('|')
	if '.' in geno:
		missing = True
		to_print = ACCinfo
		forautapo = set()
	else:
		# Looking if site coverage is sufficient
		convert_to_missing = False
		if DPPOS >= len(accession_var):
			dp_cov = -1
		elif accession_var[DPPOS] == '.':
			dp_cov = -1
		else:
			dp_cov = int(accession_var[DPPOS])
		# Not enough coverage
		if dp_cov < MINCOV:
			convert_to_missing = True
		elif dp_cov > MAXCOV:
			convert_to_missing = True
		# Enough coverage
		else:
			ad_cov = accession_var[ADPOS].split(',')
			for n in geno:
				if int(ad_cov[int(n)]) < MINAL:
					convert_to_missing = True
				if int(ad_cov[int(n)])/float(dp_cov)< MINFREQ:
					convert_to_missing = True
		# The genotype should be converted to missing
		if convert_to_missing:
			missing = True
			missing_geno = []
			for n in geno:
				missing_geno.append('.')
			accession_var[GTPOS] = '/'.join(missing_geno)
			to_print = ':'.join(accession_var)
			forautapo = set()
		else:
			to_print = ACCinfo
			# Recording alleles found for autapomorphy search
			forautapo = set(geno)
	
	return [to_print, missing, forautapo]

def IsSNP(DICO):
	
	"""
		Identify if it is a SNP alone
		
		:param DICO: A dictionary of alleles.
		:type DICO: dictionary
		:return: A bolean with true if it is a SNP variant site only
		:rtype: bolean
	"""
	
	max_size = 0
	deleted_allele = False
	for n in DICO:
		max_size = max(len(n),max_size)
		if n == "*":
			deleted_allele = True
	
	if max_size == 1 and deleted_allele == False:
		return True
	else:
		return False

def IsIndel(DICO):
	
	"""
		Identify if it is an indel
		
		:param DICO: A dictionary of alleles.
		:type DICO: dictionary
		:return: A bolean with true if it is a SNP variant site only
		:rtype: bolean
	"""
	
	max_size = 0
	for n in DICO:
		max_size = max(len(n),max_size)
	
	if max_size == 1:
		return False
	else:
		return True

def cherche_autapo(DICO):
	
	"""
		Identify autapomorphies in the set
		
		:param DICO: A dictionary counting the number of accession bearing a given allele.
		:type DICO: dictionary
		:return: A bolean with true if the variant is an autapomorphy and false if not
		:rtype: bolean
	"""
	
	trouve = 0
	for n in DICO:
		if DICO[n] > 1:
			trouve += 1
	if trouve > 1:
		return 0
	else:
		return 1

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '',	'--names',			dest='names',		default=None,			help='A one column file containing accession names to treat. [Default: %default]')
	parser.add_option( '',	'--outgroup',		dest='outgroup',	default=None,			help='Path to a one column file containing accession names not to consider for filtering but to print in the output [Default: %default]')	
	parser.add_option( '',	'--RmType',			dest='RmType',		default=None,			help='Variant status to filter out (several values can be passed, in this case they should be separated by ","). Values: PASS, DP_FILTER, QD_FILTER, SnpCluster, INDELS, SNP, AUTAPO [Default: %default]')
	parser.add_option( '',	'--RmAlAlt',		dest='RmAlAlt',		default=None,			help='Number of alleles at the site to remove the variant site (several values can be passed and should be sepatated by :). Values: 1,2,3,...,n [Default: %default]')
	parser.add_option( '',	'--MinCov',			dest='MinCov',		default='10',			help='Minimal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--MaxCov',			dest='MaxCov',		default='1000',			help='Maximal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--MinFreq',		dest='MinFreq',		default='0.05',			help='Minimal allele frequency to keep genotype calling (float). If the value is lower, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--MinAl',			dest='MinAl',		default='3',			help='Minimal allele coverage by accession to keep genotype calling (integer). If the value is lower for at least one allele, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--nMiss',			dest='nMiss',		default='0',			help='Maximal number of missing genotype in a line to keep the line (integer). [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')
	parser.add_option( '-g', '--outgzip',		dest='outgzip',		default='n',			help='Output files in gzip format. [Default: %default]')
	
	(options, args) = parser.parse_args()
	
	# Filtering vcf file
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.names == None:
		sys.exit('Please provide a name file to --names argument')
	filter_vcf(options.vcf, options.names, options.outgroup, options.prefix, options.RmType, int(options.MinCov), int(options.MinAl), int(options.nMiss), options.RmAlAlt, int(options.MaxCov),float(options.MinFreq), options.outgzip)
		
if __name__ == "__main__": __main__()
