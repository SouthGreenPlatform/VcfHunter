import argparse
import sys
import gzip


def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser()
	# Wrapper options.
	parser.add_argument( '-v', '--vcf',			dest='vcf',			default=None,		help='The vcf file')
	parser.add_argument( '-m', '--matrix',		dest='matrix',		default=None,		help='The matrix file')
	parser.add_argument( '-M', '--marker',		dest='marker',		default=None,		help='A Marker file name. It will phase only markers (if possible) in this file in addition to those provided in the matrix.')
	parser.add_argument( '-o', '--output',		dest='output',		default=None,		help='The output file name')
	options = parser.parse_args()
	
	VCF = options.vcf
	MAT = options.matrix
	OUT = options.output
	MARK = options.marker
	
	# Recording additional markers to phase (They should be homozygous)
	if MARK == None:
		MarkerNameSet = set()
	else:
		MarkerNameSet = set()
		file = open(MARK)
		for line in file:
			data = line.split()
			if data:
				MarkerNameSet.add(data[0])
		file.close()
	
	# Recording accessions the matrix and header
	DicoMat = {}
	file = open(MAT)
	Mheader = file.readline().split()
	DicoMatIndex = {}
	accessions = Mheader[4:]
	for acc in accessions:
		DicoMatIndex[acc] = Mheader.index(acc)
	# A piece of verification
	if Mheader[3] != 'rephased':
		sys.exit('Oups, the program exited without finishing... The format of the matrix file is not good!')
	for line in file:
		data = line.split()
		DicoMat[data[0]] = data
	
	# Reading vcf and calculating correspondence between genotype coding and alleles
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	outfile = open(OUT,'w')
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				IC = header.index('#CHROM')
				IP = header.index('POS')
				IREF = header.index('REF')
				IALT = header.index('ALT')
				IFOR = header.index('FORMAT')
				DIACC = {}
				for acc in accessions:
					DIACC[acc] = header.index(acc)
				# print(IC, IP, IREF, IALT, IFOR, DIACC)
				AccNumber = len(DIACC)
			# Printing headers
			elif data[0][0] == '#':
				pass
			# Working on variant
			else:
				DicoCor = {}
				SetGeno = set()
				chr =  data[IC]
				pos =  data[IP]
				MarkerName = chr+'M'+pos
				if MarkerName in DicoMat:
					ref_allele = data[IREF]
					alt_allele = data[IALT].split(',')
					AllAlleles = [ref_allele] + alt_allele
					flag_format = data[IFOR].split(':')
					for acc in DIACC:
						ACCINFO = data[DIACC[acc]].split(':')
						GT = ':'.join(sorted(list(set(ACCINFO[data[IFOR].split(':').index('GT')].replace('|','/').split('/')))))
						code = DicoMat[MarkerName][DicoMatIndex[acc]]
						if not('.') in GT:
							SetGeno.add(GT)
							if code != '--':
								if not(code in DicoCor):
									DicoCor[code] = []
								DicoCor[code].append(GT)
					DicoBest = {}
					for codes in DicoCor:
						DicoBest[codes] = []
						DicoBest[codes].append([])
						DicoBest[codes].append([])
						for geno in SetGeno:
							Value = DicoCor[codes].count(geno)
							DicoBest[codes][0].append(geno)
							DicoBest[codes][1].append(Value)
					DicoSelect = {}
					OK = True
					for codes in DicoBest:
						MaxValue =  max(DicoBest[codes][1])
						if DicoBest[codes][1].count(MaxValue) == 1:
							DicoSelect[codes] = DicoBest[codes][0][DicoBest[codes][1].index(MaxValue)]
							# print(codes, DicoSelect[codes], MaxValue/sum(DicoBest[codes][1]))
						else:
							OK = False
							sys.stdout.write('Oups, their is more than 1 maximal value for code'+codes+'. The marker '+MarkerName+' will not be reported.\n')
							sys.stdout.write(str(DicoBest)+'\n')
					if OK:
						FinalOK = True
						VerifSet = set()
						for n in DicoSelect:
							if not(DicoSelect[n] in VerifSet):
								VerifSet.add(DicoSelect[n])
							else:
								FinalOK = False
								sys.stdout.write('Oups, two phases have the same genotype. This is not a good marker... The marker '+MarkerName+' will not be reported.\n')
								for k in DicoSelect:
									sys.stdout.write(str(DicoBest)+'\n')
									sys.stdout.write(str(k)+'\t'+str(DicoSelect[n])+'\n')
						if FinalOK:
							if 'nn' in DicoSelect:
								if len(DicoSelect) != 2:
									sys.exit('Oups, their is a bug. Their is more than nn and np coding in the line for marker: '+MarkerName)
								else:
									nnGeno = set(DicoSelect['nn'].split(':'))
									npGeno = set(DicoSelect['np'].split(':'))
									if len(nnGeno) == 1 and len(npGeno) == 2:
										haplo1 = AllAlleles[int(list(nnGeno)[0])]
										Difference = npGeno - nnGeno
										if len(Difference) != 1:
											sys.exit('Oups, their is a bug. There is more than 2 alleles identified for marker: '+MarkerName)
										haplo2 = AllAlleles[int(list(Difference)[0])]
										# print(MarkerName, DicoSelect, AllAlleles, haplo1, haplo2, nnGeno, npGeno, Difference)
									elif len(nnGeno) == 2 and len(npGeno) == 1:
										haplo2 = AllAlleles[int(list(npGeno)[0])]
										Difference = nnGeno - npGeno
										if len(Difference) != 1:
											sys.exit('Oups, their is a bug. There is more than 2 alleles identified for marker: '+MarkerName)
										haplo1 = AllAlleles[int(list(Difference)[0])]
										# print(MarkerName, DicoSelect, AllAlleles, haplo1, haplo2, nnGeno, npGeno, Difference)
									else:
										sys.exit('Oups, their is a bug. There is more than 2 alleles identified for marker: '+MarkerName)
									outfile.write(MarkerName+'\t'+haplo1+'|'+haplo2+'\n')
									outfile.flush()
							elif 'hk' in DicoSelect:
								if len(DicoSelect) != 3:
									# print(MarkerName, DicoSelect, len(DicoSelect))
									if len(DicoSelect) >= 3:
										sys.exit('Oups, their is a bug. Their is more than hh, hk and kk coding in the line for marker: '+MarkerName+'\n'+str(DicoSelect)+'\n')
									else:
										sys.stdout.write('Oups, a marker has not all possible combinations. This is not allowed... The marker '+MarkerName+' will not be reported.\n')
								else:
									hhGeno = set(DicoSelect['hh'].split(':'))
									hkGeno = set(DicoSelect['hk'].split(':'))
									kkGeno = set(DicoSelect['kk'].split(':'))
									if len(hhGeno) == 1 and len(kkGeno) == 1 and len(hkGeno) == 2:
										haplo1 = AllAlleles[int(list(hhGeno)[0])]
										haplo2 = AllAlleles[int(list(kkGeno)[0])]
										if haplo1 == haplo2:
											sys.exit('Oups, their is a bug. Both haplotypes have the same allele for marker: '+MarkerName)
										# print(MarkerName, DicoSelect, AllAlleles, haplo1, haplo2, hhGeno, kkGeno)
									else:
										sys.exit('Oups, their is a bug. There is more than 2 alleles identified for marker: '+MarkerName)
									outfile.write(MarkerName+'\t'+haplo1+'|'+haplo2+'\n')
									outfile.flush()
				elif MarkerName in MarkerNameSet:
					ref_allele = data[IREF]
					alt_allele = data[IALT].split(',')
					AllAlleles = [ref_allele] + alt_allele
					flag_format = data[IFOR].split(':')
					for acc in DIACC:
						ACCINFO = data[DIACC[acc]].split(':')
						GT = ':'.join(sorted(list(set(ACCINFO[data[IFOR].split(':').index('GT')].replace('|','/').split('/')))))
						if not(GT in DicoCor):
							DicoCor[GT] = 0
						DicoCor[GT] += 1
					
					TotalSansMissing = 0 + AccNumber
					OK = True
					if '.' in DicoCor:
						TotalSansMissing -= DicoCor['.']
						if DicoCor['.']/AccNumber > 0.2:
							sys.stdout.write('Oups, to much missing data. The marker '+MarkerName+' will not be reported.\n')
							OK = False
						else:
							del DicoCor['.']
					if OK:
						RetainedGenotypes = set()
						for geno in DicoCor:
							if DicoCor[geno]/TotalSansMissing >= 0.05:
								RetainedGenotypes.add(geno)
						ListRetainedGenotypes = list(RetainedGenotypes)
						if len(ListRetainedGenotypes) == 1: # One genotype found in all our accessions according to our criteria
							if len(ListRetainedGenotypes[0]) == 1: # homozygous according to our criteria
								haplo1 = AllAlleles[int(ListRetainedGenotypes[0])]
								haplo2 = AllAlleles[int(ListRetainedGenotypes[0])]
								# print(MarkerName, DicoCor, TotalSansMissing, ListRetainedGenotypes, AllAlleles, haplo1, haplo2)
								outfile.write(MarkerName+'\t'+haplo1+'|'+haplo2+'\n')
								outfile.flush()
						elif len(ListRetainedGenotypes) == 2: #Two genotypes founds according to our criteria so one parent is heterozygous. It is likely the other one but because we are not sure, we will not report this line
							geno1 = set(ListRetainedGenotypes[0].split(":"))
							geno2 = set(ListRetainedGenotypes[1].split(":"))
							if len(geno1 & geno2) == 1 and len(geno1 | geno2) == 2:
								haplo1 = AllAlleles[int(list(geno1 & geno2)[0])]
								haplo2 = AllAlleles[int(list(geno1 & geno2)[0])]
								# print(MarkerName, DicoCor, TotalSansMissing, ListRetainedGenotypes, AllAlleles, haplo1, haplo2)
								# outfile.write(MarkerName+'\t'+haplo1+'|'+haplo2+'\n')
								# outfile.flush()
								
						else: # we are probably heterozygous in both parents and thus a phasing step is needed
							pass
							
					
				
				
	outfile.close()
	
if __name__ == "__main__": __main__()
