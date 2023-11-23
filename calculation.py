import os, sys
from Bio import SeqIO
import pandas as pd
import math

def intersection(list1, list2):
	return list(set(list1) & set(list2))
	
def modify(gene):
	contig_name = gene.split("_")
	return "_".join(contig_name[:-1])

def calculate_score(contig_file, all_data, pipline):
	#Open Fasta file
	records = list(SeqIO.parse(contig_file, "fasta"))
	nContigs = len(records)
	
	res_columns = ["score type","nContigs", "nARG", "nMGE", "nPAT", "nARG_MGE", "nARG_PAT", "nARG_MGE_PAT", "qARG", "qMGE", "qPAT", "qARG_MGE", "qARG_PAT", \
		"qARG_MGE_PAT", "distance", "Risk Score"]
		
	def risk_score(arg_data, mge_data, pat_data):
		if arg_data.empty:
			ARG = []
		else:
			arg_temp = arg_data.copy()
			arg_temp['id'] = arg_temp['id'].apply(modify)
			ARG = arg_temp.id.unique().tolist()
		
		if mge_data.empty:
			MGE = []
		else:
			mge_temp = mge_data.copy()
			mge_temp['id'] = mge_temp['id'].apply(modify)
			MGE = mge_temp.id.unique().tolist()
		if pat_data.empty:
			PAT = []
		else:    
			PAT = pat_data.id.unique().tolist()

		#get common ones
		ARG_MGE = intersection(ARG, MGE)        
		ARG_PAT = intersection(ARG, PAT)
		ARG_MGE_PAT = intersection(ARG_MGE, PAT)

		#count contigs
		nARG = len(ARG)
		nARG_MGE = len(ARG_MGE)
		nARG_PAT = len(ARG_PAT)
		nARG_MGE_PAT = len(ARG_MGE_PAT)
		
		# normalize
		if nContigs == 0:
			fARG = 0.0
			fARG_MGE = 0.0
			fARG_PAT = 0.0
			fARG_MGE_PAT = 0.0
		else:
			fARG = float(nARG)/nContigs
			fARG_MGE = float(nARG_MGE)/nContigs
			fARG_PAT = float(nARG_PAT)/nContigs
			fARG_MGE_PAT = float(nARG_MGE_PAT)/nContigs
		
		if pipeline == 1:
			distance = math.sqrt((1 - fARG)**2 + (1 - fARG_MGE)**2 + (1 - fARG_MGE_PAT)**2)
			d = math.sqrt((1 - 0)**2 + (1 - 0)**2 + (1 - 0)**2)
		elif pipeline == 2:
			distance = math.sqrt((1 - fARG)**2 + (1 - fARG_MGE)**2 + (1 - fARG_PAT)**2 + (1 - fARG_MGE_PAT)**2)
			d = math.sqrt((1 - 0)**2 + (1 - 0)**2 + (1 - 0)**2 + (1 - 0)**2)
		score = (1 - distance/d)*10000
		
		#other statics
		nMGE = len(MGE)
		nPAT = len(PAT)
		if nContigs == 0:
			fMGE = 0.0
			fPAT = 0.0
		else:
			fMGE = float(nMGE)/nContigs
			fPAT = float(nPAT)/nContigs
	
		# make a suitable format for printing
		fARG = round(fARG, 6)
		fMGE = round(fMGE, 6)
		fPAT = round(fPAT, 6)
		fARG_MGE = round(fARG_MGE, 6)
		fARG_MGE_PAT = round(fARG_MGE_PAT, 6)
		distance = round(distance, 6)
		score = round(score, 6)

		return [nContigs, nARG, nMGE, nPAT, nARG_MGE, nARG_PAT, nARG_MGE_PAT, fARG, fMGE, fPAT, fARG_MGE, fARG_PAT, fARG_MGE_PAT, distance, score]	
		
	rs = risk_score(all_data[0], all_data[1], all_data[2], pipeline)	
	if pipline == 1:
		rs.insert(0, "Ecological")	
	elif pipline == 2:
		rs.insert(0, "Human health")		
	result = pd.DataFrame(list(zip(res_columns, rs)), columns =['Name', 'quantity']) 
	return result	

       
        
