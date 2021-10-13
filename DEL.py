#!/bin/usr/python

import sys
import os
import subprocess
import pandas as pd
import re 
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt 
from subprocess import PIPE,run
import random
from statannot import add_stat_annotation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

from vcf2bed import vcf2bed 
from data_utils import extractInfo, listtoDict, filterDELGenes
from file_utils import importRNA, cellLink, getPath, findDELGenes, addExon, addIntron, importRef
from bam_proc import baminfo, annotateDELReads, importFormalBam, addFormalBam, combineCells

def MaxBetween(nestedList, regRange,chr):
	'''
	input: a nested list - [['chrom','start','end','type'],[1,3,10,'DEL'],['X',10,300,'INS']]
			a range for fuzzy match
	output: a list containing combined events  
   '''
	intArray = np.array(nestedList)[:,0:2].astype(int)
	strList = np.array(nestedList)[:,2:6].tolist()
	orig = np.array(nestedList)[:,0:2].astype(int)
	for interval in intArray : 
	#Barray: whether the other  intervals are wihtin the range of the var interval. return a boolean array whose length is 4 (the 4 ends of the ranges )
		BList = (intArray[:,0] >= interval[0] - regRange) & (intArray[:,0] <= interval[0] + regRange) & \
				(intArray[:,1] >= interval[1] - regRange) & (intArray[:,1] <= interval[1] + regRange) & \
				(intArray[:,0] < interval[1]) & (intArray[:,1] > interval[0])
		# extract the 'true' unit from nestedList  
		if intArray[BList].size!=0:
			intArray[BList, 0] = intArray[BList, 0].min()
			intArray[BList, 1] = intArray[BList, 1].max()
	# the return info :chrom start end com_link orig_start,orig_end,orig_link cell length
	intL = intArray[:,0:2].tolist()
	origL = orig[:,0:2].tolist()
	newList = []
	for i in range(0,len(intL)):
		c = []
		com_link = ':'.join([str(chr),str(intL[i][0]),str(intL[i][1])])
		c = [chr]+intL[i]+[com_link]+origL[i]+strList[i]
		newList.append(c)
	return newList

def extractInterval(dict1,regRange):
	chromList = [str(x) for x in range(1,23)] + ['X', 'Y']
	svtype = ['DEL', 'DUP', 'INS', 'INV', 'INVDUP']
	# to store all the cell intervals under one chrome and one type
	appList = []
	for chrom in chromList:
		for typ in svtype:
		# save the information
			cellInt = []	
			for cell in dict1.keys():
				for line in dict1[cell]:
					# for one type and one 
					if line[0] == chrom and line[3] == typ:
						cellInt.append([line[1],line[2],line[5],str(cell),line[4],line[7]])
			if cellInt != []:
				cellInt = np.array(cellInt)
				aList = MaxBetween(cellInt, regRange,chrom)
				# create a new list to stash type
				for i in aList:
					i.append(typ)
					appList.append(i)
	return  pd.DataFrame (appList,columns=['#chrom','comStart','comEnd','comLink','origStart','origEnd','origLink','cell','length','alt','type'])

def infoSum(df2):
	'''
	summarize the SV type info
	input: 
	output: df 
	'''
	newList = []
	for i in np.unique(df2.comLink):
		rowList = []
		comdf = df2[df2['comLink'] == i]
		rowChrLen = len(np.unique(comdf['#chrom']))
		rowTypLen = len(np.unique(comdf['type']))
		if (rowTypLen == 1) & (rowChrLen == 1) :
			times = len(np.unique(comdf.cell))
			chr = comdf['#chrom'][0:1].iloc[0]
			start = comdf.comStart[0:1].iloc[0]
			end = comdf.comEnd[0:1].iloc[0]
			comLen = end-start
			comLink = comdf.comLink[0:1].iloc[0]
			origList = []
			cellList = []
			lengthList = []
			for i in range(len(comdf)):
				origList.append(":".join([str(comdf.iloc[i].origStart),str(comdf.iloc[i].origEnd)]))
				cellList.append(comdf.iloc[i].cell)
				lengthList.append(comdf.iloc[i].length)
			origLink = ";".join(origList)
			length = ";".join(lengthList)
			cell = ";".join(cellList)
			type = comdf.type[0:1].iloc[0]
		rowList = [chr,start,end,comLink,comLen,times,cell,length,type,origLink]
		newList.append(rowList)
	return  pd.DataFrame (newList,columns=['#chrom','comStart','comEnd','comLink','comLength','times','cell','length','type','origLink'])

def createDelRef(dict1,dict2,chrom):
	'''
	create reference sequence according to SV position
	dict1 = reference info dict
	dict2 = SV info dict
	chrom : the chromosome we focus on 
	'''
	refSeq = dict1[chrom].seq
	# SV info
	SVinfo = dict2[chrom] # a nested list
	# find the positions and note down the new positions
	newRef = str()
	newPosList = []
	for i in np.arange(len(SVinfo)):
		start = int(SVinfo[i][1])
		end = int(SVinfo[i][2])
		if start < end :
			breakPoint1 = (2 * i + 1) * 500
			breakPoint2 = ((2 * i + 1) * 500)+ 1
			newStart = (i * 1000) + 1
			newEnd = (i + 1) * 1000
			tmpList = [chrom,newStart,newEnd,breakPoint1,breakPoint2] + SVinfo[i]
			front_half = refSeq[(start-500) : start]
			behind_half = refSeq[(end-1) : (end+499)]
			Ref1 = front_half + behind_half
		newRef += Ref1
		newPosList.append(tmpList)
	return newRef,newPosList

def createDELFasta(dict1,dict2,fastaFile):
	'''
	to integrate the new ref into a fasta file
	input: like the function 'createDelRef'
	'''
	chromList = list(dict2.keys())
	recordList = []
	gtf = {}
	for chrom in chromList:
		singleRef,singleList = createDelRef(dict1,dict2,chrom)
		# create biopython record object 
		record = SeqRecord(singleRef, chrom)
		recordList.append(record)
		gtf[chrom] = singleList
	SeqIO.write(recordList, fastaFile, "fasta")
	return recordList,gtf

def refineDEL(df):
	'''
	refine DEL reads to  find well-supportedly new transcripts
	df : all dels read info
	'''
	#1. find the uniquely aligned-to-new-genome reads
	df1 = df.loc[df.offiCigar == 'None']
	#2. find the reads that are almost equally spanning the bp1
	# readStart<=bp1-50 & readEnd >= bp2+50
	df2 = df1.loc[(df1.readStart.astype(int) <= (df1.bp1.astype(int)-50)) & (df1.readEnd.astype(int) >= (df1.bp2.astype(int)+50))]
	return df2

if __name__ == '__main__':
	# Notice: when upload final version, the path must be cleaned, like
	# 'path/to/XXX.file'
	'''
	expdf = importRNA('/datd/wangjun/SV/RNA/matrix/7_matrix.txt')
	celldf = cellLink('/datb/wangjun/SV/CELL.txt')
	cellList,fileList,bedList = getPath('/datd/zhouwei/01Projects/03ecDNA/Nanopore/04MINIMAPontSV','/datd/wangjun/SV/bed')
	for i in range(len(fileList)):
		vcf2bed(fileList[i],bedList[i])
	dfDict = extractInfo(cellList,bedList)
	othdf = extractInterval(dfDict,500)
	othdfsum = infoSum(othdf)
	dels = othdf.loc[othdf.type == 'DEL']
	delsum = othdfsum.loc[othdfsum.type=='DEL']
	delsum.to_csv('/datd/wangjun/SV/gittest/del.bed',index=False,sep='\t')
	delGene = findDELGenes('/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gtf.gene.bed','/datd/wangjun/SV/gittest/del.bed')
	delGene.to_csv('/datd/wangjun/SV/gittest/delGene.bed',index=False,sep='\t')
	delIntron = addIntron('/share/home/wang_jun/datb/SV/Homo_sapiens.GRCh38.100.gtf.exon.bed','/datd/wangjun/SV/gittest/delGene.bed')
	delExon = addExon('/share/home/wang_jun/datb/SV/Homo_sapiens.GRCh38.100.gtf.exon.bed','/datd/wangjun/SV/gittest/delGene.bed')
	delGene = pd.concat([delExon,delIntron])
	delGene.to_csv('/datd/wangjun/SV/gittest/delEI.bed',index=False,sep='\t')
	delDict = filterDELGenes(delGene,expdf)
	record_dict = importRef('/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
	delFasta,delGtf = createDELFasta(record_dict,delDict,"/datd/wangjun/SV/gittest/del.fasta")
	'''
	subprocess.run(["bash","bwaAlign.sh","/datd/wangjun/SV/gittest","/datd/wangjun/SV/gittest/del.fasta","/data/zhouwei/02production/20200914_1800"])
	'''
	alldels = combineCells(celldf,'/share/home/wang_jun/datb/SV/del/','/datd/zhouwei/02production/20210204_0935/',delGtf)
	delrefine = refineDEL(alldels)
	'''