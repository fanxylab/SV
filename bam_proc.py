#!/bin/usr/python

**********************************************************
* @File    : bam_proc.py					 *				
* @Author  : WANG Jun                                    *
* @Date    : 2021/03                                     *
* @E-mail  : oucwj@outlook.com                           *               
**********************************************************


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

def baminfo(bamfile):
	'''
	inport the bamfile and extract useful info
	'''
	bam_in = pysam.AlignmentFile(bamfile, 'rb')
	readInfo = []
	for read in bam_in.fetch():
		cigartuples =  read.cigartuples 
		chrom  = read.reference_name
		start, end =  read.reference_start, read.reference_end
		cigarstring = read.cigarstring
		query_name = read.query_name
		query_length  =read.query_length   
		query_alignment_length = read.query_alignment_length
		if (int(start)%1000 < 500) & (int(end)%1000 > 501) & (int(query_alignment_length)/int(query_length) >= 0.8):
	 		readInfo.append([chrom,start,end,cigarstring,query_alignment_length,query_name])
	readDict = listtoDict(readInfo)
	return readDict

def annotateDELReads(gtf,dict1):
	'''
	annotate the reads by gene names
	dict1: the dictionary containing realigned reads
	gtf: gtf info generated from the fasta creation part
	'''
	keyList = list(dict1.keys())
	annoList = []
	for key in keyList:
		# extract annotation regions
		for a in dict1[key]:
			for b in gtf[key] :
				if (int(a[1]) >= int(b[1])) & (int(a[2]) <= int(b[2])):
					annoList.append(a + list(np.array(b)[1:5]) + list(np.array(b)[6:]))
	return pd.DataFrame(annoList,columns=['chrom','readStart','readEnd','readCigar','mappedLen','readName','annoStart','annoEnd','bp1','bp2','genomeStart','genomeEnd','gene','strand','geneType','comLink','svLen','times','svType','cellLink','partType'])



def importFormalBam(form_bamfile):
	'''
	import the bam file which aligned to GRh38 reference genome using STAR
	'''
	bam_in = pysam.AlignmentFile(form_bamfile, 'rb')
	formDict = {}
	for read in bam_in.fetch(until_eof=True):
		chrom  = read.reference_name
		cigarstring = read.cigarstring
		query_name = read.query_name
		formDict[query_name] = [cigarstring,chrom]
	return formDict


def addFormalBam(df1,dict2,cell):
	'''
	df1: the annotated read Dict---delRead
	dict2: the bam file list aligned to official/formal ref genome --- 
	'''
	df1['official'] = 'Unaligned'
	df1.loc[(df1.readName.isin(list(dict2.keys()))),'official' ] = 'Aligned'
	# add official Cigar and chrom
	chrList = []
	cigarList = []
	for i in np.arange(len(df1)):
		if df1.iloc[i]['official'] == 'Aligned':
			cigarList.append(str(dict2[df1.iloc[i]['readName']][0]))
			chrList.append(str(dict2[df1.iloc[i]['readName']][1]))
		else:
			cigarList.append('unknown')
			chrList.append('unknown')
	df1['offiCigar'] = cigarList
	df1['offiChrom'] = chrList
	df1['RNAcell'] = cell
	return df1

def combineCells(celldf,delbampath,formbampath,gtf):
	'''
	combine all the info from cells in celldf
	'''
	df = pd.DataFrame()
	for i in celldf.RNA:
		bamfile = delbampath+i+'.sorted.bam'
		form_bamfile = formbampath+i+'/SS2/STAR/'+i+'__D1.Aligned.sortedByCoord.sort.markdup.bam'
		# import the official 
		formDict = importFormalBam(form_bamfile)
		# create the annoList
		readDict = baminfo(bamfile)
		readdf = annotateDELReads(gtf,readDict)
		annodf = addFormalBam(readdf,formDict,i)
		annodf.to_csv(i+'.csv')
		df = pd.concat([df,annodf])
		df.to_csv('deldf.csv')
	return df
