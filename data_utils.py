#!/bin/usr/python
**********************************************************
* @File    : data_utils.py			         *				
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

def extractInfo(cellList,bedList):
	header = ['#chrom','chromStart','chromEnd','type','length','link','quality','alt']
	dfDict = {}
	for i in range(len(bedList)):
		dfList = []
		df = pd.read_csv(bedList[i],sep='\t')[header]
		for j in range(len(df)):
			dfList.append(list(df.loc[j]))
		dfDict[cellList[i]] = dfList
	return dfDict

def listtoDict(list1):
	'''
	to split the list into dict by chromes ['chrom','.........'] the first element should be chrom and the chrom will be the key. The remaining part is the values
	'''
	##1. create the keyList
	infoArr = np.array(list1)
	keyList =list(np.unique(infoArr[:,0]))
	##2.create the dict
	infoDict = {}
	for i in keyList:
		temArr = infoArr[infoArr[:,0]==i]
		temList = []
		for j in temArr:
			temList.append(list(j))
		infoDict[i] = temList
	return infoDict

def filterDELGenes(df2,expdf):
	'''
	filter out genes
	input: insertion df with annotated genes 
	'''
	typdf = df2.loc[(df2.times.astype(int) >=3) & ((df2.geneType == 'protein_coding') | (df2.geneType == 'lncRNA'))]
	# find genes that have expression
	readEXP = expdf.loc[expdf.gene_id.isin(np.unique(typdf.gene).tolist())]
	meanExprs = readEXP.mean(axis=1).tolist()
	readEXP = readEXP.assign(meanExprs=meanExprs)
	readEXP = readEXP.loc[readEXP['meanExprs'].astype(int) > 0]
	# filter out the matrix
	infodf = typdf.loc[typdf['gene'].isin(list(readEXP.gene_id))]
	infoList = []
	for i in infodf.iloc:
		infoList.append(list(i))
	infoDict = listtoDict(infoList)
	return infoDict
