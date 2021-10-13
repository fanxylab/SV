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

def getInfo(line):
	info = line.split('\t')[0:-3]
	name = info[0]
	start_int = int(info[1])
	alt_str = str(info[4])
	#end_int = int(info[1])+1
	return name, start_int,alt_str

def getSupInfo(line):	
	# supInfo 
	supInfo = line.split('\t')[-3]
	INFO = [{pair_lst[0]: pair_lst[1] if len(pair_lst)> 1 else ""} for pair_lst in [pair.split("=") for pair in supInfo.split(";")]]
	svtype_str = str(INFO[8]['SVTYPE'])
	strands_str = str(INFO[11]['STRANDS'])
	svlen_int = int(INFO[10]['SVLEN'])
	end_int = int(INFO[3]['END'])
	for key,value in INFO[0].items():
		quality_str = str(key)
	return svtype_str,strands_str,svlen_int,end_int,quality_str


def addLink(line):
	info = line.split('\t')[0:-3]
	link1 = [info[0],info[1]]
	linkl = ':'.join(link1)
	svtype_str,strands_str,svlen_int,end_int,quality_str = getSupInfo(line)
	if svtype_str == 'BND':
		link2 = re.search('(?:\[|\])(.*)(?:\[|\])',info[4]).group(1)
		link12 = [linkl,link2]
		linkEl = ';'.join(link12) 
	else: 
		linkEl = linkl
	return(linkEl)

def filter(df): 
	#df1 = df[(df['link'].str.lower()).str.islower()==False]
	df1 = []
	for i in range(len(df)):
		link = df.iloc[i].link
		sum = 0
		for j in re.split(';',link):
			if j.split(":")[0].isdigit() or j.split(":")[0] == 'X' or j.split(":")[0] == 'Y':
				sum = sum + 1
		if sum == len(re.split(';',link)):
			df1.append(df.iloc[i])
	df1 = pd.DataFrame(df1)
	# remove imprecise rows
	df1 = df1.loc[df1['quality'] == 'PRECISE']
	# split the df into 2 sub df
	# BND and INVDUP
	dfBID = df1[(df1.type == 'BND') | (df1.type == 'INVDUP')]
	dfOther = df1[(df1.type != 'BND') & (df1.type != 'INVDUP')]
	# filter the dfOther 
	dfOther = dfOther[((dfOther.type == 'DEL') & (abs(dfOther.length) <= 1e5) & (abs(dfOther.length) >= 100)) | (dfOther.type != 'DEL')]
	dfOther = dfOther[((dfOther.type == 'DUP') & (abs(dfOther.length) <= 1e4) & (abs(dfOther.length) >= 100)) | (dfOther.type != 'DUP')]
	dfOther = dfOther[((dfOther.type == 'INV') & (abs(dfOther.length) <= 1e4) & (abs(dfOther.length) >= 100)) | (dfOther.type != 'INV')]
	dfOther = dfOther[((dfOther.type == 'INS') & (abs(dfOther.length) >= 100)) | (dfOther.type != 'INS')]
	dflist = [dfOther,dfBID]
	df2 = pd.concat(dflist)
	return df2

def vcf2bed(vcf_path,bed_path):
	vcf_f = open(vcf_path,'r')
	line = vcf_f.readline().strip()
	name_list = []
	start_list = []
	end_list = []
	svtype_list = []
	strands_list = []
	link_list = []
	svlen_list = []
	quality_list = []
	alt_list = []
	while line:
		# Some line does not contain the data
		if (len(line.split('\t'))>5) & (line.split('\t')[0][0]!='#'):
			name, start_int, alt_str = getInfo(line)
			svtype_str,strands_str,svlen_int,end_int,quality_str = getSupInfo(line)
			link_str = addLink(line)
			if svtype_str == 'BND':
				name_list.append(name)
				nameAnother = re.search('(?:\[|\])(.*)(?:\[|\])',line.split('\t')[4]).group(1).split(":")[0]
				name_list.append(nameAnother)
				start_list.append(start_int)
				end_int1 = int(start_int+1)
				start_intAnother = int(re.search('(?:\[|\])(.*)(?:\[|\])',line.split('\t')[4]).group(1).split(":")[1])
				start_list.append(start_intAnother)
				end_int2 = int(start_intAnother+1)
				end_list.append(end_int1)
				end_list.append(end_int2)
				svtype_list.append(svtype_str)
				svtype_list.append(svtype_str)
				strands_list.append(strands_str)
				strands_list.append(strands_str)
				link_list.append(link_str)
				link_list.append(link_str)
				svlen_list.append(svlen_int)
				svlen_list.append(svlen_int)
				quality_list.append(quality_str)
				quality_list.append(quality_str)
				alt_list.append(alt_str)
				alt_list.append(alt_str)
			else:
				name_list.append(name)
				start_list.append(start_int)
				end_list.append(end_int)
				svtype_list.append(svtype_str)
				strands_list.append(strands_str)
				link_list.append(link_str)
				svlen_list.append(svlen_int)
				quality_list.append(quality_str)
				alt_list.append(alt_str)
		line = vcf_f.readline().strip()
	# save bed file
	bed_df = pd.DataFrame({'#chrom':name_list,'chromStart':start_list,'chromEnd':end_list,'length':svlen_list,'type':svtype_list,'strand':strands_list,'link':link_list,'quality':quality_list,'alt':alt_list})
	# add INS length to chromEnd
	bed_df.loc[(bed_df['type'] == 'INS'),'chromEnd'] =  bed_df['chromEnd'] + bed_df['length']
	# filter SV length
	bed_df = filter(bed_df)
	bed_df.to_csv(bed_path, index=False, sep='\t')
	return bed_df