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

def importRNA(rnaPath):
	'''
	import RNA data for the cells
	link the RNA cell name with DNA cell name
	'''
	rnaExprs = pd.read_csv(rnaPath,sep='\t')
	new_sub = []
	for sub in rnaExprs.gene_id :
		new_sub.append(sub.split('_')[1])
	rnaExprs.gene_id = new_sub
	return rnaExprs

def cellLink(cellPath):
	'''
	import txt file containing cell info
	'''
	cellInfo = pd.read_csv(cellPath,sep='\t')
	return cellInfo

def getPath(inpath,outpath):
	cellList = []
	fileList = []
	bedList = []
	for i in os.listdir(inpath) :
		cell = i.split('.')[0]
		cellList.append(cell)
		fileList.append(inpath+'/'+ i)
		bedList.append(outpath+'/'+ cell +'.bed' )
	return cellList,fileList,bedList

def findDELGenes(refbed,sumfile):
	'''
	using bedtools to find the SV-related genes
	command line :
	bedtools intersect -a /share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gtf.gene.bed -b merged.bed -wb
	refbed: downloaded reference bedfile
	sumfile: the summary bedfile of SV 
	geneDF: the dataframe containing sv-related genes
	'''
	print("find genes ...")
	command = "bedtools intersect -a "
	command += refbed
	command += " -b "
	command += sumfile
	command += " -wb "
	#command += resultFile
	tempLine = run(command, stdout = PIPE, stderr = PIPE, text=True, shell=True)
	Line = tempLine.stdout
	newLine = str.rstrip(Line).split('\n')
	#header = ['chrom','start','end','type','geneLen','strand','gene','ensemblID','geneType','chrom1',"start1",'end1','comLink','svLen','times','cell','lenLink','svType','origLen']
	df = pd.DataFrame([sub.split("\t") for sub in newLine])
	dfkeep = df.loc[:,[0,10,11,6,5,8,12,13,14,17,15]]
	headers = ['#chrom','start','end','gene','strand','geneType','comLink','svLen','times','svType','cell']
	dfkeep.columns = headers
	return dfkeep

def addExon(refbed,sumfile):
	'''
	using an exon annotation file to add exon, start_codon, stop_codon information on the gene-annotated SVs
	the intron part = gene - {exon, start_codon, stop_codon}
	refbed: homo_Sapiens.gtf.exon.bed
	sumfile: the gene related bed file
	'''
	print("add exon info ...")
	command = "bedtools intersect -a "
	command += sumfile
	command += " -b "
	command += refbed
	command += " -wa "
	#command += resultFile
	tempLine = run(command, stdout = PIPE, stderr = PIPE, text=True, shell=True)
	Line = tempLine.stdout
	newLine = str.rstrip(Line).split('\n')
	df = pd.DataFrame([sub.split("\t") for sub in newLine])
	headers = ['#chrom','start','end','gene','strand','geneType','comLink','svLen','times','svType','cell']	
	df.columns = headers
	df['partType'] = 'exon'
	df = df.drop_duplicates()
	return df

def addIntron(refbed,sumfile):
	'''
	find extron : genedf extra exondf
	'''
	# 1. find the SV that has exon in it according to comLink . (how many exons in the SV/ partially or whole will be ignored)
	print("add intron info ...")
	command = "bedtools intersect -a "
	command += sumfile
	command += " -b "
	command += refbed
	command += " -v "
	#command += resultFile
	tempLine = run(command, stdout = PIPE, stderr = PIPE, text=True, shell=True)
	Line = tempLine.stdout
	newLine = str.rstrip(Line).split('\n')
	headers = ['#chrom','start','end','gene','strand','geneType','comLink','svLen','times','svType','cell']
	df = pd.DataFrame([sub.split("\t") for sub in newLine])
	df.columns = headers
	df['partType'] = 'intron'
	df = df.drop_duplicates()
	return df

def importRef(fastaPath):
	'''
	import the large fasta file into one chrom per file using biopython
	input: fastaPath - the reference fasta file
	return 
	'''
	record_dict = SeqIO.index(fastaPath, "fasta")
	return record_dict