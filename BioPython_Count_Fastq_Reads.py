#!/usr/bin/python
## import base Python packages for OS, RegExs, .csv files, statistics, etc.
import os, sys, re, csv, statistics
import argparse, logging, warnings

## import BioPython and related packages
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

## Function to validate type of input file: '.fastq'
def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.lower().endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

## Function to return filename without directory path
def getIsolateStr(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	fileString = splitStr[fileNameIdx]
	return fileString

## Initialize argparse object as 'parser'
parser = argparse.ArgumentParser(description='Show read count and base pair count of filename.fastq', usage='BioPython_Parsing_FASTQ.py filepath/filename.fastq')

## Add attribute 'filename' (the .fastq input file) to object 'parser'
parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

args = parser.parse_args()

iter = 0

## Function to iterate through the .fastq input file and count NGS reads of non-zero length
def countFastqReads(fname):
    readCount = 0
    with open(fname, 'r') as fastq:
        for header, sequence, quality in FastqGeneralIterator(fastq):
            if(len(sequence) > 0):
                readCount = readCount + 1
    print("%s has %i reads" % (getIsolateStr(fname), readCount))

## Invoke function that counts NGS reads
countFastqReads(args.filename[0].name)




        
    
