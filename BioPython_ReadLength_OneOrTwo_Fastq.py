#!/usr/bin/python
## import base Python packages for OS, RegExs, .csv files, statistics, etc.
import os, sys, re, csv, statistics
import argparse, logging, warnings, glob

## import BioPython and related packages
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

## import matplotlib
from matplotlib import pyplot as plt

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
parser = argparse.ArgumentParser(description='Show read count and base pair count of one or two filename.fastq files', usage='BioPython_Parsing_FASTQ.py filepath/filename1.fastq filepath/filename2.fastq (optional)')

## Add attribute 'filename' (the .fastq input file) to object 'parser'
parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

## Add attribute '--outputType' to object 'parser'
parser.add_argument('--outputType', '-o', default='S', choices=['S', 'L', 'P'], help="--outputType S for simple statistics, --outputType L for list of all read lengths, and --outputType P for matplotlib plots")

args = parser.parse_args()

## Function to plot read lengths of one .fastq file
def plotSingleReadLengths(readLengths, fileStr):
    SMALL_SIZE = 22
    MEDIUM_SIZE = 26
    BIG_SIZE = 32
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(14,20))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    axes1.hist(readLengths, bins = 40, color='blue')
    axes1.set_title(fileStr, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read Counts')
    axes1.set(xlabel='Read Lengths')
    fig1.savefig('/scicomp/home-pure/ydn3/nextflow_2023_for_read_mapping/SARS-CoV-2_MiSeq_VPipe_processed/readLength' + fileStr + '.png')   
## Function to plot read lengths of two .fastq files
def plotDoubleReadLengths(readLengths1, readLengths2, fileStr1, fileStr2):
    SMALL_SIZE = 20
    MEDIUM_SIZE = 24
    BIG_SIZE = 30
    fig1, axes1 = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(14,20))
    axes1[0].xaxis.label.set_size(MEDIUM_SIZE)
    axes1[0].yaxis.label.set_size(MEDIUM_SIZE)
    axes1[0].tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1[0].tick_params(axis='y', labelsize=SMALL_SIZE)
    axes1[0].hist(readLengths1, bins = 40, color='blue')
    axes1[0].set_title(fileStr1, fontsize = BIG_SIZE)
    axes1[0].set(ylabel='Read Counts')
    axes1[0].set(xlabel='Read Lengths')
    #fig1.savefig('/scicomp/home-pure/ydn3/nextflow_2023_for_read_mapping/SARS-CoV-2_MiSeq_VPipe_processed/readLength' + fileStr1 + '.png')   

    axes1[1].xaxis.label.set_size(MEDIUM_SIZE)
    axes1[1].yaxis.label.set_size(MEDIUM_SIZE)
    axes1[1].tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1[1].tick_params(axis='y', labelsize=SMALL_SIZE)
    axes1[1].hist(readLengths2, bins = 40, color='red')
    axes1[1].set_title(fileStr2, fontsize = BIG_SIZE)
    axes1[1].set(ylabel='Read Counts')
    axes1[1].set(xlabel='Read Lengths')
    fig1.savefig('/scicomp/home-pure/ydn3/nextflow_2023_for_read_mapping/SARS-CoV-2_MiSeq_VPipe_processed/readLen_' + fileStr1 + '_' + fileStr2 + '.png')   


## Function to iterate through the .fastq input file and count NGS reads of non-zero length
def lengthOfFastqReads(fnames, choice):
    readCount = 0
    readLengths = []
    readLengths2 = []
    f = 0
    if(len(fnames) < 2):
        with open(fnames[f].name, 'r') as fastq:
            for header, sequence, quality in FastqGeneralIterator(fastq):
                if(len(sequence) > 0):
                    readCount = readCount + 1
                    readLengths.append(len(sequence))
        if(choice == 'S'):
            print("%s average read length is %i base pairs" % (getIsolateStr(fnames[f].name), statistics.mean(readLengths)))
        elif(choice == 'L'):
            print("%s" % getIsolateStr(fnames[f].name))
            [print(item) for item in readLengths]
        else:
            plotSingleReadLengths(readLengths, getIsolateStr(fnames[f].name))
    elif(len(fnames) == 2):
        with open(fnames[0].name, 'r') as fastq:
            for header, sequence, quality in FastqGeneralIterator(fastq):
                if(len(sequence) > 0):
                    readLengths.append(len(sequence))
        with open(fnames[1].name, 'r') as fastq2:
            for header, sequence, quality in FastqGeneralIterator(fastq2):
                if(len(sequence) > 0):
                    readLengths2.append(len(sequence))
        if(choice == 'S'):
            print("%s average read length is %i base pairs" % (getIsolateStr(fnames[0].name), statistics.mean(readLengths)))
            print("%s average read length is %i base pairs" % (getIsolateStr(fnames[1].name), statistics.mean(readLengths2)))
        elif(choice == 'L'):
            if(len(readLengths) == len(readLengths2)):
                print("%s\t%s" % (getIsolateStr(fnames[0].name), getIsolateStr(fnames[1].name)))
                ii = 0
                while(ii < len(readLengths)):
                    print("%i\t%i" % (readLengths[ii], readLengths2[ii]))
                    ii = ii + 1
            elif(len(readLengths) > len(readLengths2)):
                print("%s\t%s" % (getIsolateStr(fnames[0].name), getIsolateStr(fnames[1].name)))
                ii = 0
                jj = len(readLengths2) - 1
                while(ii < len(readLengths2)):
                    print("%i\t%i" % (readLengths[ii], readLengths2[ii]))
                    ii = ii + 1
                while(jj < len(readLengths)):
                    print("%i\tNaN" % (readLengths[jj]))
                    jj = jj + 1
            else:
                print("%s\t%s" % (getIsolateStr(fnames[0].name), getIsolateStr(fnames[1].name)))
                ii = 0
                jj = len(readLengths) - 1
                while(ii < len(readLengths)):
                    print("%i\t%i" % (readLengths[ii], readLengths2[ii]))
                    ii = ii + 1
                while(jj < len(readLengths2)):
                    print("NaN\t%i" % (readLengths2[jj]))
                    jj = jj + 1
        else:
            plotDoubleReadLengths(readLengths, readLengths2, getIsolateStr(fnames[0].name), getIsolateStr(fnames[1].name))
                
## Exit on error if input files exceeds 2
if(len(args.filename) > 2):
    sys.exit("BioPython_LengthMultiFastq_Reads.py accepts only one or two .fastq files as input.")

## Invoke function that counts NGS reads for each file in command-line arguments

lengthOfFastqReads(args.filename, args.outputType)




