#!/usr/bin/python
## import base Python packages for OS, RegExs, .csv files, statistics, etc.
import os, sys, re, csv, statistics
import argparse, logging, warnings, glob

## import BioPython and related packages
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

## import numpy
import numpy as np

## import matplotlib
from matplotlib import pyplot as plt

## import pandas
import pandas as pd

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
    SMALL_SIZE = 28
    MEDIUM_SIZE = 32
    BIG_SIZE = 36
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00']
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(15,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [dfCols[0], dfCols[len(dfCols) -1]]
    bp = axes1.boxplot(readLengths, notch=False, sym='+', vert=True, patch_artist=True, boxprops = dict(linewidth = 5), medianprops=medianprops, whis=[15, 85], whiskerprops=whiskerprops, capprops=capprops)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    axes1.set_title(fileStr, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read Counts')
    axes1.set(xlabel='Read Lengths')
    fig1.savefig('/scicomp/home-pure/ydn3/nextflow_2023_for_read_mapping/SARS-CoV-2_MiSeq_VPipe_processed/boxLength' + fileStr + '.png')   

## Function to plot read lengths of two or more .fastq files
def plotDoubleReadLenths(readLengthDF):
    SMALL_SIZE = 28
    MEDIUM_SIZE = 32
    BIG_SIZE = 36
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00']
    fig1, axes1 = plt.subplots(figsize=(25,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    dfCols = list(readLengthDF)
    fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    readLength1 = readLengthDF[dfCols[0]]
    readLength2 = readLengthDF[dfCols[1]]
    filteredReadLength1 = readLength1[~np.isnan(readLength1)]
    filteredReadLength2 = readLength2[~np.isnan(readLength2)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [dfCols[0], dfCols[len(dfCols) -1]]
    bp = axes1.boxplot([filteredReadLength1, filteredReadLength2], labels=x_labels, notch=False, sym='+', vert=True, patch_artist=True, boxprops = dict(linewidth = 5), medianprops=medianprops, whis=[15, 85], whiskerprops=whiskerprops, capprops=capprops)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    axes1.set_title(fileTitle, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read Lengths (bp)')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/home-pure/ydn3/nextflow_2023_for_read_mapping/SARS-CoV-2_MiSeq_VPipe_processed/boxLen_' + fileTitle + '.png')   

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
## for exactly two files, invoke plotMultiReadLengths
    elif(len(fnames) == 2):
        with open(fnames[0].name, 'r') as fastq:
            for header, sequence, quality in FastqGeneralIterator(fastq):
                if(len(sequence) > 0):
                    readLengths.append(len(sequence))
        ## Convert readLengths list to pandas Series
        readLenDF = pd.Series(readLengths)
        with open(fnames[1].name, 'r') as fastq2:
            for header, sequence, quality in FastqGeneralIterator(fastq2):
                if(len(sequence) > 0):
                    readLengths2.append(len(sequence))
        ## Create column headers
        simpFileName1 = getIsolateStr(fnames[0].name)
        simpFileName2 = getIsolateStr(fnames[1].name)
        fnamesl1 = simpFileName1.split('.')
        if(re.search('R1', simpFileName1, re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R1'
        elif(re.search('R2', simpFileName1, re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R2'
        #print(fnamesl1[0])
        fnamesl2 = simpFileName2.split('.')
        if(re.search('R1', simpFileName2, re.IGNORECASE)):
            fnamesl2[0] = fnamesl2[0] + '_R1'
        elif(re.search('R2', simpFileName2, re.IGNORECASE)):
            fnamesl2[0] = fnamesl2[0] + '_R2'
        ## Convert readLengths2 list to pandas Series
        readLenDF2 = pd.Series(readLengths2)
        ## Create pandas DataFrame with truncated fnames as column headers
        readLengthDF = pd.DataFrame({fnamesl1[0] : readLenDF, fnamesl2[0] : readLenDF2})
        ## Coerce string 'Nan' to np.nan
        #readLengthDF[fnamesl2[0]] = readLengthDF[fnamesl2[0]].replace(r'Nan', np.nan, regex=True)
        if(choice == 'S'):
            print("%s average read length is %i base pairs" % (fnamesl1[0], readLengthDF[fnamesl1[0]].mean() ))
            print("%s average read length is %i base pairs" % (fnamesl2[0], readLengthDF[fnamesl2[0]].mean() ))
        elif(choice == 'L'):
            print("%s\t%s" % (fnamesl1[0], fnamesl2[0]))
            for index, row in readLengthDF.iterrows():
                print(row[fnamesl1[0]], "\t", row[fnamesl2[0]])
        else:
            plotDoubleReadLenths(readLengthDF)
    elif((len(fnames) == 3) or (len(fnames) == 4)):
        ii = 0
        readLenDF = pd.Series(dtype=int)
        readLenDF2 = pd.Series(dtype=int)
        readLengthDF = pd.DataFrame()
        simpFileName = []
        while(ii < len(fnames)):
            readLengths = []
            with open(fnames[ii].name, 'r') as fastq:
                for header, sequence, quality in FastqGeneralIterator(fastq):
                    if(len(sequence) > 0):
                        readLengths.append(len(sequence))
            ## Convert readLengths list to pandas Series
            readLenDF = pd.Series(readLengths)
            simpFileName.append(getIsolateStr(fnames[ii].name))
            fnamesl1 = simpFileName[ii].split('.')
            if(re.search('R1', simpFileName[ii], re.IGNORECASE)):
                fnamesl1[0] = fnamesl1[0] + '_R1'
            elif(re.search('R2', simpFileName[ii], re.IGNORECASE)):
                fnamesl1[0] = fnamesl1[0] + '_R2'
            readLengthDF[fnamesl1[0]] = readLenDF
            ii = ii + 1
        if(choice == 'S'):
            dfCols = list(readLengthDF)
            for cols in dfCols:
                print("%s average read length is %i base pairs" % (cols, readLengthDF[cols].mean() ))

## Exit on error if input files exceeds 2
if(len(args.filename) > 4):
    sys.exit("BioPython_LengthMultiFastq_Reads.py accepts only one or two .fastq files as input.")

## Invoke function that counts NGS reads for each file in command-line arguments

lengthOfFastqReads(args.filename, args.outputType)




