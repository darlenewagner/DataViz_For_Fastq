#!/usr/bin/python
## import base Python packages for OS, RegExs, .csv files, statistics, etc.
import os, sys, re, csv, statistics, math
import argparse, logging, warnings, glob

## import BioPython and related packages
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

## import numpy
import numpy as np
import seaborn as sns

## import matplotlib
from matplotlib import pyplot as plt

## import pandas
import pandas as pd

## Function to validate type of input file: '.fastq'
def ext_check(expected_ext1, expected_ext2, expected_ext3, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
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
parser = argparse.ArgumentParser(description='Show read count and base pair count of one to six filename.fastq files (alternate: Insert length .tsv files)', usage='BioPython_Parsing_FASTQ.py filepath/filename1.fastq filepath/filename2.fastq (optional)')

## Add attribute 'filename' (the .fastq input file) to object 'parser'
parser.add_argument('filename', nargs='+', type=ext_check('.fastq', '.tsv', '.csv', argparse.FileType('r')))

## Add attribute '--outputType' to object 'parser'
parser.add_argument('--outputType', '-o', default='S', choices=['S', 'L', 'P'], help="--outputType S for simple statistics, --outputType L for list of all read lengths, and --outputType P for matplotlib plots")

args = parser.parse_args()

## Function to extract Fastq read length (other version extracts PHRED Quality)
def extractFastqData(fnames):
    simpFileName = []
    ii = 0
    readLengthDF = pd.DataFrame()
    while(ii < len(fnames)):
        readLengths = []
        with open(fnames[ii].name, 'r') as fastq:
        ## Invoke Bio.SeqIO.QualityIO
            for header, sequence, quality in FastqGeneralIterator(fastq):
                if(len(sequence) > 0):
                    readLengths.append(len(sequence))
        readLenDF = pd.Series(readLengths)
        simpFileName.append(getIsolateStr(fnames[ii].name))
        fnamesl1 = simpFileName[ii].split('.')
        if(re.search('R1', simpFileName[ii], re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R1'
        elif(re.search('R2', simpFileName[ii], re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R2'
        readLengthDF[fnamesl1[0]] = readLenDF
        ii = ii + 1
    return(readLengthDF)

#def extractInsertLengths()


## Function to plot read lengths of one .fastq file
def plotSingleReadLengths(readLengths, fileStr):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00']
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(15,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [fileStr]
    vp = axes1.violinplot(readLengths, showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    axes1.set_xticks(np.arange(1, len(x_labels) + 1), labels=x_labels)
    xy = [[l.vertices[:,0].mean(),l.vertices[0,1]] for l in vp['cmeans'].get_paths()]
    xy = np.array(xy)
    axes1.scatter(xy[:,0], xy[:,1],s=121, c="#A020F0", marker="D", zorder=3)
    for pc, color in zip(vp['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
    vp['cmeans'].set_visible(False)
    vp['cmedians'].set_colors('black')
    vp['cmedians'].set_linewidth(3)
    axes1.set_title(fileStr, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read Lengths (bp)')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLength_' + fileStr + '.png')   

## Function to plot read lengths of two or more .fastq files
def plotDoubleReadLenths(readLengthDF):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
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
    vp = axes1.violinplot([filteredReadLength1[filteredReadLength1 < 400], filteredReadLength2[filteredReadLength2 < 400]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    axes1.set_xticks(np.arange(1, len(x_labels) + 1), labels=x_labels)
    #axes1.set_ylim(ymax=500)
    xy = [[l.vertices[:,0].mean(),l.vertices[0,1]] for l in vp['cmeans'].get_paths()]
    xy = np.array(xy)
    axes1.scatter(xy[:,0], xy[:,1],s=121, c="#A020F0", marker="D", zorder=3)
    for pc, color in zip(vp['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
    vp['cmeans'].set_visible(False)
    vp['cmedians'].set_colors('black')
    vp['cmedians'].set_linewidth(3)
    axes1.set_title(fileTitle, fontsize = BIG_SIZE)
    axes1.set(ylabel='Lengths (bp)')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLen_' + fileTitle + '.png')


def plotMultiReadLengths(readLengthDF):
    SMALL_SIZE = 30
    MEDIUM_SIZE = 36
    BIG_SIZE = 38
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00', '#00FFFF', '#FF00FF']
    fig1, axes1 = plt.subplots(figsize=(27,19))
    fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.22)
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE, labelrotation=25)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    axes1.margins(0.1)
    dfCols = list(readLengthDF)
    fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    filteredReadLength = {}
    readLength = []
    for col in dfCols:
        readLength = readLengthDF[col]
        filteredReadLength[col] = readLength[~np.isnan(readLength)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='#A020F0', markersize=16)
    x_labels = list(readLengthDF)
    #readLengthMatrix = readLengthDF.to_numpy()
    if(len(x_labels) == 3):
        vp = axes1.violinplot(filteredReadLength[x_labels[0]], filteredReadLength[x_labels[1]], filteredReadLength[x_labels[2]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 4):
        vp = axes1.violinplot([filteredReadLength[x_labels[0]], filteredReadLength[x_labels[1]], filteredReadLength[x_labels[2]], filteredReadLength[x_labels[3]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 5):
        vp = axes1.violinplot([filteredReadLength[x_labels[0]], filteredReadLength[x_labels[1]], filteredReadLength[x_labels[2]], filteredReadLength[x_labels[3]], filteredReadLength[x_labels[4]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 6):
        vp = axes1.violinplot([filteredReadLength[x_labels[0]], filteredReadLength[x_labels[1]], filteredReadLength[x_labels[2]], filteredReadLength[x_labels[3]], filteredReadLength[x_labels[4]], filteredReadLength[x_labels[5]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    axes1.set_xticks(np.arange(1, len(x_labels) + 1), labels=x_labels)
    xy = [[l.vertices[:,0].mean(),l.vertices[0,1]] for l in vp['cmeans'].get_paths()]
    xy = np.array(xy)
    axes1.scatter(xy[:,0], xy[:,1],s=121, c="#A020F0", marker="D", zorder=3)
    for pc, color in zip(vp['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
    vp['cmeans'].set_visible(False)
    vp['cmedians'].set_colors('black')
    vp['cmedians'].set_linewidth(3)
    axes1.set_title(fileTitle, fontsize = BIG_SIZE)
    axes1.set(ylabel='Lengths (bp)')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/multiViolinLen_' + fileTitle + '.png')
    

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
## for exactly two files, invoke plotDoubleReadLengths
    elif(len(fnames) == 2):
        readLengthDF = pd.DataFrame()
        ## Call function 'extractFastqData' to read multiple .fastq files
        readLengthDF = extractFastqData(fnames)
        
        if(choice == 'S'):
            print("%s average read length is %i base pairs" % ( getIsolateStr(fnames[0].name), readLengthDF[fnames[0]].mean() ))
            print("%s average read length is %i base pairs" % ( getIsolateStr(fnames[1].name), readLengthDF[fnames[1]].mean() ))
        elif(choice == 'L'):
            #print("%s\t%s" % (fnamesl1[0], fnamesl2[0]))
            #for index, row in readLengthDF.iterrows():
            #    print(row[fnamesl1[0]], "\t", row[fnamesl2[0]])
            dfCols = list(readLengthDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in readLengthDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
        ## Call funtion to plot two boxplots
            plotDoubleReadLenths(readLengthDF)
## for three or four files, invoke plotDoubleReadLengths
    elif((len(fnames) > 2) and (len(fnames) < 7)):
        readLenDF = pd.Series(dtype=int)
        readLenDF2 = pd.Series(dtype=int)
        readLengthDF = pd.DataFrame()
        ## Call function 'extractFastqData' to read multiple .fastq files
        readLengthDF = extractFastqData(fnames)
        if(choice == 'S'):
            dfCols = list(readLengthDF)
            for cols in dfCols:
                print("%s average read length is %i base pairs" % (cols, readLengthDF[cols].mean() ))
        elif(choice == 'L'):
            dfCols = list(readLengthDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in readLengthDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
        ## Call function to plot more than two violin plots
            plotMultiReadLengths(readLengthDF)

def getInsertLengths(fnames, choice):
    insertCount = 0
    insertLengths = []
    insertLengths2 = []
    insertLengthDF = pd.DataFrame()
    f = 0
    if(len(fnames) < 2):
        with open(fnames[0].name, 'r') as tsv:
            for line in csv.reader(tsv, delimiter="\t"):
                if(insertCount > 0):
                    insertLengths.append(int(line[1]))
                insertCount = insertCount + 1

        #fnames[0].close()
        if(choice == 'S'):
            print("%s average insert length is %i base pairs" % (getIsolateStr(fnames[f].name), statistics.mean(insertLengths)))
        elif(choice == 'L'):
            print("%s" % getIsolateStr(fnames[f].name))
            [print(item) for item in insertLengths]
        else:
            plotSingleReadLengths(insertLengths, getIsolateStr(fnames[f].name))
## for exactly two files, invoke plotDoubleReadLengths
    elif(len(fnames) == 2):
        with open(fnames[0].name, 'r') as tsv1:
            for line in csv.reader(tsv1, delimiter="\t"):
                if(insertCount > 0):
                    insertLengths.append(int(line[1]))
                insertCount = insertCount + 1
        insertCount = 0
        with open(fnames[1].name, 'r') as tsv2:
            for line in csv.reader(tsv2, delimiter="\t"):
                if(insertCount > 0):
                    insertLengths2.append(int(line[1]))
                insertCount = insertCount + 1
        ## Call function 'extractFastqData' to read multiple .fastq files
        #insertLengthDF = extractFastqData(fnames)
        insertLengthTemp = pd.DataFrame(insertLengths, columns=[getIsolateStr(fnames[0].name)])
        insertLengthTemp2 = pd.DataFrame(insertLengths2, columns=[getIsolateStr(fnames[1].name)])
        insertLengthDF = pd.concat([insertLengthTemp, insertLengthTemp2], axis=1)
        
        if(choice == 'S'):
            #print(insertLengthDF.tail())
            for cols in list(insertLengthDF):
                print("%s average insert length is %i base pairs" % ( cols, insertLengthDF[cols].mean() ))
                #print("%s average insert length is %i base pairs" % ( cols, insertLengthDF[cols].mean() ))
        elif(choice == 'L'):
            dfCols = list(insertLengthDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in insertLengthDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
        ## Call funtion to plot two violin plots
            plotDoubleReadLenths(insertLengthDF)
            
    elif((len(fnames) > 2) and (len(fnames) < 7)):
        tempInsertSeries = pd.Series()
        fileIdx = 0
        while(fileIdx < len(fnames)):
            with open(fnames[fileIdx].name, 'r') as tsv0:
                for line in csv.reader(tsv0, delimiter="\t"):
                    if(insertCount > 0):
                        insertLengths.append(int(line[1]))
                    insertCount = insertCount + 1
                insertCount = 0
            insertLength = [i for i in insertLengths if i < 400]
            if(fileIdx == 0):
                insertLengthDF[getIsolateStr(fnames[fileIdx].name)] = insertLength
            else:
                insertLengthTemp = pd.DataFrame(insertLength, columns=[getIsolateStr(fnames[fileIdx].name)])
                insertLengthDF = pd.concat([insertLengthDF, insertLengthTemp], axis=1)
            fileIdx = fileIdx + 1
        if(choice == 'S'):
            for cols in list(insertLengthDF):
                print("%s average insert length is %i base pairs" % ( cols, insertLengthDF[cols].mean() ))
        elif(choice == 'L'):
            dfCols = list(insertLengthDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in insertLengthDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
        ## Call funtion to plot two violin plots
            plotMultiReadLengths(insertLengthDF)
        
## Exit on error if input files exceeds 6
if(len(args.filename) > 6):
    sys.exit("BioPython_LengthMultiFastq_Reads.py accepts no more than six .fastq files as input.")

fileStr = getIsolateStr(args.filename[0].name)

if(fileStr.endswith('.fastq')):
    ## Invoke function that plots NGS reads for each file in command-line arguments
    lengthOfFastqReads(args.filename, args.outputType)
elif(fileStr.endswith('.tsv') or fileStr.endswith('.csv')):
    ## Invoke function that plots NGS read insert lengths for each file in command-line arguments
    getInsertLengths(args.filename, args.outputType)


