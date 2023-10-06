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

## Set up logger to show runtime progress to the user
logger = logging.getLogger("BioPython_PHREDViolinFastq_Reads.py")
logger.setLevel(logging.INFO)

## Usage statement
usage_text = ''' Examples:

BioPython_PHREDViolinFastq_Reads.py filepath/filename1.fastq filepath/filename2.fastq (optional) filepath/filename3.fastq (optional) --outputType S (for text summary)

BioPython_PHREDViolinFastq_Reads.py filepath/filename1.fastq filepath/filename2.fastq (optional) filepath/filename3.fastq (optional) --outputType P (for plot) --titleString "Raw Reads PHRED Scores" --showQ30 Y
 
'''

## Initialize argparse object as 'parser'
parser = argparse.ArgumentParser(description='Show read count and base pair count of one or two filename.fastq files', usage=usage_text)

## Add attribute 'filename' (the .fastq input file) to object 'parser'
parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

## Add attribute '--outputType' to object 'parser'
parser.add_argument('--outputType', '-o', default='S', choices=['S', 'L', 'P'], help="--outputType S for simple statistics, --outputType L for list of all read lengths, and --outputType P for matplotlib plots")

parser.add_argument('--titleString', '-t', default='Raw Reads PHRED', help="--titleString 'title for plot and filename'")

parser.add_argument('--showQ30', '-q', default='N', choices=['Y', 'N'], help="--showQ30 [Y/N] to print percentage where reads have PHRED scores above 30")

args = parser.parse_args()

## configuring the stream handler for logging
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)


## Function to extract Fastq read length (other version extracts PHRED Quality)
def extractData(fnames):
    simpFileName = []
    ii = 0
    readPHREDDF = pd.DataFrame()
    while(ii < len(fnames)):
        readPHRED = []
        ## Invoke Bio.SeqIO
        with open(fnames[ii].name, 'r') as fastq:
            for record in SeqIO.parse(fastq, "fastq"):
                if(len(record.seq) > 0):
                    readPHRED.append(statistics.mean(record.letter_annotations["phred_quality"]))
        readPhrDF = pd.Series(readPHRED)
        simpFileName.append(getIsolateStr(fnames[ii].name))
        fnamesl1 = simpFileName[ii].split('.')
        if(re.search('R1', simpFileName[ii], re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R1'
        elif(re.search('R2', simpFileName[ii], re.IGNORECASE)):
            fnamesl1[0] = fnamesl1[0] + '_R2'
        readPHREDDF[fnamesl1[0]] = readPhrDF
        ii = ii + 1
    return(readPHREDDF)


## Function to plot read lengths of one .fastq file
def plotSingleReadFile(readPHREDs, fileStr, title, choice, rQ30, rLen, fileMean, fileMedian):
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
    vp = axes1.violinplot(readPHREDs, showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    annotStr = "Median = %0.2f\nMean = %0.2f" % (fileMedian[0], fileMean[0])
    axes1.text(0.5, 25, annotStr, fontsize=26)
    if(choice == 'Y'):
        moreAnnotStr = "%0.2f" % (100*rQ30/rLen)
        axes1.text(0.5, 24, "Q30 = " + moreAnnotStr + "%", fontsize=26)
    axes1.set_title(title, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read PHRED Scores')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violinSinglePHRED_' + fileStr + '.png')   

## Function to plot read lengths of two or more .fastq files
def plotDoubleReadFiles(readPHREDDF, choice, title, rQ30, rLen, fileMean, fileMedian):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00']
    fig1, axes1 = plt.subplots(figsize=(25,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    dfCols = list(readPHREDDF)
    fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    readPhred1 = readPHREDDF[dfCols[0]]
    readPhred2 = readPHREDDF[dfCols[1]]
    filteredReadPhred1 = readPhred1[~np.isnan(readPhred1)]
    filteredReadPhred2 = readPhred2[~np.isnan(readPhred2)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [dfCols[0], dfCols[len(dfCols) -1]]
    vp = axes1.violinplot([filteredReadPhred1, filteredReadPhred2], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    axes1.set_title(title, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read PHRED Scores')
    axes1.set(xlabel='Fastq Files')
    annotStrA = "Median = %0.2f\nMean = %0.2f" % (fileMedian[0], fileMean[0])
    axes1.text(0.5, 25, annotStrA, fontsize=26)
    annotStrB = "Median = %0.2f\nMean = %0.2f" % (fileMedian[1], fileMean[1])
    axes1.text(1.5, 25, annotStrB, fontsize=26)
    if(choice == 'Y'):
        annotStr1 = "%0.2f" % (100*rQ30[0]/rLen[0])
        axes1.text(0.5, 24, "Q30 = " + annotStr1 + "%", fontsize=24)
        annotStr2 = "%0.2f" % (100*rQ30[1]/rLen[1])
        axes1.text(1.5, 24, "Q30 = " + annotStr2 + "%", fontsize=24)
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violin2filePHRED_' + fileTitle + '.png')

def plotTripleReadFiles(readPHREDDF, choice, title, rQ30, rLen, fileMean, fileMedian):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00']
    fig1, axes1 = plt.subplots(figsize=(25,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    dfCols = list(readPHREDDF)
    fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    readPhred1 = readPHREDDF[dfCols[0]]
    readPhred2 = readPHREDDF[dfCols[1]]
    readPhred3 = readPHREDDF[dfCols[2]]
    filteredReadPhred1 = readPhred1[~np.isnan(readPhred1)]
    filteredReadPhred2 = readPhred2[~np.isnan(readPhred2)]
    filteredReadPhred3 = readPhred3[~np.isnan(readPhred3)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [dfCols[0], dfCols[1], dfCols[len(dfCols) -1]]
    vp = axes1.violinplot([filteredReadPhred1, filteredReadPhred2, filteredReadPhred3], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    axes1.set_title(title, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read PHRED Scores')
    axes1.set(xlabel='Fastq Files')
    annotStrA = "Median = %0.2f\nMean = %0.2f" % (fileMedian[0], fileMean[0])
    axes1.text(0.5, 25, annotStrA, fontsize=26)
    annotStrB = "Median = %0.2f\nMean = %0.2f" % (fileMedian[1], fileMean[1])
    axes1.text(1.5, 25, annotStrB, fontsize=26)
    annotStrC = "Median = %0.2f\nMean = %0.2f" % (fileMedian[2], fileMean[2])
    axes1.text(2.5, 25, annotStrC, fontsize=26)
    if(choice == 'Y'):
        annotStr1 = "%0.2f" % (100*rQ30[0]/rLen[0])
        axes1.text(0.5, 24, "Q30 = " + annotStr1 + "%", fontsize=24)
        annotStr2 = "%0.2f" % (100*rQ30[1]/rLen[1])
        axes1.text(1.5, 24, "Q30 = " + annotStr2 + "%", fontsize=24)
        annotStr3 = "%0.2f" % (100*rQ30[2]/rLen[2])
        axes1.text(2.5, 24, "Q30 = " + annotStr3 + "%", fontsize=24)
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violin3filePHRED_' + fileTitle + '.png')


def plotMultiReadFiles(readPHREDDF, choice, title, rQ30, rLen, fileMean, fileMedian):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    textPositions = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00', '#00FFFF', '#FF00FF']
    fig1, axes1 = plt.subplots(figsize=(27,19))
    fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.22)
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE, labelrotation=25)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    axes1.margins(0.1)
    dfCols = list(readPHREDDF)
    fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    filteredReadPhred = {}
    readPhred = []
    for col in dfCols:
        readPhred = readPHREDDF[col]
        filteredReadPhred[col] = readPhred[~np.isnan(readPhred)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='#A020F0', markersize=16)
    x_labels = list(readPHREDDF)
    if(len(x_labels) == 3):
        vp = axes1.violinplot([filteredReadPhred[x_labels[0]], filteredReadPhred[x_labels[1]], filteredReadPhred[x_labels[2]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 4):
        vp = axes1.violinplot([filteredReadPhred[x_labels[0]], filteredReadPhred[x_labels[1]], filteredReadPhred[x_labels[2]], filteredReadPhred[x_labels[3]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 5):
        vp = axes1.violinplot([filteredReadPhred[x_labels[0]], filteredReadPhred[x_labels[1]], filteredReadPhred[x_labels[2]], filteredReadPhred[x_labels[3]], filteredReadPhred[x_labels[4]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    elif(len(x_labels) == 6):
        vp = axes1.violinplot([filteredReadPhred[x_labels[0]], filteredReadPhred[x_labels[1]], filteredReadPhred[x_labels[2]], filteredReadPhred[x_labels[3]], filteredReadPhred[x_labels[4]], filteredReadPhred[x_labels[5]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    p = 0
    while(p < len(x_labels)):
        annotStrA = "Median = %0.2f\nMean = %0.2f" % (fileMedian[p], fileMean[p])
        axes1.text(textPositions[p], 25, annotStrA, fontsize=22)
        p = p + 1
    p = 0
    if(choice == 'Y'):
        while(p < len(x_labels)):
            annotStr1 = "%0.2f" % (100*rQ30[p]/rLen[p])
            axes1.text(textPositions[p], 24.2, "Q30 = " + annotStr1 + "%", fontsize=22)
            p = p + 1
            
    axes1.set_title(title, fontsize = BIG_SIZE)
    axes1.set(ylabel='Read PHRED Scores')
    axes1.set(xlabel='Fastq Files')
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/ViolinMultiPHRED_' + fileTitle + '.png')


## Function to iterate through the .fastq input file and count NGS reads of non-zero length
def PHREDOfFastqReads(fnames, title, choice1, choice2, logger):
    readCount = 0
    fileMean = []
    fileMedian = []
    forwardAvg = []
    reverseAvg = []
    f = 0
## for exactly one file, invoke plotSingleReadFile
    if(len(fnames) < 2):
        #j = 0
        r1Q30 = 0
        r1Len = 1
        logger.info("Parsing one fastq file.")
        with open(fnames[f].name, 'r') as fastq:
            if(choice2 == 'Y'):
                logger.info("Starting calculations of Q30.")
            for record in SeqIO.parse(fastq, "fastq"):
                if(len(record.seq) > 0):
                    readCount = readCount + 1
                    forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
                    if(choice2 == 'Y'):
                        
                        r1Len = r1Len + len(record.seq)
                        j = 0
                        while( j < len(record.seq)):
                            if(record.letter_annotations["phred_quality"][j] >= 30):
                                r1Q30 = r1Q30 + 1
                                #print(r1Q30)
                            j = j + 1
            if(choice2 == 'Y'):
                logger.info("Completed calculations of Q30.")
        fileMean.append(statistics.mean(forwardAvg))
        fileMedian.append(statistics.median(forwardAvg))
        if(choice1 == 'S'):
            print("%s total average PHRED Score is %0.2f" % (getIsolateStr(fnames[f].name), statistics.mean(forwardAvg)))
            if(choice2 == 'Y'):
                print("%s PHRED Scores above Q30 is %0.2f" % (getIsolateStr(fnames[f].name), (100*r1Q30/r1Len)))
        elif(choice1 == 'L'):
            print("%s" % getIsolateStr(fnames[f].name))
            [print(item) for item in forwardAvg]
        else:
            logger.info("Plotting single series from one fastq file.")
            plotSingleReadFile(forwardAvg, getIsolateStr(fnames[f].name), title, choice2, r1Q30, r1Len, fileMean, fileMedian)
## for exactly two files, invoke plotDoubleReadFiles
    elif(len(fnames) == 2):
        rQ30 = [0, 0]
        rLen = [1, 1]
        logger.info("Parsing two fastq files.")
        readPhredDF = pd.DataFrame()
        with open(fnames[0].name, 'r') as fastq:
            if(choice2 == 'Y'):
                logger.info("Starting calculations of Q30 for " + getIsolateStr(fnames[0].name) + ".")
            for record in SeqIO.parse(fastq, "fastq"):
                if(len(record.seq) > 0):
                    forwardAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
                    if(choice2 == 'Y'):
                        rLen[0] = rLen[0] + len(record.seq)
                        j = 0
                        while( j < len(record.seq)):
                            if(record.letter_annotations["phred_quality"][j] >= 30):
                                rQ30[0] = rQ30[0] + 1
                            j = j + 1
            if(choice2 == 'Y'):
                logger.info("Completed calculations of Q30 for " + getIsolateStr(fnames[0].name) + ".")
        ## Convert read Phred scores list to pandas Series
        readPhrDF = pd.Series(forwardAvg)
        with open(fnames[1].name, 'r') as fastq2:
            if(choice2 == 'Y'):
                logger.info("Starting calculations of Q30 for " + getIsolateStr(fnames[1].name) + ".")
            for record in SeqIO.parse(fastq2, "fastq"):
                if(len(record.seq) > 0):
                    reverseAvg.append(statistics.mean(record.letter_annotations["phred_quality"]))
                    if(choice2 == 'Y'):
                        rLen[1] = rLen[1] + len(record.seq)
                        j = 0
                        while( j < len(record.seq)):
                            if(record.letter_annotations["phred_quality"][j] >= 30):
                                rQ30[1] = rQ30[1] + 1
                            j = j + 1
            if(choice2 == 'Y'):
                logger.info("Completed calculations of Q30 for " + getIsolateStr(fnames[1].name) + ".")
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
        ## Convert readPhred2 list to pandas Series
        readPhrDF2 = pd.Series(reverseAvg)
        ## Create pandas DataFrame with truncated fnames as column headers
        readPhredDF = extractData(fnames)
        #readPHREDDF = pd.DataFrame({fnamesl1[0] : readLenDF, fnamesl2[0] : readLenDF2})
        ## Coerce string 'Nan' to np.nan
        dfCols = list(readPhredDF)
        for c in dfCols:
            fileMean.append(readPhredDF[c].mean())
            fileMedian.append(readPhredDF[c].median())

        if(choice1 == 'S'):
            for cols in dfCols:
                print("%s total average PHRED is %0.2f" % (cols, readPhredDF[cols].mean() ))
            if(choice2 == 'Y'):
                ii = 0
                for entry in rQ30:
                    print("%s PHRED scores above Q30 is %0.2f" % (dfCols[ii], (100*entry/rLen[ii])))
                    ii = ii + 1
        elif(choice1 == 'L'):
            dfCols = list(readPhredDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in readPhredDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
            ## Call function to plot two violin plots
            logger.info("Plotting series from two fastq files.")
            plotDoubleReadFiles(readPhredDF, choice2, title, rQ30, rLen, fileMean, fileMedian)
## for three or four files, invoke plotMultiReadFiles
    elif((len(fnames) > 2) and (len(fnames) < 7)):
        ii = 0
        rQ30 = [0, 0, 0, 0]
        rLen = [1, 1, 1, 1]
        logger.info("Parsing " + str(len(fnames)) + " fastq files.")
        readPhredDF = pd.Series(dtype=int)
        readPhredDF2 = pd.Series(dtype=int)
        readPHREDDF = pd.DataFrame()
        simpFileName = []
        while(ii < len(fnames)):
            readPHRED = []
            if(choice2 == 'Y'):
                logger.info("Starting calculations of Q30 for " + getIsolateStr(fnames[ii].name) + ".")
            with open(fnames[ii].name, 'r') as fastq:
                for record in SeqIO.parse(fastq, "fastq"):
                    if(len(record.seq) > 0):
                        readPHRED.append(statistics.mean(record.letter_annotations["phred_quality"]))
                        if(choice2 == 'Y'):
                            rLen[ii] = rLen[ii] + len(record.seq)
                            j = 0
                            while( j < len(record.seq)):
                                if(record.letter_annotations["phred_quality"][j] >= 30):
                                    rQ30[ii] = rQ30[ii] + 1
                                j = j + 1
            if(choice2 == 'Y'):
                logger.info("Completed calculations of Q30 for " + getIsolateStr(fnames[ii].name) + ".")
            ## Convert readLengths list to pandas Series
            readPhredDF = pd.Series(readPHRED)
            simpFileName.append(getIsolateStr(fnames[ii].name))
            fnamesl1 = simpFileName[ii].split('.')
            if(re.search('R1', simpFileName[ii], re.IGNORECASE)):
                fnamesl1[0] = fnamesl1[0] + '_R1'
            elif(re.search('R2', simpFileName[ii], re.IGNORECASE)):
                fnamesl1[0] = fnamesl1[0] + '_R2'
            readPHREDDF[fnamesl1[0]] = readPhredDF
            ii = ii + 1

        dfCols = list(readPHREDDF)
        for c in dfCols:
            fileMean.append(readPHREDDF[c].mean())
            fileMedian.append(readPHREDDF[c].median())

        if(choice1 == 'S'):
            for cols in dfCols:
                print("%s total average PHRED is %0.2f" % (cols, readPHREDDF[cols].mean() ))
            if(choice2 == 'Y'):
                jj = 0
                while(jj < len(dfCols)):
                    #print(jj)
                    print("%s PHRED scores above Q30 is %0.2f" % (dfCols[jj], (100*rQ30[jj]/rLen[jj])))
                    jj = jj + 1
                
        elif(choice1 == 'L'):
            dfCols = list(readPHREDDF)
            for header in dfCols:
                print("%s\t" % (header), end="")
            print()
            ii = 0
            for index, row in readPHREDDF.iterrows():
                for cell in dfCols:
                    print(str(row[cell]) + "\t", end="")
                print()
        else:
        ## Call function to plot more than two boxplots
            logger.info("Plotting series from " + str(len(fnames)) + " fastq files.")
            #print(list(readPHREDDF))
            plotMultiReadFiles(readPHREDDF, choice2, title, rQ30, rLen, fileMean, fileMedian)
            
## Exit on error if input files exceeds 6
if(len(args.filename) > 6):
    sys.exit("BioPython_LengthMultiFastq_Reads.py accepts no more than six fastq files as input.")

verifyUniqueFiles = [getIsolateStr(f.name) for f in args.filename]

#Check for redundant/duplicated file input, print(verifyUniqueFiles)
if(len(verifyUniqueFiles) > len(set(verifyUniqueFiles))):
    sys.exit("Identical filenames detected, check input.")

## Invoke function that counts NGS reads for each file in command-line arguments
PHREDOfFastqReads(args.filename, args.titleString, args.outputType, args.showQ30, logger)

