import os, sys, re, csv, statistics
import argparse, logging, warnings

## import numpy
import numpy as np

## import matplotlib
from matplotlib import pyplot as plt

## import pandas
import pandas as pd

## import defaultdict
from collections import defaultdict

def ext_check(expected_ext1, expected_ext2, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2) ):
                        raise ValueError()
                return openner(filename)
        return extension

def getIsolateStr(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	fileString = splitStr[fileNameIdx]
	return fileString

parser = argparse.ArgumentParser(description='Show life expectancy statistics or histogram from US_A_USALEEP.csv or equivalent life table .csv or .txt', usage='Pandas_Parsing_CSV.py filepath/data_table.csv')

parser.add_argument('filename', nargs='+', type=ext_check('.csv', '.txt', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'P'], help="--outputType S for simple statistics and --outputType P for single histogram plot")

args = parser.parse_args()

iter = 0

## Read the life expectancy column, 'e(0)'
def CSV_File_Processor(fname, choice):
    columns = defaultdict(list)
    finExpect = []
    with open(fname[0].name) as csvfile:
        lineReader = csv.DictReader(csvfile, delimiter=',')
        for row in lineReader:
            for (k,v) in row.items():
                columns[k].append(v)
    lifeExpect = columns['e(0)']
    finExpect = [float(i) for i in lifeExpect]
    return(finExpect)

def plotSingleHistOfExp(readLengths, mean, stdDev, fileStr):
    SMALL_SIZE = 30
    MEDIUM_SIZE = 36
    BIG_SIZE = 38
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(15,16))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    smean = "{:.2f}".format(mean)
    sstdDev = "{:.2f}".format(stdDev)
    annotStr = "mean = " + smean + ", sd = " + sstdDev
    axes1.text(58, 5000, annotStr, fontsize=24)
    axes1.hist(readLengths, bins = 50, color='green')
    axes1.set_title('U.S. ' + fileStr + ' (2010-2015)', fontsize = BIG_SIZE)
    axes1.set(ylabel='Census Tract Counts')
    axes1.set(xlabel='Age in Years')
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/lifeExpHist_' + fileStr + '.png')


lifeExpectancy = CSV_File_Processor(args.filename, args.outputType)

## compute mean
mean = statistics.mean(lifeExpectancy)

## comput standard deviation
stdDev = statistics.stdev(lifeExpectancy)

if(args.outputType == 'S'):  ## User selection by --outputType
    print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean, stdDev))
else:
    plotSingleHistOfExp(lifeExpectancy, mean, stdDev, "Life Expectancy")

