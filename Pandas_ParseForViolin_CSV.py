import os, sys, re, csv, statistics
import argparse, logging, warnings

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

## Data Frame of life expectancies broken down by state
statesLifeExpect = pd.DataFrame()

parser = argparse.ArgumentParser(description='Show life expectancy statistics or histogram from US_A_USALEEP.csv or equivalent life table .csv or .txt', usage='Pandas_Parsing_CSV.py filepath/data_table.csv')

parser.add_argument('filename', nargs='+', type=ext_check('.csv', '.txt', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'P'], help="--outputType S for simple statistics and --outputType P for single histogram plot")

parser.add_argument('--statePlots', '-s', default='0', choices=['0', '1', '2', '3', '4'], help="--statePlots count, 0 for total only, 1 for total and one state, 2 for two states, 3 for three states, and 4 for four states")

#parser.add_argument('--stateCode', '-o', default='N', choices=['Y', 'N'], help="--outputType S for simple statistics and --outputType P for single histogram plot")

args = parser.parse_args()

iter = 0

## Read the life expectancy column, 'e(0)'
def simple_CSV_File_Processor(fname, choice):
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

## Read the life expectancy column and break down my state code
def state_CSV_File_Processor(fname, states):
    columns = defaultdict(list)
    firstExpect = []
    secondExpect = []
    thirdExpect = []
    fourthExpect = []
    ##expectDataFrame = pd.DataFrame()
    with open(fname[0].name) as csvfile:
        lineReader = csv.DictReader(csvfile, delimiter=',')
        for row in lineReader:
            for (k,v) in row.items():
                if(row['STATE2KX'] == '13'):
                    #print(row)
                    columns[k].append(v)
    firstExpect = columns['e(0)']        
    finExpect = [float(i) for i in firstExpect]
    expectDataFrame = pd.DataFrame(finExpect, columns=['Georgia'])
    return(expectDataFrame)


## plot a single violin plot annotated with mean and standard deviation
def plotSingleViolinLifeExp(lifeExp, mean, stdDev, fileStr):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00', '#FFFF00', '#FFA500', '#FF0000']
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(15,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    smean = "{:.2f}".format(mean)
    sstdDev = "{:.2f}".format(stdDev)
    annotStr = "mean = " + smean + ", sd = " + sstdDev
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [fileStr]
    vp = axes1.violinplot(lifeExp, showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    axes1.set_title('U.S. ' + fileStr + ' (2010-2015)', fontsize = BIG_SIZE)
    axes1.set(ylabel='Age in Years')
    axes1.set(xlabel=annotStr)
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLifeExpr_' + fileStr + '.png')


allLifeExpectancy = simple_CSV_File_Processor(args.filename, args.outputType)

statesLifeExpect = state_CSV_File_Processor(args.filename, args.statePlots)

## compute mean
mean = statistics.mean(allLifeExpectancy)
gmean = statistics.mean(statesLifeExpect['Georgia'])

## comput standard deviation
stdDev = statistics.stdev(allLifeExpectancy)
gstdDev = statistics.stdev(statesLifeExpect['Georgia'])

if(args.outputType == 'S'):  ## User selection by --outputType
    if(args.statePlots == '0'):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean, stdDev))
    elif(args.statePlots == '1'):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean, stdDev))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (list(statesLifeExpect)[0], gmean, gstdDev))
else:
    plotSingleViolinLifeExp(allLifeExpectancy, mean, stdDev, getIsolateStr(args.filename[0].name))

