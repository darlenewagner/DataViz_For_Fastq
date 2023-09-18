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

## Data Frame of life expectancies broken down by state
statesLifeExpect = pd.DataFrame()

parser = argparse.ArgumentParser(description='Show life expectancy statistics or histogram from US_A_USALEEP.csv or equivalent life table .csv or .txt', usage='Pandas_Parsing_CSV.py filepath/data_table.csv')

parser.add_argument('filename', nargs='+', type=ext_check('.csv', '.txt', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'P'], help="--outputType S for simple statistics and --outputType P for single histogram plot")

#parser.add_argument('--statePlots', '-s', default='3', choices=['0', '1', '2', '3'], help="--statePlots count, 0 for total only, 1 for total and one state, 2 for total and two states, 3 for total and three states")

parser.add_argument('--titleString', '-t', default='United States', help="--outputType S for simple statistics and --outputType P for single histogram plot")

parser.add_argument('--stateCodes', '-c', nargs='+', type=str, default=['13', '54'], help='--stateCodes expects one to three integers separated by a space, integer must be less than 55.')

args = parser.parse_args()

## Dictionary to convert stateCodes to state names
stateCodesDict = {'1':'Alabama', '2':'Alaska', '4':'Arizona', '5':'Arkansas', '6':'California', '8':'Colorado',
                  '9':'Connecticut', '10':'Delaware', '11':'District of Columbia', '12':'Florida', '13':'Georgia', '15':'Hawaii',
                  '16':'Idaho', '17':'Illinois', '18':'Indiana', '19':'Iowa', '20':'Kansas', '21':'Kentucky', '22':'Louisiana',
                  '23':'Maine', '24':'Maryland', '25':'Massachusetts', '26':'Michigan', '27':'Minnesota', '28':'Mississippi',
                  '29':'Missouri', '30':'Montana', '54':'West Virginia'}

## Check for existence of --stateCodes and stateCodesDict
codes = args.stateCodes
for code in codes:
    if(code not in stateCodesDict):
        print("Error: No matching state code in %s !" % (getIsolateStr(args.filename[0].name)))
        sys.exit()

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

## Read the life expectancy column and break down by stateCodes
def state_CSV_File_Processor(fname, stateCodes, stateCodesDict):
    columns = defaultdict(list)
    firstExpect = []
    secondExpect = []
    thirdExpect = []
    fourthExpect = []
    expectDataFrame = pd.DataFrame()
    with open(fname[0].name) as csvfile:
        lineReader = csv.DictReader(csvfile, delimiter=',')
        for row in lineReader:
            for (k,v) in row.items():
                if(row['STATE2KX'] == stateCodes[0]):
                    #print(row)
                    columns[k].append(v)
        firstExpect = columns['e(0)']        
        firstFinExpect = [float(i) for i in firstExpect]
        columns = defaultdict(list)
        csvfile.seek(0)
        for row in lineReader:
            for (k,v) in row.items():
                if(row['STATE2KX'] == stateCodes[1]):
                    columns[k].append(v)
        secondExpect = columns['e(0)']        
        secondFinExpect = [float(i) for i in secondExpect]
        columns = defaultdict(list)
        print(len(firstFinExpect))
        print(len(secondFinExpect))
    if(len(secondFinExpect) > len(firstFinExpect)):
        extend_length = len(secondFinExpect) - len(firstFinExpect)
        pad = 0
        while(pad < extend_length):
            firstFinExpect.append(0)
            pad = pad + 1
    expectDataFrame[stateCodesDict[stateCodes[0]]] = firstFinExpect
    
    if(len(secondFinExpect) < len(firstFinExpect)):
        extend_length = len(firstFinExpect) - len(secondFinExpect)
        pad = 0
        while(pad < extend_length):
            secondFinExpect.append(0)
            pad = pad + 1
    expectDataFrame[stateCodesDict[stateCodes[1]]] = secondFinExpect
    return(expectDataFrame)


## plot a single violin plot of life expectancy annotated with mean and standard deviation
def plotSingleViolinLifeExp(lifeExp, mean, stdDev, fileStr, titleStr):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00', '#FFFF00', '#FFA500', '#FF0000']
    fig1, axes1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(15,19))
    axes1.xaxis.label.set_size(MEDIUM_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=SMALL_SIZE)
    axes1.tick_params(axis='y', labelsize=SMALL_SIZE)
    smean = "{:.2f}".format(mean[0])
    sstdDev = "{:.2f}".format(stdDev[0])
    annotStr = "mean = " + smean + ", sd = " + sstdDev
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = [titleStr]
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
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLifeExpr_' + titleStr + '.png')

## Function to plot read lengths of two life expectancy series
def plotDoubleViolinLifeExp(allLifeExp, stateLifeExp, mean, stdDev, fileStr, titleStr, stateCodes, stateCodesDict):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00', '#FFFF00', '#FFA500', '#FF0000']
    fig1, axes1 = plt.subplots(figsize=(25,19))
    axes1.xaxis.label.set_size(SMALL_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=MEDIUM_SIZE)
    axes1.tick_params(axis='y', labelsize=MEDIUM_SIZE)
    tmean = "{:.2f}".format(mean[0])
    fmean = "{:.2f}".format(mean[1])
    annotStr = "U.S. mean = " + tmean + " and Georgia mean = " + fmean
    #fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    #readLength1 = readLengthDF[dfCols[0]]
    #readLength2 = readLengthDF[dfCols[1]]
    #filteredReadLength1 = readLength1[~np.isnan(readLength1)]
    #filteredReadLength2 = readLength2[~np.isnan(readLength2)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = ['U.S.', stateCodesDict[stateCodes[0]]]
    vp = axes1.violinplot([allLifeExp, stateLifeExp], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    axes1.set_title( titleStr + ' (2010-2015)', fontsize = BIG_SIZE)
    axes1.set(ylabel='Age in Years')
    axes1.set(xlabel=annotStr)
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLen_' + titleStr + '.png')


## Function to plot read lengths of two life expectancy series
def plotTripleViolinLifeExp(allLifeExp, stateLifeExp1, stateLifeExp2, mean, stdDev, fileStr, titleStr, stateCodes, stateCodesDict):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00', '#FFFF00', '#FFA500', '#FF0000']
    fig1, axes1 = plt.subplots(figsize=(26,19))
    axes1.xaxis.label.set_size(SMALL_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=MEDIUM_SIZE)
    axes1.tick_params(axis='y', labelsize=MEDIUM_SIZE)
    tmean = "{:.2f}".format(mean[0])
    fmean = "{:.2f}".format(mean[1])
    smean = "{:.2f}".format(mean[2])
    annotStr = "U.S. mean = " + tmean + ", " + stateCodesDict[stateCodes[0]]  + " mean = " + fmean + ", and " + stateCodesDict[stateCodes[1]] + " mean = " + smean
    #fileTitle = dfCols[0] + "_thru_" + dfCols[len(dfCols) - 1]
    #readLength1 = readLengthDF[dfCols[0]]
    #readLength2 = readLengthDF[dfCols[1]]
    #filteredReadLength1 = readLength1[~np.isnan(readLength1)]
    #filteredReadLength2 = readLength2[~np.isnan(readLength2)]
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = ['U.S.', stateCodesDict[stateCodes[0]], stateCodesDict[stateCodes[1]] ]
    vp = axes1.violinplot([allLifeExp, stateLifeExp1, stateLifeExp2], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    axes1.set_title(titleStr + ' (2010-2015)', fontsize = BIG_SIZE)
    axes1.set(ylabel='Age in Years')
    axes1.set(xlabel=annotStr)
    titleStrMod = re.sub(r',', '', titleStr)
    fig1.savefig('/scicomp/home-pure/ydn3/output_of_DataViz_For_Fastq/violinLen_' + titleStrMod + '.png')


    

allLifeExpectancy = simple_CSV_File_Processor(args.filename, args.outputType)

statesLifeExpect = state_CSV_File_Processor(args.filename, codes, stateCodesDict)

print(statesLifeExpect.tail())

mean = []
stdDev = []
temp1Series = []
temp2Series = []
temp3Series = []

## compute mean
mean.append(statistics.mean(allLifeExpectancy))
temp1Series = statesLifeExpect[stateCodesDict[codes[0]]]
mean.append(statistics.mean(temp1Series[temp1Series!=0]))
temp2Series = statesLifeExpect[stateCodesDict[codes[1]]]
mean.append(statistics.mean(temp2Series[temp2Series!=0]))

## compute standard deviation
stdDev.append(statistics.stdev(allLifeExpectancy))
stdDev.append(statistics.stdev(temp1Series[temp1Series!=0]))
stdDev.append(statistics.stdev(temp2Series[temp2Series!=0]))


if(args.outputType == 'S'):  ## User selection by --outputType
    if(len(codes) == 0):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (args.titleString, mean[0], stdDev[0]))
    elif(len(codes) == 1):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean[0], stdDev[0]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (list(statesLifeExpect)[0], mean[1], stdDev[1]))
    elif(len(codes) == 2):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean[0], stdDev[0]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (list(statesLifeExpect)[0], mean[1], stdDev[1]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (list(statesLifeExpect)[1], mean[2], stdDev[2]))
else:
    if(len(codes) == 0):
        plotSingleViolinLifeExp(allLifeExpectancy, mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString)
    elif(len(codes) == 1):
        plotDoubleViolinLifeExp(allLifeExpectancy, temp1Series[temp1Series!=0], mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString)
    elif(len(codes) == 2):
        plotTripleViolinLifeExp(allLifeExpectancy, temp1Series[temp1Series!=0], temp2Series[temp2Series!=0], mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString, codes, stateCodesDict)
