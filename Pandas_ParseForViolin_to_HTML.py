import os, sys, re, csv, statistics
import argparse, logging, warnings
from itertools import repeat

## import numpy
import numpy as np, seaborn as sns

## import matplotlib and supporting packages
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib import image, patches, colors
from matplotlib.colors import colorConverter

## import pandas
import pandas as pd

## import defaultdict
from collections import defaultdict

## import mpld3 to embed figure in a webpage
import mpld3

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

parser.add_argument('--stateCodes', '-c', nargs='+', type=str, default=['13', '54'], help='--stateCodes expects one to three integers separated by a space, integer must be less than 57.')

args = parser.parse_args()

## Dictionary to convert stateCodes to state names
stateCodesDict = {'1':'Alabama', '2':'Alaska', '4':'Arizona', '5':'Arkansas', '6':'California', '8':'Colorado',
                  '9':'Connecticut', '10':'Delaware', '11':'District of Columbia', '12':'Florida', '13':'Georgia', '15':'Hawaii',
                  '16':'Idaho', '17':'Illinois', '18':'Indiana', '19':'Iowa', '20':'Kansas', '21':'Kentucky', '22':'Louisiana',
                  '23':'Maine', '24':'Maryland', '25':'Massachusetts', '26':'Michigan', '27':'Minnesota', '28':'Mississippi',
                  '29':'Missouri', '30':'Montana', '31':'Nebraska', '33':'New Hampshire', '34':'New Jersey', '35':'New Mexico',
                  '36':'New York', '37':'North Carolina', '38':'North Dakota', '39':'Ohio', '40':'Oklahoma', '41':'Oregon',
                  '42':'Pennsylvania', '44':'Rhode Island', '45':'South Carolina', '46':'South Dakota', '47':'Tennessee',
                  '48':'Texas', '49':'Utah', '50':'Vermont', '51':'Virginia', '53':'Washington', 
                  '54':'West Virginia', '55':'Wisconsin', '56':'Wyoming' }

## Check for existence of --stateCodes and stateCodesDict
codes = args.stateCodes
for code in codes:
    if(code not in stateCodesDict):
        print("Error: No matching state code in %s ! Try an integer between 1 and 56." % (getIsolateStr(args.filename[0].name)))
        sys.exit()

iter = 0

def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))


#def make_violinplot(x1, x2, simpDataFrame, my_labels, my_colors):
    #horizontal_concat = pd.concat([x1, x2], axis=1)
    #fig1 = plt.subplots(figsize=(26,19))
    #print(list(simpDataFrame))
    #fig = sns.violinplot(data=simpDataFrame, cut=0, inner='box', showmeans=True, showmedians=True).get_figure
    #fig.tick_params(axis='x', labelsize=10)
    #fig.set_xticks(np.arange(1, len(['USA', 'DC'])+1), ['USA', 'DC'])
    #fig.set_xticklabels(['USA', 'DC'])

    #for body, color in zip(parts['bodies'], my_colors):
    #    body.set_facecolor(color)
    #    parts['cmeans'].set_color('red')
    #    parts['cmedians'].set_color('blue')
    
    #


    #axes1.title("U.S. and Washington D.C. Life Expectancy, (2010-2015)")
    #plt.rcParams["figure.figsize"] = [8,8]
    #return fig

## Read the life expectancy column, 'e(0)'
def simple_CSV_File_Processor(fname, choice):
    columns = defaultdict(list)
    finExpect = []
    finExpectDataFrame = pd.DataFrame()
    with open(fname[0].name) as csvfile:
        lineReader = csv.DictReader(csvfile, delimiter=',')
        for row in lineReader:
            for (k,v) in row.items():
                columns[k].append(v)
    lifeExpect = columns['e(0)']
    finExpect = [float(i) for i in lifeExpect]
    finExpectDataFrame['USA'] = finExpect
    return(finExpectDataFrame)

## Read the life expectancy column and break down by stateCodes
def state_CSV_File_Processor(fname, stateCodes, stateCodesDict):
    columns = defaultdict(list)
    firstExpect = []
    firstFinExpect = []
    secondExpect = []
    secondFinExpect = []
    thirdExpect = []
    thirdFinExpect = []
    
    expectDataFrame = pd.DataFrame()
    with open(fname[0].name) as csvfile:
        csvfile.seek(0)
        lineReader = csv.DictReader(csvfile, delimiter=',')
        for row in lineReader:
            for (k,v) in row.items():
                if(row['STATE2KX'] == stateCodes[0]):
                    columns[k].append(v)
        firstExpect = columns['e(0)']        
        firstFinExpect = [float(i) for i in firstExpect]
        columns = defaultdict(list)
        if(len(stateCodes) > 1):
            csvfile.seek(0)
            for row in lineReader:
                for (k,v) in row.items():
                    if(row['STATE2KX'] == stateCodes[1]):
                        columns[k].append(v)
            secondExpect = columns['e(0)']        
            secondFinExpect = [float(i) for i in secondExpect]
        columns = defaultdict(list)
        if(len(stateCodes) > 2):
            csvfile.seek(0)
            for row in lineReader:
                for (k,v) in row.items():
                    if(row['STATE2KX'] == stateCodes[2]):
                        columns[k].append(v)
            thirdExpect = columns['e(0)']
            thirdFinExpect = [float(i) for i in thirdExpect]

    print(len(firstFinExpect))
    
    if(len(stateCodes) > 1):
        print(len(secondFinExpect))
        if(len(secondFinExpect) > len(firstFinExpect)):
            extend_length = len(secondFinExpect) - len(firstFinExpect)
            pad = 0
            while(pad < extend_length):
                firstFinExpect.append(0)
                pad = pad + 1
        if(len(secondFinExpect) < len(firstFinExpect)):
            extend_length = len(firstFinExpect) - len(secondFinExpect)
            pad = 0
            while(pad < extend_length):
                secondFinExpect.append(0)
                pad = pad + 1
    if(len(stateCodes) > 2):
        print(len(thirdFinExpect))
        if(len(thirdFinExpect) > len(secondFinExpect)):
            extend_length = len(thirdFinExpect) - len(secondFinExpect)
            pad = 0
            while(pad < extend_length):
                firstFinExpect.append(0)
                secondFinExpect.append(0)
                pad = pad + 1
        if(len(thirdFinExpect) < len(secondFinExpect)):
            extend_length = len(secondFinExpect) - len(thirdFinExpect)
            pad = 0
            while(pad < extend_length):
                thirdFinExpect.append(0)
                pad = pad + 1
                
    expectDataFrame[stateCodesDict[stateCodes[0]]] = firstFinExpect
    if(len(stateCodes) > 1):
        expectDataFrame[stateCodesDict[stateCodes[1]]] = secondFinExpect
    if(len(stateCodes) > 2):
        expectDataFrame[stateCodesDict[stateCodes[2]]] = thirdFinExpect
        
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
    vp = sns.violinplot(lifeExp, showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violinLifeExp_' + titleStr + '.png')

## Function to plot read lengths of two life expectancy series
def plotDoubleViolinLifeExp(allLifeExp, stateLifeExp, mean, stdDev, fileStr, titleStr, stateCodes, stateCodesDict):
    #plt.style.use('seaborn')
    #my_labels = ['U.S.', 'District of Columbia']
    #my_colors = ['Green', 'Yellow']
    SMALL_SIZE = 34
    MEDIUM_SIZE = 40
    BIG_SIZE = 44
    colors = ['#00FF00', '#FFFF00', '#FFA500', '#FF0000']
    fig1 = plt.figure(figsize=(25,19))
    axes1 = plt.subplot(111)
    axes1.xaxis.label.set_size(SMALL_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=MEDIUM_SIZE)
    axes1.tick_params(axis='y', labelsize=MEDIUM_SIZE)
    tmean = "{:.2f}".format(mean[0])
    fmean = "{:.2f}".format(mean[1])
    annotStr = "U.S. mean = " + tmean + " and " + stateCodesDict[stateCodes[0]]  + " mean = " + fmean 
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = ['U.S.', stateCodesDict[stateCodes[0]]]
    pos = [0, 1]
    green_patch = mpatches.Patch(color=colors[0])
    yellow_patch = mpatches.Patch(color=colors[1])
    fake_handles = [green_patch, yellow_patch]
    vp = plt.violinplot([allLifeExp['USA'], stateLifeExp[stateCodesDict[stateCodes[0]]]], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
    #vp.xticks([0, 1], x_labels)
    font_dict = {'fontsize': 20, 'fontweight': 'heavy', 'verticalalignment':'top'}
    #axes1.legend(x_labels)
    axes1.legend(fake_handles, x_labels, bbox_to_anchor=(1,0), prop={'size': SMALL_SIZE})
    axes1.set_xticks([])
    axes1.set_xticklabels([])
    
    xy = [[l.vertices[:,0].mean(),l.vertices[0,1]] for l in vp['cmeans'].get_paths()]
    xy = np.array(xy)
    axes1.scatter(xy[:,0], xy[:,1],s=121, c="#A020F0", marker="D", zorder=3)
    for pc, color in zip(vp['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_linewidth(4)
        pc.set_alpha(0.70)
    vp['cmeans'].set_visible(False)
    vp['cmedians'].set_colors('black')
    vp['cmedians'].set_linewidth(3)
#    
#    vp.set_xticklabels(x_labels)
    axes1.set_title( titleStr + ' (2010-2015)', fontsize = BIG_SIZE)
    axes1.set(ylabel='Age in Years')
    axes1.set_xlabel(annotStr, ha='right')
    ##fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violinLifeExp_1state_' + titleStr + '.png')
    #htm_str = mpld3.fig_to_html(make_violinplot(allLifeExp['USA'], stateLifeExp[stateCodesDict[stateCodes[0]]], pd.concat([allLifeExp, stateLifeExp], axis=0), my_labels, my_colors))
    htm_str = mpld3.fig_to_html(fig1)
    htm_file = open("/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/seabornFig_1state.html", "w")
    htm_file.write(htm_str)
    htm_file.close()


## Function to plot read lengths of three life expectancy series
def plotTripleViolinLifeExp(allLifeExp, stateLifeExp1, stateLifeExp2, mean, stdDev, fileStr, titleStr, stateCodes, stateCodesDict):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00', '#006400', '#FFFF00', '#FFA500']
    fig1, axes1 = plt.subplots(figsize=(26,19))
    axes1.xaxis.label.set_size(SMALL_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=MEDIUM_SIZE)
    axes1.tick_params(axis='y', labelsize=MEDIUM_SIZE)
    tmean = "{:.2f}".format(mean[0])
    fmean = "{:.2f}".format(mean[1])
    smean = "{:.2f}".format(mean[2])
    annotStr = "U.S. mean = " + tmean + ", " + stateCodesDict[stateCodes[0]]  + " mean = " + fmean + ", and " + stateCodesDict[stateCodes[1]] + " mean = " + smean
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
    ## Old command for saving figure without embedding in html
    ## fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violinLifeExp_2state_' + titleStrMod + '.png')
    htm_str = mpld3.fig_to_html(fig1)
    htm_file = open("/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/indexFig_2state.html", "w")
    htm_file.write(htm_str)
    htm_file.close()
    

## Function to plot read lengths of three life expectancy series
def plotQuadrupleViolinLifeExp(allLifeExp, stateLifeExp1, stateLifeExp2, stateLifeExp3, mean, stdDev, fileStr, titleStr, stateCodes, stateCodesDict):
    SMALL_SIZE = 30
    MEDIUM_SIZE = 36
    BIG_SIZE = 38
    colors = ['#00FF00', '#006400', '#FFA500', '#FF0000']
    fig1, axes1 = plt.subplots(figsize=(27,19))
    axes1.xaxis.label.set_size(SMALL_SIZE)
    axes1.yaxis.label.set_size(MEDIUM_SIZE)
    axes1.tick_params(axis='x', labelsize=MEDIUM_SIZE)
    axes1.tick_params(axis='y', labelsize=MEDIUM_SIZE)
    omean = "{:.2f}".format(mean[0])
    fmean = "{:.2f}".format(mean[1])
    smean = "{:.2f}".format(mean[2])
    tmean = "{:.2f}".format(mean[3])
    annotStr = "U.S. mean = " + omean + ", " + stateCodesDict[stateCodes[0]]  + " mean = " + fmean + ", " + stateCodesDict[stateCodes[1]] + " mean = " + smean + ", and " + stateCodesDict[stateCodes[2]] + " mean = " + tmean
    medianprops = dict(linewidth=6, color='black')
    whiskerprops = dict(linewidth=5, color='black')
    capprops = dict(linewidth=5, color='black')
    x_labels = ['U.S.', stateCodesDict[stateCodes[0]], stateCodesDict[stateCodes[1]], stateCodesDict[stateCodes[2]] ]
    vp = axes1.violinplot([allLifeExp, stateLifeExp1, stateLifeExp2, stateLifeExp3], showmedians=True, showmeans=True, widths=0.95, showextrema=False)
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
    #htm = fig1.savefig('/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/violinLifeExp_3state_' + titleStrMod + '.png')
    htm_str = mpld3.fig_to_html(fig1)
    htm_file = open("/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/indexFig_3state.html", "w")
    htm_file.write(htm_str)
    htm_file.close()

allLifeExpectancy = simple_CSV_File_Processor(args.filename, args.outputType)

statesLifeExpect = state_CSV_File_Processor(args.filename, codes, stateCodesDict)

#print(statesLifeExpect.tail())

mean = []
stdDev = []
temp0Series =[]
temp1Series = []
temp2Series = []
temp3Series = []

## compute mean
temp0Series = allLifeExpectancy['USA']
mean.append(statistics.mean(temp0Series[temp0Series!=0]))
temp1Series = statesLifeExpect[stateCodesDict[codes[0]]]
mean.append(statistics.mean(temp1Series[temp1Series!=0]))
if(len(codes) > 1):
    temp2Series = statesLifeExpect[stateCodesDict[codes[1]]]
    mean.append(statistics.mean(temp2Series[temp2Series!=0]))
if(len(codes) > 2):
    temp3Series = statesLifeExpect[stateCodesDict[codes[2]]]
    mean.append(statistics.mean(temp3Series[temp3Series!=0]))

## compute standard deviation
stdDev.append(statistics.stdev(temp0Series[temp0Series!=0]))
stdDev.append(statistics.stdev(temp1Series[temp1Series!=0]))
if(len(codes) > 1):
    stdDev.append(statistics.stdev(temp2Series[temp2Series!=0]))
if(len(codes) > 2):
    stdDev.append(statistics.stdev(temp3Series[temp3Series!=0]))

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
    elif(len(codes) == 3):
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (getIsolateStr(args.filename[0].name), mean[0], stdDev[0]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (stateCodesDict[codes[0]], mean[1], stdDev[1]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (stateCodesDict[codes[1]], mean[2], stdDev[2]))
        print("Mean life expectancy from %s is %0.2f and standard deviation is %0.2f" % (stateCodesDict[codes[2]], mean[3], stdDev[3]))
else:
    if(len(codes) == 0):
        plotSingleViolinLifeExp(allLifeExpectancy, mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString)
    elif(len(codes) == 1):
        plotDoubleViolinLifeExp(allLifeExpectancy, statesLifeExpect, mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString, codes, stateCodesDict)
    elif(len(codes) == 2):
        plotTripleViolinLifeExp(allLifeExpectancy, temp1Series[temp1Series!=0], temp2Series[temp2Series!=0], mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString, codes, stateCodesDict)
    elif(len(codes) == 3):
        plotQuadrupleViolinLifeExp(allLifeExpectancy, temp1Series[temp1Series!=0], temp2Series[temp2Series!=0], temp3Series[temp3Series!=0], mean, stdDev, getIsolateStr(args.filename[0].name), args.titleString, codes, stateCodesDict)
