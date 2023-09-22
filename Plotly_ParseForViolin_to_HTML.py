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

## import plotly
import plotly.express as px
#import plotly.plotly as plt
import plotly as plt
import plotly.graph_objects as go
from plotly.offline import iplot

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

#parser.add_argument('--outputType', '-o', default='S', choices=['S', 'P'], help="--outputType S for simple statistics and --outputType P for single histogram plot")

#parser.add_argument('--statePlots', '-s', default='3', choices=['0', '1', '2', '3'], help="--statePlots count, 0 for total only, 1 for total and one state, 2 for total and two states, 3 for total and three states")

parser.add_argument('--titleString', '-t', default='United_States', help="--outputType S for simple statistics and --outputType P for single histogram plot")

parser.add_argument('--stateCodes', '-c', nargs='+', type=str, default=['13'], help='--stateCodes expects one to three integers separated by a space, integer must be less than 57.')

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


## Read the life expectancy column, 'e(0)'
def simple_CSV_File_Processor(fname):
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
    #print(len(firstFinExpect))
    
    if(len(stateCodes) > 1):
        #print(len(secondFinExpect))
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
        #print(len(thirdFinExpect))
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


## plot a single violin plot of life expectancy with interactive annotation
def plotSingleViolinLifeExp(lifeExp, titleStr):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    fig = go.Figure(data=go.Violin(y=lifeExp['USA'], name='USA', box_visible=True, line_color='black', meanline_visible=True, fillcolor='darkseagreen'))
    fig.update_layout(title="U.S. Life Expectancy (2010-2015)", xaxis_title="Census Tracts", yaxis_title="Age in Years")
    plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyViolin_" + titleStr + ".html")

## plot two violin plots of life expectancies, U.S. plus 1 state, with interactive annotation
def plotDoubleViolinLifeExp(allLifeExp, stateLifeExp, titleStr, codes, stateCodesDict):
    bothDF = pd.concat([allLifeExp, stateLifeExp])
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    fig = go.Figure()
    colors = ['darkseagreen', 'goldenrod']
    ii = 0
    for item in list(bothDF):
        fig = fig.add_trace(go.Violin(y=bothDF[item], name=item, box_visible=True, line_color='black', meanline_visible=True, fillcolor=colors[ii]))
        ii = ii + 1
    fig.update_layout(title="U.S. and Georgia Life Expectancy (2010-2015)", xaxis_title="Census Tract Groupings", yaxis_title="Age in Years")
    plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyViolin_1state_" + titleStr + ".html")

## plot three violin plots of life expectancies, U.S. plus 2 states, with interactive annotation
def plotTripleViolinLifeExp(allLifeExp, stateLifeExp, titleStr, codes, stateCodesDict):
    bothDF = pd.concat([allLifeExp, stateLifeExp])
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    fig = go.Figure()
    colors = ['darkseagreen', 'yellow', 'orange']
    ii = 0
    for item in list(bothDF):
        tempSeries = bothDF[item].squeeze()
        fig = fig.add_trace(go.Violin(y=tempSeries[tempSeries>0.0], name=item, box_visible=True, line_color='black', meanline_visible=True, fillcolor=colors[ii]))
        ii = ii + 1
    fig.update_layout(title="%s Life Expectancy (2010-2015)" % (titleStr), xaxis_title="Census Tract Groupings", yaxis_title="Age in Years")
    plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyViolin_2state_" + titleStr + ".html")
    

## plot four violin plots of life expectancies, U.S. plus 3 states, with interactive annotation
def plotQuadrupleViolinLifeExp(allLifeExp, stateLifeExp, titleStr, codes, stateCodesDict):
    bothDF = pd.concat([allLifeExp, stateLifeExp])
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    fig = go.Figure()
    colors = ['darkseagreen', 'green', 'yellow', 'orange']
    ii = 0
    for item in list(bothDF):
        tempSeries = bothDF[item].squeeze()
        fig = fig.add_trace(go.Violin(y=tempSeries[tempSeries>0.0], name=item, box_visible=True, line_color='black', meanline_visible=True, fillcolor=colors[ii]))
        ii = ii + 1
    fig.update_layout(title="%s Life Expectancy (2010-2015)" % (titleStr), xaxis_title="Census Tract Groupings", yaxis_title="Age in Years")
    plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyViolin_2state_" + titleStr + ".html")


allLifeExpectancy = simple_CSV_File_Processor(args.filename)

statesLifeExpect = state_CSV_File_Processor(args.filename, codes, stateCodesDict)

#print(statesLifeExpect.tail())

## Decision statements for number of state codes supplied by user
if(len(codes) == 0):
    plotSingleViolinLifeExp(allLifeExpectancy, args.titleString)
elif(len(codes) == 1):
    plotDoubleViolinLifeExp(allLifeExpectancy, statesLifeExpect, args.titleString, codes, stateCodesDict)
elif(len(codes) == 2):
    plotTripleViolinLifeExp(allLifeExpectancy, statesLifeExpect, args.titleString, codes, stateCodesDict)
elif(len(codes) == 3):
    plotQuadrupleViolinLifeExp(allLifeExpectancy, statesLifeExpect, args.titleString, codes, stateCodesDict)
else:
    print("EXCEPTION: Current version of Plotly_ParseForViolin_to_HTML.py not able to show more than 4 violin plots.\nPlease select one to three numbers between 1 and 56.\n.\n.\n.")
    sys.exit()


