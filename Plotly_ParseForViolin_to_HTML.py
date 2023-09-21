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

#parser.add_argument('--titleString', '-t', default='United States', help="--outputType S for simple statistics and --outputType P for single histogram plot")

#parser.add_argument('--stateCodes', '-c', nargs='+', type=str, default=['13', '54'], help='--stateCodes expects one to three integers separated by a space, integer must be less than 57.')

args = parser.parse_args()

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

## plot a single violin plot of life expectancy annotated with mean and standard deviation
def plotSingleViolinLifeExp(lifeExp):
    SMALL_SIZE = 32
    MEDIUM_SIZE = 36
    BIG_SIZE = 40
    colors = ['#00FF00']
    fig = go.Figure(data=go.Violin(y=lifeExp, box_visible=True, line_color='black', meanline_visible=True, fillcolor='darkseagreen'))
    fig.update_layout(title="U.S. Life Expectancy (2010-2015)", xaxis_title="Jurisdiction", yaxis_title="Age in Years")
    plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyViolin.html")


allLifeExpectancy = simple_CSV_File_Processor(args.filename)

plotSingleViolinLifeExp(allLifeExpectancy)

#fig = px.line(x=[1,2,3,4], y=[1,2,3,5])
#plt.offline.plot(fig, filename="/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/Temp/Darlene/plotlyLine.html")

