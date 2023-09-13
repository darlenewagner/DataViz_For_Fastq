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

parser = argparse.ArgumentParser(description='Show life expectancy from US_A_USALEEP.csv or equivalent life table .csv or .txt', usage='Pandas_Parsing_CSV.py filepath/data_table.csv')

parser.add_argument('filename', nargs='+', type=ext_check('.csv', '.txt', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='S', choices=['S', 'C'], help="--outputType S for simple statistics and --outputType C for condensed table")

args = parser.parse_args()

iter = 0

#myFastq = open(args.filename[0].name, "r")

def CSV_File_Processor(fname, choice):
    lineCount = 0
    with open(fname[0].name) as csvfile:
        lineReader = csv.reader(csvfile, delimiter=',')
        csv_headers = next(lineReader)
        first_line = next(lineReader)
    #print("%s has %i reads" % (getIsolateStr(fname), readCount))
    print(csv_headers)

CSV_File_Processor(args.filename, args.outputType)




        
    
