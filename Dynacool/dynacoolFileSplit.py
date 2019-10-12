#!/usr/bin/python

#This program takes a raw Dynacool resistance measurement .dat file and
#   1. separates non-empty columns into files for channel 1 and channel 2 (separate samples)

import sys
import string
import os.path
import numpy as np

#READ IN ARGUMENTS - [filename]
#first argument is the filename to read in
filename = sys.argv[1]
head, tail = os.path.splitext(filename)

print "Reading File: " + filename
print "Unstacking Channel 1 and Channel 2"

#SEPARATE HEADER FROM DATA
headerline = ''
header = '' #header will be added to the beginning of any output files
nline = 0
inputfile = open(filename,'rU')
while headerline != "[Data]\n":
    headerline = inputfile.readline()
    nline += 1
    header = header + headerline
columnLabels = inputfile.readline().strip().split(',') #list of columns, as formatted in the text file (with parentheses and spaces)
inputfile.close()

#READ IN DATA FILE, masked. Not all of them will be filled. The mask has values of True where data is missing
data = np.genfromtxt(filename, delimiter=',', skip_header=nline, names=True, usemask=True, dtype=None)
notmissing = np.logical_not(data.mask.view(bool).reshape(len(data),-1)) #turns the mask (an array of tuples) into a 2D array where True is a valid cell

#FIGURE OUT WHICH COLUMNS TO WRITE OUT
#find empty columns
columnsWithData = np.logical_or.reduce(notmissing)
#find out which columns belong to which channel
columnNames = np.array(data.dtype.names)
columnsInCh1 = np.array([a.find('Ch2') == -1 for a in columnNames])
np.logical_and(columnsWithData, columnsInCh1, columnsInCh1)
columnsInCh2 = np.array([a.find('Ch1') == -1 for a in columnNames])
np.logical_and(columnsWithData, columnsInCh2, columnsInCh2)
#find out which columns to unstack
columnsToSplit = np.logical_xor(columnsInCh1, columnsInCh2)
columnsToKeep = np.logical_and(columnsInCh1, columnsInCh2)
#figure out which columns are unique to a specific measurement
colsUniquetoCh1 = np.logical_and(columnsInCh1, np.logical_not(columnsInCh2))
colsUniquetoCh2 = np.logical_and(columnsInCh2, np.logical_not(columnsInCh1))

#FIGURE OUT WHICH ROWS CORRESPOND TO WHICH COLUMNS
rowsInCh1 = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh1), axis=1))
rowsInCh2 = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh2), axis=1))

fill_vals = {'<i4':-1, '<f8':np.nan, '<i8':-1}

for chConfig in ('Ch1', 'Ch2'):
    rows = eval('rowsIn' + chConfig)
    if len(rows[0]) > 0:
        outfilename = head+'-'+chConfig+tail
        outputfile = open(outfilename, 'wb')
        print "Creating File: " + outfilename
        columnsAll = eval('columnsIn' + chConfig)
        #create the fully reduced data table's dtype and column labels
        colToSplit = np.logical_and(columnsAll,columnsToSplit)
        colToKeep = np.logical_and(columnsAll,columnsToKeep)
        typeList = []
        redColLabels = []
        for i in xrange(0,len(data.dtype.descr)):
            typeTuple = data.dtype.descr[i]
            if colToKeep[i] or colToSplit[i]:
                typeList.append(typeTuple)
                redColLabels.append(columnLabels[i])
        datareduced = np.ma.masked_all(len(rows[0]), np.dtype(typeList))
        fill_values = [fill_vals.get(t[1],'') for t in typeList]
            
        #UNSTACK DATA
        for column in columnNames[colToKeep]:
            datareduced[column] = data[column][rows]
            datareduced.mask[column] = data.mask[column][rows]
        for column in columnNames[colToSplit]:
            pickedRows = rows[0]
            datareduced[column] = data[column][pickedRows]
            datareduced.mask[column] = data.mask[column][pickedRows]
        #fill empty cells
        datareduced.set_fill_value(fill_values)
        np.savetxt(outputfile, datareduced.filled(), fmt='%s', delimiter=',', header = header + ','.join(redColLabels), comments='')
        outputfile.close()