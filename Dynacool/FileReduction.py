#!/usr/bin/python

#This program takes a raw Dynacool resistance measurement .dat file and
#   1. separates non-empty columns into files for channel 1 and channel 2 (separate samples) and unstacks the VDP and Hall channels (separate measurements)
#   2. outputs a "concise" file containing only resistance, phase, and gain data
#   3. calculates resistance for values where resistance was not written, but voltage and current were

import sys
import string
import os.path
import numpy as np

#READ IN ARGUMENTS - [filename] [configuration #1 - default VDPandHall] [# pts - default 1] [configuration #2] [#pts #2]
#first argument is the filename to read in
filename = sys.argv[1]
head, tail = os.path.splitext(filename)
modes = {"VDPandHall":['VDPA', 'VDPB', 'HallA', 'HallB'] , "VDP":['VDPA', 'VDPB'], "Hall":['HallA', 'HallB'], "VDPA":['VDPA'], "VDPB":['VDPB'], "HallA":['HallA'], "HallB":['HallB'], "VDPBandHallB":['VDPB', 'HallB'], "4wire":['4wire']}
#if there are no more arguments, here are the default settings
modeCh1 = modes['VDPandHall']
modeCh2 = modeCh1
numPtsCh1 = 1
numPtsCh2 = numPtsCh1
#second argument is the unstacking mode for Ch1 (order points were taken in)
if len(sys.argv) == 3:
    modeCh1 = modes[sys.argv[2]]
    modeCh2 = modeCh1
#third argument is either the number of times each channel is repeated or the unstacking mode for Ch2
elif len(sys.argv) == 4:
    modeCh1 = modes[sys.argv[2]]
    try:
        numPtsCh1 = int(sys.argv[3])
        modeCh2 = modeCh1
        numPtsCh2 = numPtsCh1
    except ValueError as a:
        modeCh2 = modes[sys.argv[3]]
elif len(sys.argv) == 5:
    modeCh1 = modes[sys.argv[2]]
    try: #VDPA 2 VDPB would have ch2 have 1 pt per config
        numPtsCh1 = int(sys.argv[3])
        modeCh2 = modes[sys.argv[4]]
    except ValueError as a: #VDPA VDPB 2 would have chs 1 and 2 have 2 pts per 
        modeCh2 = modes[sys.argv[3]]
        numPtsCh1 = int(sys.argv[4])
        numPtsCh2 = numPtsCh1
elif len(sys.argv) > 5:
    modeCh1 = modes[sys.argv[2]]
    numPtsCh1 = int(sys.argv[3])
    modeCh2 = modes[sys.argv[4]]
    numPtsCh2 = int(sys.argv[5])
print "Reading File: " + filename
print "Unstacking Channel 1: " + ','.join(modeCh1) + ', ' + str(numPtsCh1) + " points each"
print "Unstacking Channel 2: " + ','.join(modeCh2) + ', ' + str(numPtsCh2) + " points each"

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
#find out which columns belong to Resistance/dV-dI or IV Curve measurements
columnsInRes = np.array([a.find('IV') == -1 for a in columnNames])
columnsIndVdI = columnsInRes
columnsInIV = np.array([a.find('Resistance') == -1 and a.find('Phase') == -1 and a.find('Harmonic') == -1 for a in columnNames])
#combine columns
columnsInCh1Res = np.logical_and(columnsInCh1, columnsInRes)
columnsInCh2Res = np.logical_and(columnsInCh2, columnsInRes)
columnsInCh1dVdI = columnsInCh1Res
columnsInCh2dVdI = columnsInCh2Res
columnsInCh1IV = np.logical_and(columnsInCh1, columnsInIV)
columnsInCh2IV = np.logical_and(columnsInCh2, columnsInIV)
#figure out which columns are unique to a specific measurement
colsUniquetoCh1Res = np.logical_and(columnsInCh1Res, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1IV),columnsInCh2IV)))
colsUniquetoCh2Res = np.logical_and(columnsInCh2Res, np.logical_not(np.logical_or(np.logical_or(columnsInCh1Res, columnsInCh1IV),columnsInCh2IV)))
colsUniquetoCh1dVdI = colsUniquetoCh1Res
colsUniquetoCh2dVdI = colsUniquetoCh2Res
colsUniquetoCh1IV = np.logical_and(columnsInCh1IV, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1Res),columnsInCh2IV)))
colsUniquetoCh2IV = np.logical_and(columnsInCh2IV, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1IV),columnsInCh1Res)))

#FIGURE OUT WHICH ROWS CORRESPOND TO WHICH COLUMNS
rowsInCh1Res = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh1Res), axis=1))
rowsInCh2Res = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh2Res), axis=1))
#dVdI has the same occupied data columns as Res, so the user must determine whether the data was res or dI/dV (unlike PPMS, T, H are taken at every point so will not be exactly the same)
rowsInCh1dVdI = rowsInCh1Res
rowsInCh2dVdI = rowsInCh2Res
rowsInCh1IV = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh1IV), axis=1))
rowsInCh2IV = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh2IV), axis=1))

fill_vals = {'<i4':-1, '<f8':np.nan, '<i8':-1}

#OUTPUT REDUCED FILES
ResSimpleColumns = ['Temperature_K', 'Field_Oe','Sample_Position_deg','Resistance_Ch1_Ohms', 'Phase_Angle_Ch1_deg', 'Gain_Ch1', 'Resistance_Ch2_Ohms', 'Phase_Angle_Ch2_deg', 'Gain_Ch2']
ResSimpleLabels = ['Temperature (K)', 'Field (Oe)', 'Sample Position (deg)', 'Resistance Ch1 (Ohms)', 'Phase Angle Ch1 (deg)', 'Gain Ch1', 'Resistance Ch2 (Ohms)', 'Phase Angle Ch2 (deg)', 'Gain Ch2']
dVdISimpleColumns = ['Temperature_K', 'Field_Oe','Sample_Position_deg','Resistance_Ch1_Ohms', 'Phase_Angle_Ch1_deg', 'AC_Current_Ch1_mA', 'DC_Current_Ch1_mA', 'Gain_Ch1', 'Resistance_Ch2_Ohms', 'Phase_Angle_Ch2_deg', 'AC_Current_Ch2_mA', 'DC_Current_Ch2_mA', 'Gain_Ch2']
dVdISimpleLabels = ['Temperature (K)', 'Field (Oe)', 'Sample Position (deg)', 'Resistance Ch1 (Ohms)', 'Phase Angle Ch1 (deg)', 'AC Current Ch1 (mA)', 'DC Current Ch1 (mA)', 'Gain Ch1', 'Resistance Ch2 (Ohms)', 'Phase Angle Ch2 (deg)', 'AC Current Ch2 (mA)', 'DC Current Ch2 (mA)', 'Gain Ch2']
IVSimpleColumns = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'IV_Current_Ch1_mA', 'IV_Voltage_Ch1_V', 'Gain_Ch1', 'IV_Current_Ch2_mA', 'IV_Voltage_Ch2_V', 'Gain_Ch2']
IVSimpleLabels = ['Temperature (K)', 'Field (Oe)', 'Sample Position (deg)', 'I-V Current Ch1 (mA)', 'I-V Voltage Ch1 (V)', 'Gain Ch1', 'I-V Current Ch2 (mA)', 'I-V Voltage Ch2 (V)', 'Gain Ch2']
for chConfig in ('Ch1', 'Ch2'):
    mode = eval('mode' + chConfig)
    numPts = eval('numPts' + chConfig)
    for measConfig in ('Res', 'dVdI', 'IV'):
        if measConfig == 'IV':
            ptsPerCurve = 1025 #256 per quadrant, 4 quadrants
        else:
            ptsPerCurve = 1
        rows = eval('rowsIn' + chConfig + measConfig)
        if len(rows[0]) > 0:
            outfilename = head+'-'+chConfig+measConfig+tail
            outputfile = open(outfilename, 'wb')
            print "Creating File: " + outfilename
            columnsAll = eval('columnsIn' + chConfig + measConfig)
            #create the fully reduced data table's dtype and column labels
            colToSplit = np.logical_and(columnsAll, columnsToSplit)
            colToKeep = np.logical_and(columnsAll,columnsToKeep)
            typeList = []
            redColLabels = []
            for i in xrange(0,len(data.dtype.descr)):
                typeTuple = data.dtype.descr[i]
                if colToKeep[i]:
                    typeList.append(typeTuple)
                    redColLabels.append(columnLabels[i])
                if colToSplit[i]:
                    for wireConfig in mode:
                        typeList.append((typeTuple[0].replace(chConfig, wireConfig), typeTuple[1]))
                        redColLabels.append(columnLabels[i].replace(chConfig, wireConfig))
            datareduced = np.ma.masked_all(len(rows[0]), np.dtype(typeList))
            fill_values = [fill_vals.get(t[1],'') for t in typeList]
                
            #UNSTACK DATA
            for column in columnNames[colToKeep]:
                datareduced[column] = data[column][rows]
                datareduced.mask[column] = data.mask[column][rows]
            for column in columnNames[colToSplit]:
                for i in xrange(0, len(mode)):
                    for j in xrange(0, numPts*ptsPerCurve):
                        start = i*numPts*ptsPerCurve+j
                        stride = len(mode)*numPts*ptsPerCurve
                        pickedRows = rows[0][start::stride]
                        newCol = column.replace(chConfig, mode[i])
                        datareduced[newCol][start::stride] = data[column][pickedRows]
                        datareduced.mask[newCol][start::stride] = data.mask[column][pickedRows]
            #fill empty cells
            datareduced.set_fill_value(fill_values)
            np.savetxt(outputfile, datareduced.filled(), fmt='%s', delimiter=',', header = header + ','.join(redColLabels), comments='')
            outputfile.close()
            
            outfilename = head + '-' + chConfig + measConfig + 'simp' + tail
            outputfile = open(outfilename, 'wb')
            print "Creating File: " + outfilename
            simpleHeader = 'Data obtained from: ' + os.path.basename(filename) + '\n'
            #create the simplified data table's dtype and column labels
            typeList = []
            simpColLabels = []
            simpCols = eval(measConfig+'SimpleColumns')
            simpLabels = eval(measConfig+'SimpleLabels')
            for i in xrange(0, len(simpCols)):
                column = simpCols[i]
                if column.find('Ch') == -1:
                    try:
                        typeList.append((column, datareduced[column].dtype.descr[0][1]))
                    except ValueError as a:
                        typeList.append((column, '<f8'))
                    simpColLabels.append(simpLabels[i])
                else:
                    for wireConfig in mode:
                        if column.find(chConfig) != -1:
                            newCol = column.replace(chConfig, wireConfig)
                            try:
                                typeList.append((newCol, datareduced[newCol].dtype.descr[0][1]))
                            except ValueError as a:
                                typeList.append((newCol, '<f8'))
                            simpColLabels.append(simpLabels[i].replace(chConfig, wireConfig))
            datasimplified = np.ma.masked_all(len(rows[0]), np.dtype(typeList))
            fill_values = [fill_vals.get(t[1],'') for t in typeList]
            datasimplified.set_fill_value(fill_values)
            #simplify data
            for colType in typeList:
                try:
                    datasimplified[colType[0]] = datareduced[colType[0]]
                    datasimplified.mask[colType[0]] = datareduced.mask[colType[0]]
                except ValueError as a:
                    datasimplified.mask[colType[0]].fill(True)
            #if resistance is blank, but current and voltage are available, calculate
            if measConfig == 'Res':
                count = 0
                for k in xrange(0, len(datareduced)):
                    for wireConfig in mode:
                        resCol = 'Resistance_'+wireConfig+'_Ohms'
                        currCol = 'AC_Current_'+wireConfig+'_mA'
                        voltCol = 'Voltage_Ampl_'+wireConfig+'_V'
                        if (datareduced.dtype.names.count(resCol) == 0 or datareduced[k].mask[resCol] or datareduced[k][resCol] == 0) and not datareduced[k].mask[currCol] and not datareduced[k].mask[voltCol]:
                            datasimplified[k][resCol] = datareduced[k][voltCol]/(1000*datareduced[k][currCol])
                            datasimplified[k].mask[resCol].fill(datasimplified[k].mask[resCol] == 0)
                            count = count + 1
                print str(count) + ' resistance points calculated'
            np.savetxt(outputfile, datasimplified.filled(), fmt='%s', delimiter=',', header = simpleHeader + ','.join(simpColLabels), comments='#')
            outputfile.close()
