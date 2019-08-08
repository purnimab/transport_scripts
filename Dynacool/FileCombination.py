#!/usr/bin/python

#This program takes a set of raw Dynacool resistance measurement .dat files and
#   1. clear empty columns of data
#   2. combines files by creating new columns with the VDP and Hall channels
#   3. calculates resistance for values where resistance was not written, but voltage and current were

import sys
import string
import os.path
import numpy as np

#READ IN ARGUMENTS
filename = sys.argv[1]
sorts = {'RvsT':'Temperature_K', 'RvsH':'Field_Oe', 'RvsPos':'Sample_Position_deg'}
if len(sys.argv) > 2:
    sortCol = sorts[sys.argv[2]]
else:
    sortCol = sorts['RvsH']
modes = {'VDPandHall':['VDPA', 'VDPB', 'HallA', 'HallB'], 'VDP':['VDPA', 'VDPB'], 'Hall':['HallA', 'HallB']}
if len(sys.argv) > 3:
    mode = modes[sys.argv[3]]
else:
    mode = modes['VDPandHall']
if os.path.exists(filename):
    #if the input name is actually a file, figure out the other filenames
    prefix = filename.rpartition('-')[0]
else:
    prefix = filename.strip('-')

#define output variables
dataConc = [[np.ma.masked_array([]), np.ma.masked_array([])], [np.ma.masked_array([]), np.ma.masked_array([])]]
fill_vals = {'<i4':-1, '<f8':np.nan, '<i8':-1}

#LOOP OVER INPUT FILES
for wireConfig in mode:
    fname = prefix + '-' + wireConfig + '.dat'
    if os.path.exists(fname):
        print "Reading File: " + fname
        inputfile = open(fname, 'rU')
        #SEPARATE EACH HEADER FROM DATA
        headerline = ''
        header = '' #only the header for the last read file will be kept
        nline = 0
        while headerline != "[Data]\n":
            headerline = inputfile.readline()
            nline += 1
            header = header + headerline
        columnLabels = inputfile.readline().strip().split(',') #list of columns, as formatted in the text file (with parentheses and spaces) - only the last will be kept
        inputfile.close()
        #READ IN DATA FILE, masked. Not all of them will be filled. The mask has values of True where data is missing
        data = np.genfromtxt(fname, delimiter=',', skip_header=nline, names=True, usemask=True, dtype=None)
        notmissing = np.logical_not(data.mask.view(bool).reshape(len(data),-1)) #turns the mask (an array of tuples) into a 2D array where True is a valid cell
        
        #FIGURE OUT WHICH COLUMNS TO WRITE OUT
        #find empty columns
        columnsWithData = np.logical_or.reduce(notmissing)
        #find out which columns belong to which channel
        columnNames = np.array(data.dtype.names)
        columnsInCh1 = np.array([a.find('Ch2') == -1 for a in columnNames])
        columnsInCh2 = np.array([a.find('Ch1') == -1 for a in columnNames])
        #find out which columns to unstack
        columnsToSplit = np.logical_xor(columnsInCh1, columnsInCh2)
        columnsToKeep = np.logical_and(columnsInCh1, columnsInCh2)
        #find out which columns belong to Resistance/dI-dV or IV Curve measurements
        columnsInRes = np.array([a.find('IV') == -1 for a in columnNames])
        columnsInIV = np.array([a.find('Resistance') == -1 and a.find('Phase') == -1 and a.find('Harmonic') == -1 for a in columnNames])
        #combine columns
        columnsInCh1Res = np.logical_and(columnsInCh1, columnsInRes)
        columnsInCh2Res = np.logical_and(columnsInCh2, columnsInRes)
        columnsInCh1IV = np.logical_and(columnsInCh1, columnsInIV)
        columnsInCh2IV = np.logical_and(columnsInCh2, columnsInIV)
        #figure out which columns are unique to a specific measurement
        colsUniquetoCh1Res = np.logical_and(columnsInCh1Res, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1IV),columnsInCh2IV)))
        colsUniquetoCh2Res = np.logical_and(columnsInCh2Res, np.logical_not(np.logical_or(np.logical_or(columnsInCh1Res, columnsInCh1IV),columnsInCh2IV)))
        colsUniquetoCh1IV = np.logical_and(columnsInCh1IV, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1Res),columnsInCh2IV)))
        colsUniquetoCh2IV = np.logical_and(columnsInCh2IV, np.logical_not(np.logical_or(np.logical_or(columnsInCh2Res, columnsInCh1IV),columnsInCh1Res)))
        
        #FIGURE OUT WHICH ROWS CORRESPOND TO WHICH COLUMNS
        rowsInCh1Res = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh1Res), axis=1))
        rowsInCh2Res = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh2Res), axis=1))
        rowsInCh1IV = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh1IV), axis=1))
        rowsInCh2IV = np.where(np.any(np.logical_and(notmissing, colsUniquetoCh2IV), axis=1))
        
        #ADD DATA TO RESPECTIVE CHANNELS
        for i in xrange(0,2):
            chConfig = ('Ch1', 'Ch2')[i]
            for j in xrange(0,2):
                measConfig = ('Res', 'IV')[j]
                rows = eval('rowsIn' + chConfig + measConfig)
                if len(rows[0]) > 0:
                    datareduced = dataConc[i][j]
                    start = datareduced.size
                    #determine dtype of concatenated data if this is the first file being read
                    if start == 0:
                        typeList = []
                        redColLabels = []
                        for k in xrange(0,len(data.dtype.descr)):
                            typeTuple = data.dtype.descr[k]
                            if typeTuple[1] == '|b1':
                                typeTuple = (typeTuple[0],'<f8')
                            if columnsToKeep[k]:
                                typeList.append(typeTuple)
                                redColLabels.append(columnLabels[k])
                            if not columnsInCh2[k]:
                                for wireConfig2 in mode:
                                    typeList.append((typeTuple[0].replace('Ch1', wireConfig2), typeTuple[1]))
                                    redColLabels.append(columnLabels[k].replace('Ch1', wireConfig2))
                        datareduced = np.ma.masked_all(len(rows[0]), np.dtype(typeList))
                        dataConc[i][j] = datareduced
                        fill_values = [fill_vals.get(t[1],'') for t in typeList]
                        #fill empty cells
                        datareduced.set_fill_value(fill_values)
                    else:
                        dataConc[i][j] = np.ma.resize(datareduced, start+len(rows[0]))
                        dataConc[i][j].mask[start:].fill(True)
                        dataConc[i][j].set_fill_value(fill_values)
                    
                    datareduced = dataConc[i][j]
                    #SEPARATE DATA
                    for column in columnNames[np.logical_and(columnsToKeep, columnsWithData)]:
                        try:
                            datareduced[column][start:] = data[column][rows]
                            datareduced.mask[column][start:] = data.mask[column][rows]
                        except ValueError as a:
                            datareduced.mask[column][start:].fill(True)
                    #TODO: fix how this deals with rows where there is data for both channels (not using the switchbox)
                    for column in columnNames[np.logical_and(columnsToSplit, columnsWithData)]:
                        newCol = column.replace(chConfig, wireConfig)
                        try:
                            datareduced[newCol][start:] = data[column][rows]
                            datareduced.mask[newCol][start:] = data.mask[column][rows]
                        except ValueError as a:
                            datareduced.mask[newCol][start:].fill(True)

#figure out which columns to write to a simplified file
ResSimpleColumns = ['Temperature_K', 'Field_Oe','Sample_Position_deg','Resistance_Ch1_Ohms', 'Phase_Angle_Ch1_deg', 'Gain_Ch1', 'Resistance_Ch2_Ohms', 'Phase_Angle_Ch2_deg', 'Gain_Ch2']
ResSimpleLabels = ['Temperature (K)', 'Field (Oe)', 'Sample Position (deg)', 'Resistance Ch1 (Ohms)', 'Phase Angle Ch1 (deg)', 'Gain Ch1', 'Resistance Ch2 (Ohms)', 'Phase Angle Ch2 (deg)', 'Gain Ch2']
IVSimpleColumns = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'IV_Current_Ch1_mA', 'IV_Voltage_Ch1_V', 'Gain_Ch1', 'IV_Current_Ch2_mA', 'IV_Voltage_Ch2_V', 'Gain_Ch2']
IVSimpleLabels = ['Temperature (K)', 'Field (Oe)', 'Sample Position (deg)', 'I-V Current Ch1 (mA)', 'I-V Voltage Ch1 (V)', 'Gain Ch1', 'I-V Current Ch2 (mA)', 'I-V Voltage Ch2 (V)', 'Gain Ch2']

#OUTPUT REDUCED FILES
for i in xrange(0,2):
    chConfig = ('Ch1', 'Ch2')[i]
    for j in xrange(0,2):
        measConfig = ('Res', 'IV')[j]
        datareduced = dataConc[i][j]
        if datareduced.size > 0:
            #find empty columns
            notmissing = np.logical_not(datareduced.mask.view(bool).reshape(len(datareduced),-1)) #turns the mask (an array of tuples) into a 2D array where True is a valid cell            
            columnsWithData = np.logical_or.reduce(notmissing)
            columnNames = np.array(datareduced.dtype.names)
            outfilename = prefix+'-'+chConfig+measConfig+'.dat'
            outputfile = open(outfilename, 'wb')
            print "Creating File: " + outfilename
            redLabels = list(np.array(redColLabels)[columnsWithData])
            dataout = np.sort(datareduced[columnNames[columnsWithData]].filled(), order=[sortCol, 'Time_Stamp_s'], kind='mergesort')
            np.savetxt(outputfile, dataout, fmt='%s', delimiter=',', header = header + ','.join(redLabels), comments='')
            outputfile.close()
            
            #SIMPLIFY FILE
          
            outfilename = prefix + '-'+chConfig + measConfig + 'simp' + '.dat'
            outputfile = open(outfilename, 'wb')
            print "Creating File: " + outfilename
            filenames = [os.path.basename(prefix)+'-'+wireConfig for wireConfig in mode]
            simpleHeader = 'Data obtained from: ' + ','.join(filenames) + '\n'
            #create the simplified data table's dtype and column labels
            typeList = []
            simpColLabels = []
            simpCols = eval(measConfig+'SimpleColumns')
            simpLabels = eval(measConfig+'SimpleLabels')
            for k in xrange(0, len(simpCols)):
                column = simpCols[k]
                if column.find('Ch') == -1:
                    try:
                        typeList.append((column, datareduced[column].dtype.descr[0][1]))
                    except ValueError as a:
                        typeList.append((column, '<f8'))
                    simpColLabels.append(simpLabels[k])
                else:
                    for wireConfig in mode:
                        if column.find(chConfig) != -1:
                            newCol = column.replace(chConfig, wireConfig)
                            typeList.append((newCol, datareduced[newCol].dtype.descr[0][1]))
                            simpColLabels.append(simpLabels[k].replace(chConfig, wireConfig))
            datasimplified = np.ma.masked_all(len(datareduced), np.dtype(typeList))
            fill_values = [fill_vals.get(t[1],'') for t in typeList]
            datasimplified.set_fill_value(fill_values)
            #simplify data
            for colType in typeList:
                try:
                    datasimplified[colType[0]] = datareduced[colType[0]]
                    datasimplified.mask[colType[0]] = datareduced.mask[colType[0]]
                except ValueError as a: #if resistance was completely blank
                    datasimplified.mask[colType[0]].fill(True)
            #if resistance is blank, but current and voltage are available, calculate
            if measConfig == 'Res':
                for wireConfig in mode:
                    count = 0
                    resCol = 'Resistance_'+wireConfig+'_Ohms'
                    currCol = 'AC_Current_'+wireConfig+'_mA'
                    voltCol = 'Voltage_Ampl_'+wireConfig+'_V'
                    for k in xrange(0, len(datareduced)):
                        if (datareduced[k].mask[resCol] or datareduced[k][resCol] == 0) and not datareduced[k].mask[currCol] and not datareduced[k].mask[voltCol]:
                            datasimplified[k][resCol] = datareduced[k][voltCol]/(1000*datareduced[k][currCol])
                            datasimplified[k].mask[resCol].fill(datasimplified[k].mask[resCol] == 0)
                            count = count + 1
                    if count > 0:
                        print str(count) + ' resistance points calculated for channel ' + wireConfig
            np.savetxt(outputfile, np.sort(datasimplified.filled(), order=sortCol, kind='mergesort')[::-1], fmt='%s', delimiter=',', header = simpleHeader + ','.join(simpColLabels), comments='#')
            outputfile.close()
