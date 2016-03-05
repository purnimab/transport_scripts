#!/usr/bin/python

#not written yet
#This script takes a Dynacool file or prefix, and combines all files with the same prefix of data HallA, HallB, VDPA, and VDPB

import sys
import string
import os.path
import numpy as np

ROWSPERIV = 1025 #256 per quadrant, 4 quadrants
channelNames = ['HallA', 'HallB', 'VDPA', 'VDPB']
resColNames = [('Resistance', 'Ohms'), ('Resistance Std Dev', 'Ohms'), ('Phase Angle', 'deg'), ('Gain',''), ('2nd Harmonic','dB'), ('3rd Harmonic','dB')]
ivColNames = [('IV Current', 'mA'), ('IV Voltage', 'V'), ('Gain', '')]
allColNames = [('Temperature', 'K'), ('Field', 'Oe'), ('Sample Position', 'deg')]

filename = sys.argv[1]
if os.path.exists(filename):
    #if the input name is actually a file, figure out the other filenames
    prefix = filename.rpartition('-')[0]
else:
    prefix = filename
filenames = []
for i in xrange(0,4):
    filenames.append(prefix+'-'+channelNames[i]+'.dat')

def read_data(fname):
    headerline = ''
    nline = 0
    inputfile = open(fname,'r')
    while headerline != "[Data]\n":
        headerline = inputfile.readline()
        nline += 1
    inputfile.close()
    #read in data file, masked. not all of them will be filled. the mask has values of true where data is missing
    return np.genfromtxt(fname, delimiter=',', skip_header=nline, filling_values=np.NaN, usecols=(2,3,4,6,7,8,9,10,23,24,25,26,27,28,29,30,43,44,45), names=True, usemask=True)

#data to put into the resistance or IV files, respectively - each row can only be one of these
dtyperes = np.dtype([('Temperature (K)', float), ('Field (Oe)', float), ('Sample Position (deg)', float)] + [(resColNames[i][0]+' '+channelNames[j]+(' ('+resColNames[i][1]+')' if resColNames[i][1] != '' else ''), float) for i in xrange(0,6) for j in xrange(0,4)] + [('Measured Channel', 'S5')])
resCh1Mask = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'Resistance_Ch1_Ohms', 'Resistance_Std_Dev_Ch1_Ohms', 'Phase_Angle_Ch1_deg', 'Gain_Ch1', '2nd_Harmonic_Ch1_dB', '3rd_Harmonic_Ch1_dB']
rch1ind = []
rch1 = None
#add Hall/VDP names to dtypeiv, define maskHallA etc.
dtypeiv = np.dtype([('Temperature (K)', float), ('Field (Oe)', float), ('Sample Position (deg)', float), ('IV Current (mA)', float), ('IV Voltage (V)', float), ('Gain', float)])
ivCh1Mask = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'IV_Current_Ch1_mA', 'IV_Voltage_Ch1_V', 'Gain_Ch1']
ivch1ind = []
ivch1 = None
resCh2Mask = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'Resistance_Ch2_Ohms', 'Resistance_Std_Dev_Ch2_Ohms', 'Phase_Angle_Ch2_deg', 'Gain_Ch2', '2nd_Harmonic_Ch2_dB', '3rd_Harmonic_Ch2_dB']
rch2ind = []
rch2 = None
ivCh2Mask = ['Temperature_K', 'Field_Oe', 'Sample_Position_deg', 'IV_Current_Ch2_mA', 'IV_Voltage_Ch2_V', 'Gain_Ch2']
ivch2ind = []
ivch2 = None

#initialize 

#this file will be all Ch1 or Ch2, but could have both IV curves and resistance measurements
for i in xrange(0, 4): #loop through each data file
    data = read_data(filenames[i])
    
    for j in xrange(0, data.shape[0]): #loop through each row of data
        #add the index of ___ind of whatever measurement type it is
        if not data.mask['Resistance_Ch1_Ohms'][0]: #res ch1 measurement exist
            rch1ind.append(j)
        elif not data.mask['Resistance_Ch2_Ohms'][0]: #res ch2 measurement exists
            rch2ind.append(j)
        elif not (data.mask['IV_Current_Ch1_mA'][0] or data.mask['IV_Voltage_Ch1_V'][0]): #IV ch1 measurement exists
            ivch1ind.append(j)
        elif not (data.mask['IV_Current_Ch2_mA'][0] or data.mask['IV_Voltage_Ch2_V'][0]): #IV ch2 measurement exists
            ivch2ind.append(j)
    #TODO: put the values in the right columns by changing the size of the array
    if len(rch1ind) > 0:
        if rch1 is not None:
            rch1 = np.append(rch1,data.data[resCh1Mask][rch1ind])
        else:
            rch1 = data.data[resCh1Mask][rch1ind]
        rch1ind = []
    if len(rch2ind) > 0:
        if rch2 is not None:
            rch2 = np.append(rch2,data.data[resCh2Mask][rch2ind])
        else:
            rch2 = data.data[resCh2Mask][rch2ind]
        rch2ind = []
    if len(ivch1ind) > 0:
        if ivch1 is not None:
            ivch1 = np.append(ivch1,data.data[ivCh1Mask][ivch1ind])
        else:
            ivch1 = data.data[ivCh1Mask][ivch1ind]
        ivch1ind = []
    if len(ivch2ind) > 0:
        if ivch2 is not None:
            ivch2 = np.append(ivch2,data.data[ivCh2Mask][ivch2ind])
        else:
            ivch2 = data.data[ivCh2Mask][ivch2ind]
        ivch2ind = []

headerEnd = " measurements concatenated from:\n"+'\n'.join(filenames)+'\n[Data]\n'
headerRes = "[Header]\nThis file contains Dynacool ETO Resistance" + headerEnd + "Temperature (K)\tField (Oe)\tSample Position (deg)\tResistance (Ohms)\tResistance Std Dev (Ohms)\tPhase Angle (deg)\tGain\t2nd Harmonic (dB)\t3rd Harmonic (dB)"
headerIV = "[Header]\nThis file contains Dynacool ETO IV curve" + headerEnd + "Temperature (K)\tField (Oe)\tSample Position (deg)\tIV Current (mA)\tIV Voltage (V)\tGain"

#print out files
outputs = [prefix+'-ResCh1.dat', prefix+'-IVCh1.dat', prefix+'-ResCh2.dat', prefix+'-IVCh2.dat']
if rch1 is not None:
    np.savetxt(outputs[0], rch1, header = headerRes, comments='')
if ivch1 is not None:
    np.savetxt(outputs[1], ivch1, header = headerIV, comments='')
if rch2 is not None:
    np.savetxt(outputs[2], rch2, header = headerRes, comments='')
if ivch2 is not None:
    np.savetxt(outputs[3], ivch2, header = headerIV, comments='')