# -*- coding: utf-8 -*-
import pylab as py
from matplotlib.widgets import Slider, Button#, CheckButtons, RadioButtons
import sys
import numpy as np
import os.path

#VARIABLES SETTING COLUMN POSITIONS
names = ['HallA', 'HallB', 'VDPA', 'VDPB']
plotpos = [1,3,2,4]

#LOAD FILE from an IVCurves.txt file
filename = sys.argv[1]
#check that it's an IVCurve
if os.path.split(filename)[1].count('IVCurves') == 0:
    print filename + " is not an IVCurves file"
    sys.exit(1)
filetypes = ['RvsT', 'RvsH']
#check file type - they have different orderings of columns
if filename.count('RvsT') > 0:
    filetype = 0
    print filetypes[filetype] + ' Sweep'
elif filename.count('RvsH') > 0:
    filetype = 1
    print filetypes[filetype] + ' Sweep'
elif len(sys.argv) > 1:
    filetype = int(sys.argv[2] == 'RvsH')
    print filetypes[filetype] + ' Sweep'
else:
    print filename + " is neither RvsH nor RvsT"
    sys.exit(1)

IVCurves = np.genfromtxt(filename, comments='\\*', delimiter='\t', names=('Temperature', 'Field', 'VHallA', 'IHallA', 'VHallB', 'IHallB', 'VVDPA', 'IVDPA', 'VVDPB', 'IVDPB')) #load text into a structured array, ignoring header lines
print "Reading " + filename

s = IVCurves.shape
#determine number of points per IV curve
pointsPerCurve = 1
if filetype==0:
    while pointsPerCurve < s[0] and IVCurves['Temperature'][pointsPerCurve] == IVCurves['Temperature'][0]: #check same temperature value
        pointsPerCurve += 1
else:
    while pointsPerCurve < s[0] and IVCurves[pointsPerCurve]['Field'] == IVCurves[0]['Field']: #check same field value
        pointsPerCurve += 1
print str(pointsPerCurve) + " points per IV curve"

#reshape IV curves so that the first index is the curve #, then each point in the curve, then temp, current, voltage, etc.
IVCurves.shape = (s[0]/pointsPerCurve, pointsPerCurve)
s = IVCurves.shape

#create inlier and outlier masks for each curve for RANSAC robust linear regression, to remove bad data points
#inlier_mask = np.array(

#save output file
def saveoutput(event):
    if os.path.exists(linregfilename):
        print "Overwriting " + linregfilename
    else:
        print "Saving " + linregfilename
    #columns = 'Temperature (K)\tField (Oe)\tRes ' + ' (Ohm)\tRes '.join(names) + ' (Ohm)\tErr ' + ' (Ohm)\tErr '.join(names) + ' (Ohm)\tDC Off ' + ' (V)\tDC Off '.join(names) + ' (V)'
    names.insert(0,'') #insert an empty string for easy naming of things
    columns = 'Temperature\tField' + '\tRes '.join(names) + '\tErr '.join(names) + '\tDC Off '.join(names) + '\nK\tOe' + '\tOhm'*8 + '\tV'*4
    names.pop(0) #remove the empty string
    np.savetxt(linregfilename, results, delimiter='\t', header=columns, comments='\\*')
    sys.stdout.flush()

#Get Linear regression, either calculated or from existing file
linregfilename = 'LinFitPy'.join(filename.rsplit('IVCurves',1)) # replaces the last occurrence of IVCurves with LinFitPy
if os.path.exists(linregfilename):
    results = np.genfromtxt(linregfilename, comments='\\*', delimiter='\t', names=['Temperature', 'Field'] + ['Res' + name for name in names] + ['Err' + name for name in names] + ['DCOff' + name for name in names]) #load text into structured array
    print "Reading " + linregfilename
else:
    #Calculate linear fit for each IV curve
    s = IVCurves.shape
    results = np.zeros(s[0], dtype=[('Temperature', '<f8'), ('Field', '<f8')] + [('Res'+name, '<f8') for name in names] + [('Err'+name, '<f8') for name in names] + [('DCOff'+name, '<f8') for name in names])
    results['Temperature'] = IVCurves['Temperature'][:,0] #T, H are the same for each curve, so arbitrarily choose the first value
    results['Field'] = IVCurves['Field'][:,0]
    
    for curve in xrange(0, s[0]):
        for i in xrange(0,4):
            p, cov = np.polyfit(IVCurves['I'+names[i]][curve,:], IVCurves['V'+names[i]][curve,:], 1, cov=True)
            
            #save resistance, error in resistance, and DC offset
            results['Res'+names[i]][curve] = p[0]
            results['DCOff'+names[i]][curve] = p[1]
            results['Err'+names[i]][curve] = cov[0,0]
    
#values for slider
tempindex = 0
minindex = 0
maxindex = s[0]-1

#create a figure
fig = py.figure(figsize=(20,15))
py.subplots_adjust(left=0.25, bottom=0.25)
py.hold(True)

#plot curves - called every time the IVcurve index is changed
def plotcurves(index):
    #pick out data for a specific temperature/field index
    curve = IVCurves[index]
    
    #plot each curve
    for i in xrange(0,4):
        py.subplot(2,2,plotpos[i])
        py.cla()
        #plot data
        x = curve['V'+names[i]]
        y = curve['I'+names[i]]
        py.plot(x, y,'.')
        #plot linear fit
        fitcurve = np.poly1d([results['Res'+names[i]][index], results['DCOff'+names[i]][index]])
        py.plot(fitcurve(y), y, 'r')
        py.title(names[i])
        py.xlabel('Voltage (V)')
        py.ylabel('Current (A)')
    
    py.suptitle(str(curve['Temperature'][0])+" K, "+str(curve['Field'][0])+" Oe", size='x-large', x = .55, y = .2)
    py.draw()
    #py.hold(False)
plotcurves(tempindex)

#create a slider
axindex = py.axes([0.25, 0.1, 0.65, 0.03], axisbg='lightgoldenrodyellow')
sliderindex = Slider(axindex, 'IV Curve #', minindex, maxindex, valinit=tempindex, dragging=True, valfmt='%1.0f')

#when the slider value changes, the index of the IV curve also changes
def update(val):
    tempindex = round(sliderindex.val)
    if sliderindex.val != tempindex: #otherwise, set_val calls update, and you get infinite recursion
        sliderindex.set_val(tempindex)
    plotcurves(tempindex)
sliderindex.on_changed(update)

#bind left/right arrow keys to also change the index of the IV curve
def leftright(event):
    if event.key == 'left':
        if sliderindex.val > minindex:
            tempindex = sliderindex.val-1
            sliderindex.set_val(tempindex)
    elif event.key == 'right':
        if sliderindex.val < maxindex:
            tempindex = sliderindex.val+1
            sliderindex.set_val(tempindex)
arrows = fig.canvas.mpl_connect('key_release_event', leftright)

#create a save fits button
axbutt = py.axes([0.025, 0.1, 0.1, 0.03])
button = Button(axbutt, 'Save Fits')
#button.on_clicked(saveoutput(names, linregfilename, results))
button.on_clicked(saveoutput)
py.show()