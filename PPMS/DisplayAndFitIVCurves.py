# -*- coding: utf-8 -*-
import pylab as py
from matplotlib.widgets import Slider, Button, RectangleSelector#, CheckButtons, RadioButtons
import sys
import numpy as np
import os.path
from scipy import stats

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

#create a range identifier to store which range of the IV curve is "good"
values = np.zeros((2, s[0]), dtype = [(name, '<i4') for name in names])
for name in names:
    values[name][1].fill(pointsPerCurve)
#values = [[0]*s[0], [pointsPerCurve]*s[0]] #first column is lower limit, second column is upper limit

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

#linear fit throwing out outliers
#starting from ends
def find_best_polyfit(index, name):
    x = IVCurves['I'+name][index]
    y = IVCurves['V'+name][index]
    low, high = 0, pointsPerCurve
    p_cutoff = 0.01
    debounce = 2 #need two points to have seen a change in slope

    seen = 0
    for j in xrange(0,pointsPerCurve):
        slope, intercept, r, p, stderr = stats.linregress(x[:j+1], y[:j+1])
        if p < p_cutoff:
            seen += 1
        else:
            seen = max(0, seen-1) #the two points don't have to be in a row

        if seen == debounce:
            low = max(j - debounce, 0)
            break

    seen = 0
    for j in xrange(1, pointsPerCurve-low):
        slope, intercept, r, p, stderr = stats.linregress(x[-j-1:], y[-j-1:])
        if p < p_cutoff:
            seen += 1
        else:
            seen = max(0, seen-1)

        if seen == debounce:
            high = pointsPerCurve - (j - debounce)
            break

    values[name][0][index] = low
    values[name][1][index] = high
    #slope, intercept, r, p, stderr = stats.linregress(x[low:high], y[low:high])
    #return slope, intercept, stderr
    return stats.linregress(x[low:high], y[low:high])
    
#starting from middle
def fitresfindindex(index, name):
    i = IVCurves['I'+name][index]
    v = IVCurves['V'+name][index]
    lower = int(len(i)/2-1)
    upper = int(len(i)/2+2)
    r2 = 1
    up = False
    lowerFound = False
    upperFound = False
    slope, intercept, r2, p, err = stats.linregress(i[lower:upper], v[lower:upper])
    perc = 9
    r2limit = .98#(1+perc)*abs(r2)-perc
    while (not (upperFound and lowerFound) and abs(r2) > r2limit):
        if up and not upperFound:
            slope, intercept, r2, p, err = stats.linregress(i[lower:upper+1], v[lower:upper+1])
            if abs(r2) <= r2limit:
                upperFound = True
            elif upper == len(i) - 1:
                upperFound = True
                upper += 1
            else:
                upper += 1
            up = lowerFound
        elif not up and not lowerFound:
            slope, intercept, r2, p, err = stats.linregress(i[lower-1:upper], v[lower-1:upper])
            if abs(r2) <= r2limit:
                lowerFound = True
            elif lower == 1:
                lowerFound = True
                lower -= 1
            else:
                lower -= 1
            up = not upperFound
        #p, cov = np.polyfit(i[lower:upper], v[lower:upper], 1, cov=True)
    values[name][0][index] = lower
    values[name][1][index] = upper
    #return np.polyfit(i[lower:upper], v[lower:upper], 1, cov=True)
    return stats.linregress(i[lower:upper], v[lower:upper])
            
def linearfit(index, name, findIndex):
    if findIndex == 0:
        i = IVCurves['I'+name][index]
        v = IVCurves['V'+name][index]
        start = values[name][0][index]
        end = values[name][1][index]
        slope, intercept, r2, p, err = stats.linregress(i[start:end], v[start:end])
    elif findIndex == 1:
        slope, intercept, r2, p, err = find_best_polyfit(index, name)
    elif findIndex == 2:
        slope, intercept, r2, p, err = fitresfindindex(index, name)
        
    results['Res'+name][index] = slope #p[0]
    results['DCOff'+name][index] = intercept #p[1]
    results['Err'+name][index] = err #cov[0,0]

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
            linearfit(curve, names[i], 1) #1 means to iterate through the fit to find outliers from outside, 2 from inside
            #p, cov = fitresfindindex(curve, names[i])
            #p, cov = np.polyfit(IVCurves['I'+names[i]][curve], IVCurves['V'+names[i]][curve], 1, cov=True)

            #save resistance, error in resistance, and DC offset
            #results['Res'+names[i]][curve] = p[0]
            #results['DCOff'+names[i]][curve] = p[1]
            #results['Err'+names[i]][curve] = cov[0,0]
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
        #plotcurve(index, i)
        py.subplot(2,2,plotpos[i])
        py.cla()
        #plot data
        name = names[i]
        x = curve['V'+name]
        y = curve['I'+name]
        #plot unused points in green
        lower = values[name][0][index]
        if lower > 0:
            py.plot(x[:lower], y[:lower], '.g')
        upper = values[name][1][index]
        if upper < pointsPerCurve:
            py.plot(x[upper:], y[upper:], '.g')
        #plot used points in blue
        py.plot(x[lower:upper], y[lower:upper],'.')
        #plot linear fit
        fitcurve = np.poly1d([results['Res'+name][index], results['DCOff'+name][index]])
        #plot fit in red
        py.plot(fitcurve(y), y, 'r')
        py.title(name)
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
    plotcurves(int(tempindex))
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

#bind rectangle selector to select points in IV curve to use for fitting
def line_select_callback0(eclick, erelease):
    line_select_callback(eclick, erelease, 0)
def line_select_callback1(eclick, erelease):
    line_select_callback(eclick, erelease, 1)
def line_select_callback2(eclick, erelease):
    line_select_callback(eclick, erelease, 2)
def line_select_callback3(eclick, erelease):
    line_select_callback(eclick, erelease, 3)
def line_select_callback(eclick, erelease, plotIndex):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    curveIndex = int(sliderindex.val)
    name = names[plotIndex]
    x = IVCurves['V'+name][curveIndex]
    y = IVCurves['I'+name][curveIndex]
    start = 0
    end = pointsPerCurve
    #determine which points in the IV curve are within the rectangle
    if y[0] < y[1]: #increasing current
        while y[start] < min(y1,y2) and start < pointsPerCurve-1:
            start += 1
        while y[end-1] > max(y1,y2) and end > 0:
            end -= 1
    else: #decreasing current
        while y[start] > max(y1,y2) and start < pointsPerCurve-1:
            start += 1
        while y[end-1] < min(y1,y2) and end > 0:
            end -= 1
    if x[0] < x[1]: #increasing voltage
        while x[start] < min(x1,x2) and start < pointsPerCurve-1:
            start += 1
        while x[end-1] > max(x1,x2) and end > 0:
            end -= 1
    else: #decreasing voltage
        while x[start] > max(x1,x2) and start < pointsPerCurve-1:
            start += 1
        while x[end-1] < min(x1,x2) and end > 0:
            end -= 1
    #update the fit and plot
    values[name][0][curveIndex] = start
    values[name][1][curveIndex] = end
    linearfit(curveIndex, name, 0) #0 means use precalculated edges of the dataset
    plotcurves(curveIndex)

#instantiate a separate rectangle selector for each plot
RS0 = RectangleSelector(py.subplot(2,2,plotpos[0]), line_select_callback0, drawtype='box', useblit=True, interactive=True)
RS1 = RectangleSelector(py.subplot(2,2,plotpos[1]), line_select_callback1, drawtype='box', useblit=True, interactive=True)
RS2 = RectangleSelector(py.subplot(2,2,plotpos[2]), line_select_callback2, drawtype='box', useblit=True, interactive=True)
RS3 = RectangleSelector(py.subplot(2,2,plotpos[3]), line_select_callback3, drawtype='box', useblit=True, interactive=True)
