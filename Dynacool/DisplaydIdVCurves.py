# -*- coding: utf-8 -*-
import pylab as py
from matplotlib.widgets import Slider#, Button, RectangleSelector#, CheckButtons, RadioButtons
import sys
import numpy as np
import os.path
from scipy.integrate import trapz, simps

#VARIABLES SETTING COLUMN POSITIONS
names = ['VDPA', 'VDPB', 'HallA', 'HallB']
plotpos = [2,4,1,3]

#LOAD FILE from a dIdVsimp.dat file
filename = sys.argv[1]
#check that it's an IVCurve
if os.path.split(filename)[1].count('dVdIsimp') == 0:
    print filename + " is not a dIdVsimp file"
    sys.exit(1)
numPtsPerQuadrant = int(sys.argv[2])
numQuadrants = int(sys.argv[3])
print str(numPtsPerQuadrant) + ' points per quadrant * ' + str(numQuadrants) + ' quadrants'
filetypes = ['RvsT', 'RvsH']
#check file type - they have different orderings of columns
if filename.count('RvsT') > 0:
    filetype = 0
    print filetypes[filetype] + ' Sweep'
elif filename.count('RvsH') > 0:
    filetype = 1
    print filetypes[filetype] + ' Sweep'
elif len(sys.argv) > 4:
    filetype = int(sys.argv[-1] == 'RvsH') #check last item
    print filetypes[filetype] + ' Sweep'
else:
    print filename + " is neither RvsH nor RvsT"
    sys.exit(1)

#SEPARATE HEADER FROM DATA
fileheader = '' #header will be added to the beginning of any output files
inputfile = open(filename,'r')
headerline = inputfile.readline()
nline = 0
while headerline[0] == '#' and headerline.find('Temperature') == -1:
    fileheader = fileheader + headerline[1:]
    nline += 1
    headerline = inputfile.readline()
print fileheader
inputfile.close()

dIdVCurves = np.genfromtxt(filename, comments='//*', delimiter=',', skip_header = nline, names=True) #load text into a structured array, ignoring header lines but not ignoring the first one
print "Reading " + filename

s = dIdVCurves.shape
#determine number of points per IV curve
pointsPerCurve = numPtsPerQuadrant*numQuadrants+ 1 #e.g. 4 quadrants with 20 points each would be 81 points (including the first one)
print str(pointsPerCurve) + " points per IV curve"

#reshape IV curves so that the first index is the curve #, then each point in the curve, then temp, current, voltage, etc.
newCurves = np.copy(dIdVCurves)
newCurves.resize((s[0]/(pointsPerCurve*4), 4, pointsPerCurve))
dIdVCurves = newCurves
s = dIdVCurves.shape

#values for slider
tempindex = 0
minindex = 0
maxindex = s[0]-1

#create a figure
fig = py.figure(figsize=(20,15))
py.subplots_adjust(left=0.25, bottom=0.25)
py.hold(True)

#calculate IV curves
#def calculateIV(tempIndex, channelIndex):
#    curve = dIdVCurves[tempIndex, channelIndex]
#    name = names[channelIndex]
#    


#plot curves - called every time the IVcurve index is changed
def plotcurves(index):
    #pick out data for a specific temperature/field index
    #curve = dIdVCurves[index]
    
    #plot each curve
    for i in xrange(0,4):
        curve = dIdVCurves[index,i]
        #plotcurve(index, i)
        py.subplot(2,2,plotpos[i])
        py.cla()
        #plot data
        name = names[i]
        colors = ('b', 'g', 'r', 'm')
        matchv = 0
        matchi = 0
        for quadrant in xrange(numQuadrants):
            nPts = numPtsPerQuadrant
            start = numPtsPerQuadrant*quadrant + 1
            if quadrant == 0:
                nPts += 1
                start -= 1
                matchi = curve[start]['DC_Current_'+name+'_mA']
                if matchi != 0:
                    matchv = matchi * curve[start]['Resistance_'+name+'_Ohms']
            i = list(curve[start:start+nPts]['DC_Current_'+name+'_mA'])
            ac = list(curve[start:start+nPts]['AC_Current_'+name+'_mA'])
            #di = curve[start:start+nPts]['AC_Current_'+name+'_mA']
            dvdi = list(curve[start:start+nPts]['Resistance_'+name+'_Ohms'])
            phase = list(curve[start:start+nPts]['Phase_Angle_'+name+'_deg'])
            v = []
            ipoint = 0
            while ipoint < len(dvdi):
                if dvdi[ipoint] == 0: #remove R=0 points
                    i.pop(ipoint)
                    ac.pop(ipoint)
                    dvdi.pop(ipoint)
                    phase.pop(ipoint)
                else:
                    ipoint += 1
            if i[0] != matchi:
                i.insert(0, matchi)
                ac.insert(0,ac[0])
                #v.insert(0, matchv)
                dvdi.insert(0,dvdi[0])
                phase.insert(0,0)
                #matchv += (i[0] - matchi)*dvdi[0]
                #match = 1
            v.append(matchv)
            for point in xrange(2,len(dvdi)+1):
                #v.append(simps(dvdi[:point], x=i[:point],even='last')+matchv)
                v.append(trapz(dvdi[:point], x=i[:point])+matchv)
#            v = simps(dvdi, x=i)
            
            if quadrant == 3 and name == 'VDPB':
                print i
                print len(i)
                print dvdi
                print len(dvdi)
                print v
                print len(v)
                
            #plot points
            py.plot(v, i, '.'+colors[quadrant]+'-')
            #py.plot(v, phase)
            #for point in xrange(0,nPts): #plot slopes
            matchi = i[-1]
            matchv = v[-1]
#        x = curve['V'+name]
#        y = curve['I'+name]
        #plot used points in black
        #py.plot(x[lower:upper], y[lower:upper],'.k')
        #plot linear fit
        #fitcurve = np.poly1d([results['Res'+name][index], results['DCOff'+name][index]])

        py.title(name)
        py.xlabel('Voltage (V)')
        py.ylabel('Current (mA)')
    
    py.suptitle(str(curve['Temperature_K'][0])+" K, "+str(curve['Field_Oe'][0])+" Oe", size='x-large', x = .55, y = .2)
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

py.show()