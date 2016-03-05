#!/usr/bin/python

#not written yet
#This script takes a Dynacool .dat file and splits it into files only containing channel 1 or channel 2

import sys
import string

ROWSPERIV = 1025 #256 per quadrant, 4 quadrants

filename = sys.argv[1]
modes = {"VDPandHall":1, "VDP":2, "VDPA":3, "VDPB":4, "HallA":5, "HallB":6, "VDPBandHallB":7}
if len(sys.argv) == 2:
    modech1 = modech2 = 1
    print "Splitting " + filename + " into VDPandHall\n"
elif len(sys.argv) == 3:
    modech1 = modech2 = modes.get(sys.argv[2],1)
    print "Splitting " + filename + " into " + sys.argv[2] + '\n'
else:
    modech1 = modes.get(sys.argv[2],1)
    modech2 = modes.get(sys.argv[3],1)
    print "Splitting " + filename + " into: Ch1: " + sys.argv[2] + ", Ch2: " + sys.argv[3] + '\n'

inputfile = open(filename,'r')
outputfile = open(filename + "Split",'w')