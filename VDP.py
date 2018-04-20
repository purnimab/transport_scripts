import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

#read from a nanovoltmeter linear fit file (order of columns: temp, field, HallA, HallB, VDPA, VDPB)
def readNVLinFitResfile(filename):
    RvsTdata = np.loadtxt(filename, comments='\*')
    temperature = RvsTdata[:,0]
    field = RvsTdata[:,1]
    R_HallA = RvsTdata[:,2]
    R_HallB = RvsTdata[:,3]
    R_VDPA = RvsTdata[:,4]
    R_VDPB = RvsTdata[:,5]
    return temperature, field, R_VDPA, R_VDPB, R_HallA, R_HallB

#read from SMU RvsT linear fit file (order of columns: temp, Field, Max Current Ch1, Max Current Ch2, VDPA, VDPB, HallA, HallB)
def readSMULinFitResfile(filename):
    RvsTdata = np.loadtxt(filename, comments='\*')
    temperature = RvsTdata[:,0]
    field = RvsTdata[:,1]
    R_HallA = RvsTdata[:,6]
    R_HallB = RvsTdata[:,7]
    R_VDPA = RvsTdata[:,4]
    R_VDPB = RvsTdata[:,5]
    return temperature, field, R_VDPA, R_VDPB, R_HallA, R_HallB

#interpolate room temperature property
def roomtemp(temperature,RS):
    RvsT = interp1d(temperature, RS)
    RS_RT = RvsT(298.15)
    return RS_RT

##From Van der Pauw original paper

#function to minimize for VDP solution
def vdpfunc(rhod, RA, RB):
    return np.exp(-np.pi*RA/rhod) + np.exp(-np.pi*RB/rhod) - 1.

#solve for VDP sheet resistance
def vdp(R_VDPA, R_VDPB):
    guess = (R_VDPA+R_VDPB)/2*np.pi/np.log(2)
    RS = fsolve(vdpfunc, guess, (R_VDPA, R_VDPB), xtol=1.49012e-08)
    return RS

#alternate function to minimize for VDP solution
def vdpfunc2(e, RA, RB):
    return np.cosh(e*(RA/RB-1.)/(RA/RB+1.)) - np.exp(e)/2.

#alternate solution for VDP sheet resistance
def vdp2(R_VDPA, R_VDPB):
    guess2 = np.ones_like(R_VDPA)
    f = np.log(2)/fsolve(vdpfunc2, guess2, (R_VDPA, R_VDPB), xtol=1.49012e-08)
    RS2 = np.pi/np.log(2)*f*(R_VDPA+R_VDPB)/2.
    return RS2

#errors for a circular sample from the VDP paper
def vdpErrors(length, width, bond_distance):
    size = min(width,length)
    error_drho = -bond_distance**2/(2.*size**2*np.log(2.))*4.
    error_dRH = -2.*bond_distance/(np.pi*size)*4.
    return error_drho, error_dRH


## from S.H.N. Lim et al., Rev. Sci. Inst. 80, 075109 (2009) - "van der Pauw method for measuring resistivity of a plane sample with distant boundaries"

#function to minimize for VDP solution in the infinite plane limit
def vdpfarfunc(rhod, RA, RB):
    return np.exp(-2*np.pi*RA/rhod) + np.exp(-2*np.pi*RB/rhod) - 1.

#solve for VDP sheet resistance in the infinite plane limit
def vdpfar(R_VDPA, R_VDPB):
    guess = (R_VDPA+R_VDPB)*np.pi/np.log(2)
    RS = fsolve(vdpfarfunc, guess, (R_VDPA, R_VDPB), xtol=1.49012e-08)
    return RS