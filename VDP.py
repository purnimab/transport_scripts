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

#read from a Dynacool nanovoltmeter linear fit file (order of columns: temp, field, SheetRes)
def readDynNVSheetResfile(filename):
    RvsTdata = np.loadtxt(filename, comments='\*')
    temperature = RvsTdata[:,0]
    field = RvsTdata[:,1]
    SheetRes = RvsTdata[:,2]
    return temperature, field, SheetRes

#read from a Dynacool nanovoltmeter linear fit file (order of columns: temp, field, HallRes)
def readDynNVHallResfile(filename):
    RvsTdata = np.loadtxt(filename, comments='\*')
    temperature = RvsTdata[:,0]
    field = RvsTdata[:,1]
    HallRes = RvsTdata[:,2]
    return temperature, field, HallRes

#read from Dynacool reduced file, no interpolation - but returns masked arrays
def readDynResSimpdatfileNoInterp(filename):
    RvsTdata = np.genfromtxt(filename, delimiter=',', skip_header=1, names=True, usemask=True, dtype=None)
    temperature = RvsTdata['Temperature_K']
    field = RvsTdata['Field_Oe']
    if 'Resistance_VDPA_Ohms' in RvsTdata.dtype.names:
        R_VDPA = RvsTdata['Resistance_VDPA_Ohms']
        R_VDPA.mask = np.logical_or(R_VDPA.mask,np.isnan(R_VDPA))
        R_VDPB = RvsTdata['Resistance_VDPB_Ohms']
        R_VDPB.mask = np.logical_or(R_VDPB.mask,np.isnan(R_VDPB))
    else:
        R_VDPA = np.ma.array(np.zeros_like(temperature), mask=np.ones_like(temperature))
        R_VDPB = np.ma.array(np.zeros_like(temperature), mask=np.ones_like(temperature))
    if 'Resistance_HallA_Ohms' in RvsTdata.dtype.names:
        R_HallA = RvsTdata['Resistance_HallA_Ohms']
        R_HallA.mask = np.logical_or(R_HallA.mask,np.isnan(R_HallA))
        R_HallB = RvsTdata['Resistance_HallB_Ohms']
        R_HallB.mask = np.logical_or(R_HallB.mask,np.isnan(R_HallB))
    else:
        R_HallA = np.ma.array(np.zeros_like(temperature), mask=np.ones_like(temperature))
        R_HallB = np.ma.array(np.zeros_like(temperature), mask=np.ones_like(temperature))
    return temperature, field, R_VDPA, R_VDPB, R_HallA, R_HallB

#read from Dynacool reduced file, all values in a masked array
def readDynResSimpdatfileFull(filename):
    RvsTdata = np.genfromtxt(filename, delimiter=',', skip_header=1, names=True, usemask=True, dtype=None)
    for column in RvsTdata.dtype.names:
        data = RvsTdata[column]
        data.mask = np.logical_or(data.mask, np.isnan(data))
        if 'Gain' in column:
            data.mask = np.logical_or(data.mask, np.equal(data,-1))
    return RvsTdata

#interpolates resistances taken at slightly different temperatures or fields using a masked array
def interpolateSwitchboxScan(x, R_list):
    return (interpSwitchbox(x,r) for r in R_list)
def interpSwitchbox(x,r):
    mask = np.logical_and(np.logical_not(r.mask),np.logical_not(x.mask))
    try:
        rvsx = interp1d(x.data[mask],r.data[mask],bounds_error=False)
        rfull = rvsx(x.data)
    except ValueError as err:
        rfull = np.zeros_like(x.data)
        print err
    return rfull

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