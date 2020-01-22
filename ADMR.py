import numpy as np
from scipy.optimize import curve_fit, differential_evolution
#from scipy.optimize import minimize, basinhopping

def readDynResSimpdatfile(filename):
    return np.genfromtxt(filename, names=True, skip_header=1, delimiter=',')

def M(thetaM, phi, ms): #magnetization for given angles
    return ms*np.array([np.sin(thetaM)*np.cos(phi),
                        np.sin(thetaM)*np.sin(phi),
                        np.cos(thetaM)])

def H(theta,h,beta): #field for given angles
    if beta:
        return h*np.array([0,-np.sin(theta),np.cos(theta)])
    else:
        return h*np.array([-np.sin(theta),0,np.cos(theta)])

def F(angsM, theta, h, beta, ms, meff, hcub, phi0): #free energy function to be minimized
    #thetaM = angsM[0]
    #phi = angsM[1]
    thetaM, phi = angsM #unpacks the angle tuple

    MM = M(thetaM,phi,ms)
    HH = H(theta,h,beta)

    return -np.dot(HH,MM) + ms/2*(meff*np.cos(thetaM)**2 - hcub*np.sin(thetaM)**4*(3.+np.cos(4.*(phi+phi0)))/8.)

angsM = [(0,0)]
fit = []
ms = 320.*4*np.pi
meff = 298.*4*np.pi
hcub = -20.

def setMag(msl,meffl,hcubl):
    global ms, meff, hcub
    ms = msl
    meff = meffl
    hcub = hcubl
    return

def calcAngsM(theta, h, beta, ms, meff, hcub, phi0):
    global angsM, fit
    if beta:
        angsM = [(theta[0],-np.pi/2.)]
    else:
        angsM = [(theta[0],np.pi)]
    fit = []
    for t in theta:
        #fit.append(minimize(F, x0=angsM[-1], args=(t, h, beta, ms, meff, hcub, phi0)))
        fit.append(differential_evolution(F, bounds=[(0,np.pi),(0,2*np.pi)], args=(t, h, beta, ms, meff, hcub, phi0)))
        #fit.append(basinhopping(lambda x: F(x, t, h, beta, ms, meff, hcub, phi0), x0=angsM[-1]))#, args=(t, h, beta, ms, meff, hcub, phi0)))
        angsM.append(fit[-1].x)
    angsM.pop(0)
    return angsM, fit

def ADMR(h, beta, rho0, drho1, drho2): #assumes you've already calculated the angles
    thetaM = [a[0] for a in angsM]
    phi = [a[1] for a in angsM]
    mm = M(thetaM, phi, 1.)
    return rho0 + drho1*mm[1]**2 + drho2*mm[0]**2

def ADMRoff(t,h,beta,rho0,drho1,drho2,toff,phi0): #calculates the angles each time
    global ms, meff, hcub
    angsM, f = calcAngsM(t+toff,h,beta,ms,meff,hcub,phi0)
    thetaM = [a[0] for a in angsM]
    phi = [a[1] for a in angsM]
    mm = M(thetaM, phi, 1.)
    return rho0 + drho1*mm[1]**2 + drho2*mm[0]**2

def ADMRshape(t,h,beta,rho0,drho1,drho2,ms,meff,hcub,phi0): #calculates the angles each time
    angsM, f = calcAngsM(t,h,beta,ms,meff,hcub,phi0)
    thetaM = [a[0] for a in angsM]
    phi = [a[1] for a in angsM]
    mm = M(thetaM, phi, 1.)
    return rho0 + drho1*mm[1]**2 + drho2*mm[0]**2

def fitADMRshape(thetas, rhos, r0, dr, h, beta, phi0):
    global meff, hcub
    guess = [r0, dr, meff, hcub]
    if beta:
        fit = curve_fit(lambda x,r0,r1,meff,hcub: ADMRshape(x,h,beta,r0,r1,0,ms,meff,hcub,phi0), thetas, rhos, guess)
    else:
        fit = curve_fit(lambda x,r0,r2,meff,hcub: ADMRshape(x,h,beta,r0,0,r2,ms,meff,hcub,phi0), thetas, rhos, guess)
    return fit

def fitADMRoff(thetas, rhos, h, beta, phi0):
    switch = np.where(np.sign(thetas[1:-1]-thetas[:-2]) != np.sign(thetas[2:]-thetas[1:-1]))[0]+1 #endpoint where the sweep direction switches
    slices = [slice(None,switch[0]+1)]+[slice(x+1,switch[i+1]+1) for i,x in enumerate(switch[:-1])]+[slice(switch[-1]+1,None)] #separate segments of the sweep to separate output
    guess = [rhos[0], max(rhos[0]-np.min(rhos), rhos[0]-np.max(rhos),key=abs)]+[0.]*len(slices)
    if beta:
        fit = curve_fit(lambda x,r0,r1,*off: np.concatenate([ADMRoff(x[s],h,beta,r0,r1,0,off[i],phi0) for i,s in enumerate(slices)]), thetas, rhos, guess)
    else:
        fit = curve_fit(lambda x,r0,r2,*off: np.concatenate([ADMRoff(x[s],h,beta,r0,0,r2,off[i],phi0) for i,s in enumerate(slices)]), thetas, rhos, guess)
    return fit

def fitADMR(thetas, rhos, h, beta, phi0):
    calcAngsM(thetas, h, beta, ms, meff, hcub, phi0)
    guess = [rhos[0], max(rhos[0]-np.min(rhos), rhos[0]-np.max(rhos),key=abs)]
    if beta:
        fit = curve_fit(lambda x,r0,r1: ADMR(h,beta,r0,r1,0), thetas, rhos, guess)
    else:
        fit = curve_fit(lambda x,r0,r2: ADMR(h,beta,r0,0,r2), thetas, rhos, guess)
    return fit
