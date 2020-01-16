import numpy as np
from scipy.optimize import curve_fit, minimize

def readDynResSimpdatfile(filename):
    return np.genfromtxt(filename, names=True, skip_header=1, delimiter=',')

def M(thetaM, phi, ms):
    return ms*np.array([np.sin(thetaM)*np.cos(phi),
                        np.sin(thetaM)*np.sin(phi),
                        np.cos(thetaM)])

def H(theta,h,beta):
    if beta:
        return h*np.array([0,-np.sin(theta),np.cos(theta)])
    else:
        return h*np.array([-np.sin(theta),0,np.cos(theta)])

def F(angsM, theta, h, beta, ms, meff, hcub, phi0):
    thetaM = angsM[0]
    phi = angsM[1]

    MM = M(thetaM,phi,ms)
    HH = H(theta,h,beta)

    return -np.dot(HH,MM) + ms/2*(meff*np.cos(thetaM)**2-hcub*np.sin(thetaM)**4*(3.+np.cos(4.*(phi+phi0)))/8.)

angsM = [(0,0)]
fit = []
ms = 320.*4*np.pi
meff = 298.
hcub = -20.

def calcAngsM(theta, h, beta, ms, meff, hcub, phi0):
    global angsM, fit
    if beta:
        angsM = [(theta[0],-np.pi/2.)]
    else:
        angsM = [(theta[0],np.pi)]
    fit = []
    for t in theta:
        fit.append(minimize(F, x0=angsM[-1], args=(t, h, beta, ms, meff, hcub, phi0)))
        angsM.append(fit[-1].x)
    angsM.pop(0)
    return angsM, fit

def ADMR(h, beta, rho0, drho1, drho2): #assumes you've already calculated the angles
    thetaM = [a[0] for a in angsM]
    phi = [a[1] for a in angsM]
    mm = M(thetaM, phi, 1.)
    return rho0 + drho1*mm[1]**2 + drho2*mm[0]**2

def ADMRshape(t,h,beta,rho0,drho1,drho2,ms,meff,hcub,phi0):
    angsM, f = calcAngsM(t,h,beta,ms,meff,hcub,phi0)
    thetaM = [a[0] for a in angsM]
    phi = [a[1] for a in angsM]
    mm = M(thetaM, phi, 1.)
    return rho0 + drho1*mm[1]**2 + drho2*mm[0]**2

def fitADMRshape(thetas, rhos, h, beta, phi0):
    guess = [rhos[0], max(rhos[0]-np.min(rhos), rhos[0]-np.max(rhos),key=abs),320.*4*np.pi,298.,-20.]
    if beta:
        fit = curve_fit(lambda x,r0,r1,ms,meff,hcub: ADMRshape(x,h,beta,r0,r1,0,ms,meff,hcub,phi0), thetas, rhos, guess)
    else:
        fit = curve_fit(lambda x,r0,r2,ms,meff,hcub: ADMRshape(x,h,beta,r0,0,r2,ms,meff,hcub,phi0), thetas, rhos, guess)
    return fit

def setMag(msl,meffl,hcubl):
    global ms, meff, hcub
    ms = msl
    meff = meffl
    hcub = hcubl
    return

def fitADMR(thetas, rhos, h, beta, phi0):
    calcAngsM(thetas, h, beta, ms, meff, hcub, phi0)
    guess = [rhos[0], max(rhos[0]-np.min(rhos), rhos[0]-np.max(rhos),key=abs)]
    if beta:
        fit = curve_fit(lambda x,r0,r1: ADMR(h,beta,r0,r1,0), thetas, rhos, guess)
    else:
        fit = curve_fit(lambda x,r0,r2: ADMR(h,beta,r0,0,r2), thetas, rhos, guess)
    return fit
