# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:12:28 2015

@author: user
"""

from __future__ import division
import scipy as sp
import numpy as np
from holograms.ApplyCalibration import loadcal
import holograms as holo


def maketheta(h, w, xc, yc):
        
    x = np.arange(w)
    y = np.arange(h)
    xmat, ymat = np.meshgrid(x - w/2 - xc, y - h/2 - yc)
    
    theta = sp.arctan2(ymat, xmat)
    
    return theta
    
    
def vortex(ell, dim=None, wavelength=None, xc=None, yc=None):
    
    if dim is None:
        w=holo.w
        h=holo.h
    else:
        w=dim[0]
        h=dim[1]
    
    if xc is None and yc is None:
        tempcal = None
    else:
        tempcal = {'xc': xc, 'yc': yc}
    
    cal = loadcal(cal=tempcal)
    
    theta = maketheta(w, h, cal['xc'], cal['yc'])
    phi = np.mod(theta*ell, 2*sp.pi)
    
    return phi
    
    
