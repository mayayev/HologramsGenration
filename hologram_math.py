# -*- coding: utf-8 -*-
"""
Created on Sun May 17 21:23:40 2015

@author: M&M
"""
from __future__ import division
import numpy as np


def grid(n, ny=None, nz=None, triangular=False, theta=0, rad=False):
    
    ndims = 3
    if ny is None:
        ny=n
    if nz is None:
        nz=1
        ndims = 2
        
    npts = n*ny*nz
    
    xc = (n-1)/2
    yc = (ny-1)/2
    zc = (nz-1)/2
    
    x = np.arange(npts)
    y = np.floor(x/n)
    if nz > 1:
        y = np.mod(y, ny)
        z = np.floor(x/(n*ny))
    x = np.mod(x, n)
    
    # triangular lattice
    if triangular is not False:
        x = x - 0.25 + 0.5*np.mod(y,2)
        y = y * np.sqrt(3)/2
        yc = yc * np.sqrt(3)/2
    
    g = np.vstack([x-xc, y-yc])
    if nz > 1:
        g = np.vstack([g, z-zc])
    
    # rotation    
    if rad is True:
        theta = theta*np.pi/180
    g_rot = np.zeros([ndims, npts])    
    g_rot[0,] = g[0,] * np.cos(theta) + g[1,] * np.sin(theta)
    g_rot[1,] = -1*g[0,] * np.sin(theta) + g[1,] * np.cos(theta)

    return g_rot

    
def circle(n, r=1):
    theta = np.arange(0, 2*np.pi, 2*np.pi/n)
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    circ = np.vstack([x, y])
    
    return circ