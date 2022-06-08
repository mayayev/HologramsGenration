"""
Created on Thu Apr 16 15:14:14 2015

@author: Maya Yevnin
"""
from __future__ import division

def loadcal(cal=None):
    from scipy.io import loadmat

    if cal is None:
        mat = loadmat('E:\SLM calibration\cal.mat')
#        xc = mat['xc']
#        yc = mat['yc']
#        xfac = mat['xfac']
#        rfac_0 = mat['rfac_0']
#        rfac_slope = mat['rfac_slope']
#        thetafac = mat['thetafac']
#        x0_0 = mat['x0_0']
#        x0_slope = mat['x0_slope']
#        y0_0 = mat['y0_0']
#        y0_slope = mat['y0_slope']
        newcal = mat

    elif cal=='default' or cal=='Default':
        xc = 0          # Center of phase mask on SLM
        yc = 0          # Center of phase mask on SLM
        xfac = 1        # Scale factor for square pixels
        rfac_0 = 1      # Projection scale factor at z=0
        rfac_slope = 0  # Projection scale factor z correction
        thetafac = 0    # Orientation
        x0_0 = 300      # center of CCD image at z=0
        x0_slope = 0	 # center of CCD image z correction
        y0_0 = 400      # center of CCD image at z=0
        y0_slope = 0    # center of CCD image z correction
        newcal = {'xc': xc, 'yc': yc, 'xfac': xfac, 'rfac_0': rfac_0, 'rfac_slope': rfac_slope, 'thetafac': thetafac, 'x0_0': x0_0, 'x0_slope': x0_slope, 'y0_0': y0_0, 'y0_slope': y0_slope}
        
    else:
        newcal=cal
        
    return newcal


def applycal(points, cal, centered, screen_res):

    import numpy as np

    [ndims, npts] = points.shape
    corrected_points = np.ones((ndims, npts))
    thetafac = cal['thetafac']
    # apply rotation
    corrected_points[0,:] = points[0,:] * np.cos(thetafac) + points[1,:] * np.sin(thetafac)
    corrected_points[1,:] = points[1,:] * np.cos(thetafac) - points[0,:] * np.sin(thetafac)


    # center screen - usually used only for online apps such as point n click
    if centered:
        delta_x = screen_res[0]/2 - cal['x0_0']
        delta_y = screen_res[1]/2 - cal['y0_0']
        corrected_points[0,:] -= delta_x
        corrected_points[1,:] -= delta_y

    # z corrections
    if ndims==3:
        corrected_points[2,:] = points[2,:].copy()
        rfac = cal['rfac_0'] + cal['rfac_slope'] * corrected_points[2,:]
    else:
        rfac = cal['rfac_0']
    corrected_points[0,:] = cal['xfac'] * rfac * corrected_points[0,:]
    corrected_points[1,:] = rfac * corrected_points[1,:]

    return corrected_points



