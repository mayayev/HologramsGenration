# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from holograms.ApplyCalibration import loadcal, applycal

def fastphase(points, wavelength=None, dim=None, cal=None, screen=[640,480], centered=False):
    """
    Calculates the phase hologram encoding a desired trapping
    pattern by superposition as fast as possible.
    
    CALLING SEQUENCE:
    phi = fastphase(points)
    
    INPUTS:
          points: [2,npts] or [3,npts] array of (x,y) coordinates of points in the plane,
             relative to the center of the field of view.
             Coordinates are measured in arbitrary pixel units, unless
             the CAL keyword is set.
    
    KEYWORD PARAMETERS:
    	wavelength: 	of the trapping laser. Default according to computer
    	dim: 			[nx] or [nx,ny] array giving dimensions of DOE, in pixels.
    						Default: 600 x 800         
    	cal: 			Spatial calibration factors returned by CALIBRATE.
    						Default: square pixels, no rotation, arbitrary scale.
    	screen:			camera resolution
    						Default: 640x480
    	centered:		Is the hologram centered relative to the camera. Default is False
    	
               
    OUTPUTS:
          phi: phase pattern encoding pattern of traps described by
               points with values ranging from [0,2 pi].
    
    RESTRICTIONS:
          Can be very memory intensive for large numbers of points.
    
    PROCEDURE:
          Initial estimate is obtained by superposing the fields of the
          specified beams.
    
     NOTES:
          If we calibrate the SLM's phase transfer function, then we
          should be able to pass the lookup table to this function, and
          return a hologram of indices into the look-up table.
    
    MODIFICATION HISTORY:
     Created by David G. Grier, New York University, 12/12/2004.
     1/19/2005: DGG.  Implemented 3D.
     1/29/2005: DGG.  Major code clean-up leading to improved speed.
                      Implemented SUBSAMPLE for major speedup.
     11/5/2011: Maya Yevnin -  added wavelength parameterand calibration files
     Mar 2013: Maya Yevnin -   Moved calibration parameters from the xmat, ymat
                               vectors to new vector corrected_points, and
                               added zfac parameter
     May 2015: Maya Yevnin, translated to python
    """

    [ndims,npts] = points.shape
    
    # calibration
    cal = loadcal(cal=cal)
    corrected_points = applycal(points, cal, centered, screen)
    
    # SLM coordinates
    if dim is None:
        from holograms import h,w
    else:
        w = dim[0]
        h = dim[1]

    x = np.arange(w)
    y = np.arange(h)
    xmat, ymat = np.meshgrid(x - w/2 - cal['xc'], y - h/2 - cal['yc'])        

    # wavevectors associated with trap position (times i)
    ikx = 1j * (2*np.pi/w) * corrected_points[0,]
    iky = 1j * (2*np.pi/h) * corrected_points[1,]
    
    ikz = np.zeros(npts)
    rsq = np.zeros((h,w))
    if ndims > 2:
        aperture = 5000                            # 5mm aperture
        if wavelength is None:
            from Holograms import wavelength
        wavelength = (wavelength/1000)*w/aperture; # wavelength in aperture pixels
        na = 1.4                                   # numerical aperture of objective
        f = aperture/na                            # focal length in pixels
        ikz = 1j * (2*np.pi) * corrected_points[2,]/(wavelength*f^2)
        rsq = xmat**2 + ymat**2
    
    iphase = 1j * 2*np.pi * np.random.rand(npts)
    psi = np.zeros((h,w))
    
    for n in xrange(0, npts):
        ikxx = ikx[n] * x + iphase[n]
        ikyy = iky[n] * y
        ikzz = ikz[n] * rsq
        ex = np.exp(-ikxx).reshape([1,w])
        ey = np.exp(-ikyy).reshape([1,h])
        ez = np.exp(-ikzz)
        
        psi = psi + ex*ey.transpose() * ez 
    phi = np.angle(psi) + np.pi
    return phi



def phasemask(points, wavelength=None, dim=None, cal=None, screen=[640,480], centered=False, ell=None, win=None, alpha=1):
    """
    Calculates the phase hologram encoding a desired trapping
    pattern by superposition including vortices, windowing effect and trap intensities.
    This is the first stage of dsphase.
    
    CALLING SEQUENCE:
    phi = fastphase(points)
    
    INPUTS:
          points: [2,npts] or [3,npts] array of (x,y) coordinates of points in the plane,
             relative to the center of the field of view.
             Coordinates are measured in arbitrary pixel units, unless
             the CAL keyword is set.
    
    KEYWORD PARAMETERS:
    wavelength:  of the trapping laser. Default according to computer
    dim: 		[nx] or [nx,ny] array giving dimensions of DOE, in pixels.
    				Default: 600 x 800         
    cal: 		Spatial calibration factors returned by CALIBRATE.
    				Default: square pixels, no rotation, arbitrary scale.
    screen:	      camera resolution
    				Default: 640x480
    centered:	Is the hologram centered relative to the camera. Default is False
    ell:         Topological chrge of each point in the array. Default: 0
    win:         sinc effect. Default is False (0)
    alpha:       Relative trap intensities. Default is 1.
               
    OUTPUTS:
          phi: phase pattern encoding pattern of traps described by
               points with values ranging from [0,2 pi].
    
    RESTRICTIONS:
          Can be very memory intensive for large numbers of points.
    
    PROCEDURE:
          Initial estimate is obtained by superposing the fields of the
          specified beams.
    
     NOTES:
          If we calibrate the SLM's phase transfer function, then we
          should be able to pass the lookup table to this function, and
          return a hologram of indices into the look-up table.
    
    MODIFICATION HISTORY:
     Created by David G. Grier, New York University, 12/12/2004.
     1/19/2005: DGG.  Implemented 3D.
     1/29/2005: DGG.  Major code clean-up leading to improved speed.
                      Implemented SUBSAMPLE for major speedup.
     11/5/2011: Maya Yevnin -  added wavelength parameterand calibration files
     Mar 2013: Maya Yevnin -   Moved calibration parameters from the xmat, ymat
                               vectors to new vector corrected_points, and
                               added zfac parameter
     May 2015: Maya Yevnin, translated to python
    """
    [ndims,npts] = points.shape
    
    # calibration
    cal = loadcal(cal=cal)
    corrected_points = applycal(points, cal, centered, screen)
    
    # SLM coordinates
    if dim is None:
        from Holograms import h,w
    else:
        w = dim[0]
        h = dim[1]

    x = np.arange(w)
    y = np.arange(h)
    xmat, ymat = np.meshgrid(x - w/2 - cal['xc'], y - h/2 - cal['yc'])        

    # wavevectors associated with trap position (times i)
    ikx = 1j * (2*np.pi/w) * corrected_points[0,]
    iky = 1j * (2*np.pi/h) * corrected_points[1,]
    
    ikz = np.zeros(npts)
    rsq = np.zeros((h,w))
    if ndims > 2:
        aperture = 5000                            # 5mm aperture
        if wavelength is None:
            from Holograms import wavelength
        wavelength = (wavelength/1000)*w/aperture; # wavelength in aperture pixels
        na = 1.4                                   # numerical aperture of objective
        f = aperture/na                            # focal length in pixels
        ikz = 1j * (2*np.pi) * corrected_points[2,]/(wavelength*f^2)
        
        rsq = xmat**2 + ymat**2
    
    # vortex creation
    if ell is None:
        itheta = 0
        ell = np.ones(npts)
    else:
        itheta = 1j * np.arctan2(xmat, ymat)
        
    # correct sinc    
    if win is not None and w!=0:
        qxx = np.logical_not(corrected_points[0,])*1e-6 + corrected_points[0,]
        qyy = np.logical_not(corrected_points[1,])*1e-6 + corrected_points[1,]
        sinc = np.sin(win*np.pi/w) / (win*np.pi*qxx/w) * np.sin(win*np.pi/h) / (win*np.pi*qyy/h)
        alpha = alpha/(sinc**2)
    
    # trap intenseties
    if np.size(alpha)<2:
        alpha = alpha * np.ones(npts)
    alpha = alpha/np.sum(alpha)    # normalize alpha
    
    iphase = 1j * 2*np.pi * np.random.rand(npts)
    psi = np.zeros((h,w))
    
    for n in xrange(0, npts-1):
        ikxx = ikx[n] * x + iphase[n]
        ikyy = iky[n] * y
        ikzz = ikz[n] * rsq
        ex = np.exp(-ikxx).reshape([1,w])
        ey = np.exp(-ikyy).reshape([1,h])
        ez = np.exp(-ikzz)
        etheta = np.exp(-itheta * ell[n])
        
        psi = psi + np.sqrt(alpha[n]) * ex*ey.transpose() * ez * etheta
    phi = np.angle(psi) + np.pi
    return phi
        
    