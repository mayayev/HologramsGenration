from __future__ import division


def displace(dx, dy, dz=0, dim=None, wavelength=None, cal=None, screen=[640,480], centered=False):
    """
    Calculates a phase grating which displaces an optical trap from the center of the field of view to any desired
    location in three dimensions. Can be used to displace the traps in a previously calculated HOT DOE.

    CALLING SEQUENCE:
    phi = offset( dx, dy, [dz] )

    INPUTS:
    dx, dy: in-plane offsets measured in pixels.

    OPTIONAL INPUT:
    dz:         out-of-plane offset in pixels.  Default: dz=0
    w,h:        dimensions of phase hologram.   Default: w=600, h=800
    cal:        calibration parameters.
    wavelength: laser wavelength.               Default: wavelength=1085
    screen:     camera resolution.              Default: 480x640
    centered:   if set to True, then the hologram will be shifted to the center of the screen. Default is False.


    OUTPUTS:
    phi: Real-valued phase grating in the range 0 to 2 pi.

    PROCEDURE:
    Appropriately scaled plane wave and Fresnel lens.

    Created on Thu Apr 16 14:29:31 2015

    @author: Maya Yevnin
    """

    import numpy as np
    from Holograms.ApplyCalibration import loadcal, applycal
    import Holograms as holo
    
    cal = loadcal(cal=cal)
    corrected_displace = applycal(np.vstack([dx,dy,dz]), cal, centered, screen)
    dx = corrected_displace[0,0]
    dy = corrected_displace[1,0]
    dz = corrected_displace[2,0]
    if dim is None:
        w = holo.w
        h = holo.h
    else:
        w = dim[0]
        h = dim[1]
        
    x = np.arange(w)
    y = np.arange(h)
    xmat, ymat = np.meshgrid(x - w/2 - cal['xc'], y - h/2 - cal['yc'])

    phi = dx * xmat / w + dy * ymat / h

    if (dz != 0):
        aperture = 5000.                                    # 5 mm aperture on lens
        if wavelength is None:
            wavelength = holo.wavelength
        wavelength = (wavelength/1000.) * w / aperture      # wavelength in pixels at aperture
        na = 1.4                                            # numerical aperture
        f = aperture/na                                     # focal length in pixels
        rmat = xmat**2. + ymat**2.
        phi = phi - rmat * dz / (wavelength * f**2.)

    phi = 2. * np.pi * (phi - np.amin(phi))

    phi = np.mod(phi, 2*np.pi)

    return phi

