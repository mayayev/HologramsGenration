from __future__ import division
import numpy as np
import scipy as sp
from holograms.ApplyCalibration import loadcal, applycal


def dsphase(points, wavelength=None, dim=None, cal=None, screen=[640,480], centered=False, ell=None, win=None, alpha=1, rho=0.5, thresh=0.00005):

    ## Satge 1: initial hologram. Identical to phasemask    
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
            from holograms import wavelength
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
#    iphase = sp.io.loadmat('D:\Dropbox\Matlab\Holograms\dsphase\iphase.mat')['iphase']
#    iphase = iphase.squeeze()
    psi = np.zeros((h,w))
    
    for n in xrange(0, npts):
        ikxx = ikx[n] * x + iphase[n]
        ikyy = iky[n] * y
        ikzz = ikz[n] * rsq
        ex = np.exp(-ikxx).reshape([1,w])
        ey = np.exp(-ikyy).reshape([1,h])
        ez = np.exp(-ikzz)
        etheta = np.exp(-itheta * ell[n])
        psi = psi + np.sqrt(alpha[n]) * ex*ey.transpose() * ez * etheta
    phi_initial = np.angle(psi) + np.pi
	
    # quantize
    phi_initial = (np.round(phi_initial/(2*np.pi) * 256))*2*np.pi/256
    phaselut = 1j*phi_initial
    expphaselut = np.exp(phaselut)


    ## stage 2: Compute initial fields & intensities of traps using the initial phase estimate from stage 1
    exp_all = np.zeros((npts,h,w)) + 0j   
    E = np.zeros(npts) + 0j
    for n in xrange(0, npts):
        ikxx = ikx[n] * x
        ikyy = iky[n] * y
        ikzz = ikz[n] * rsq
        # these are the propagators from each pixel to each trap
        ex = np.exp(ikxx).reshape([1,w])
        ey = np.exp(ikyy).reshape([1,h])
        ez = np.exp(ikzz)
        etheta = np.exp(itheta * ell[n])
        exp_all[n,:,:] = ex * ey.transpose() * ez * etheta
        E[n] = np.sum(np.sum(exp_all[n,:,:]*expphaselut))
    I = E*E.conj()    # Initial intensities
    
    # calc initial convergence factor & other trap array parameters
    gamma = np.sum(alpha*I) / np.sum(alpha**2)
    sigma = np.sqrt(np.mean((I - np.dot(gamma,alpha))**2))
    initial_conv = rho*sigma - (1-rho)*np.mean(I)
    conv = initial_conv
    efficiency = np.sum(I)/np.sum(alpha)
    rmserror= sigma/np.max(I)
    max_I = np.max(I/alpha)
    min_I = np.min(I/alpha)
    uniformity = (max_I-min_I)/(max_I+min_I)
    
    
    ## stage 3: direct search
    acc=0;
    phaselut_new = np.zeros((h,w)) + 0j
    
    # while (rmserror > thresh) or (uniformity > thresh):
    rms_maya = np.arange(w*h)+0j
    unif_maya = np.arange(w*h)+0j
    conv_maya = np.arange(w*h, dtype=float)+0j
    for r in xrange(1):
        for t in xrange(0, w*h):

            # choose a random pixel
            l = np.ceil(np.random.rand()* (h-1))
            m = np.ceil(np.random.rand()* (w-1))
    
            # change current pixel to a random phase (quantized)
    
            deltaphi = np.round(np.random.rand()*256)*2*np.pi/256
            phaselut_new[l,m] = deltaphi 
    
            # calc the new fields (maybe this can be done in matrix form to
            # save the for loop)
            E_new = np.zeros(npts) + 0j
            for q in xrange(0,npts):
                # subtract the old phase and add the new
                E_new[q] = E[q] + exp_all[q,l,m]*(-expphaselut[l,m] + np.exp(phaselut_new[l,m])) 
    
            # calc new intensities, conv factor and other trap array
            # parameters
            I_new = E_new*E_new.conj()
            I_new_norm = I_new/np.max(I_new)
    
            gamma_new = np.sum(alpha*I_new)/np.sum(alpha**2)
            sigma_new = np.sqrt(np.mean((I_new-gamma_new*alpha)**2))
            conv_new = rho*sigma_new - (1-rho)*np.mean(I_new)
            efficiency_new = np.sum(I_new_norm)/np.sum(alpha)
            rmserror_new = sigma/np.max(I_new)
            max_I_new = np.max(I_new/alpha)
            min_I_new = np.min(I_new/alpha)
            uniformity_new = (max_I_new-min_I_new)/(max_I_new+min_I_new)
            
            print t
            # print "sigma =", sigma_new
            # print "rms =", rmserror_new
            # print "uniformity =", uniformity_new
            # print "conv =", conv_new
            # print "\n"
            rms_maya[t] = rmserror_new
            unif_maya[t] = uniformity_new
            conv_maya[t] = conv_new
            
            if conv_new < conv: # if convergence has improved accept new values
                acc=acc+1;
                phaselut[l,m] = phaselut_new[l,m]
                expphaselut[l,m] = np.exp(phaselut_new[l,m])
                E = E_new;
                I = I_new;
                I_norm = I_new_norm;
                gamma = gamma_new;
                sigma = sigma_new;
                conv = conv_new;
                efficiency = efficiency_new;
                rmserror = rmserror_new;
                uniformity = uniformity_new;
    
                if  (rmserror < thresh) and (uniformity < thresh):
                    break
    # return np.abs(phaselut)
	return {'rms': rms_maya, 'uniformity':unif_maya, 'conv':conv_maya}

    

    
        
