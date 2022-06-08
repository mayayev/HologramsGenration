# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:42:06 2015

@author: user
"""
from __future__ import division
from displace import displace
from slm import slm
from vortex import vortex
from math import *
from fastphase import fastphase, phasemask
from dsphase import dsphase

# default wavelength
wavelength = 532
w = 800
h = 600

# default SLM calibration
from scipy.ndimage.io import imread
import numpy as np
import os

path = os.path.dirname(os.path.realpath(__file__))
if wavelength == 532:
    correction_im = imread(path + '\HOTs_correction_images\CAL_LSH0400106_532nm.bmp')
    alpha = 208
#elif wavelength == 1085:
#    correction_im=imread('E:\Dropbox\Matlab\Holograms\SLM\HOTs_correction_images\CAL_LSH0300104_1085nm.bmp')
#    alpha = 215
else:
    correction_im = np.zeros((600,792))
    alpha = 255

correction_im = np.pad(correction_im, ((0,0),(8,0)), 'constant', constant_values=0)
correction_im = correction_im * (2*np.pi)/255


del imread
del np
del os
del path