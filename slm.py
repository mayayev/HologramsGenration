# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 15:55:06 2015

@author: user
"""

from __future__ import division, unicode_literals

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def fullscreen(final_image):
    mpl.rcParams['toolbar'] = 'None'
    fig = plt.figure('slm')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    plt.imshow(final_image, cmap='gray')
    plt.show()
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    


def slm(phi, alpha=None, correction_im=None, SLM_screen=2):


    if correction_im is None:
        from Holograms import correction_im

    if alpha is None:
        from Holograms import alpha

    final_image = phi + correction_im
    final_image = np.mod(final_image, 2*np.pi)
    final_image = final_image * (alpha/255)
    fullscreen(final_image)




