#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:33:01 2017

@author: wrightad
"""
import spectral
import numpy as N
import scipy.interpolate as sci_interp
from skimage import exposure
import matplotlib.pyplot as plt

#import glob
#import h5py
#import numpy as N
#import scipy.interpolate as sci_interp
#from skimage import exposure
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from datetime import datetime
#from mpl_toolkits.mplot3d import Axes3D

#import bitget

class datacube(object):
    def __init__(self,data_cube,resp_func,nav,flags,rgb):
        self.data_cube = data_cube
        self.resp_func = resp_func
        self.nav = nav
        self.flags = flags
        self.rgb = rgb

def load_nis(file):
    img = spectral.io.envi.open(file + '.hdr',file)
    data = img.load()
    data[data == -9999] = N.NaN
    # Mask Invalid Data and Divide by 100 to convert to W/sr/m^2/nm
    scaled_data = N.ma.masked_invalid(data)/100

    flags = None
    nav = None
    bandinfo = data.bands
    resp = dict([('units','W/m^2/sr/nm'),('wvl',bandinfo.centers),('fwhm',bandinfo.bandwidths)])

    # Load NIS Metadata
    # Load metadata
    img = spectral.io.envi.open(file + '.hdr',file)
    data = img.load()
    nav = dict([('gps_time',metadata[:,:,10])])
    # data[data == -9999] = N.NaN

    # Produce and RGB "Quasi-truecolor" preview image
    # Apply 2% Linear Stretch to "True Color"
    R = scaled_data[:,:,50]/N.amax(scaled_data[:,:,50])
    p2, p98 = N.nanpercentile(R, (1, 98))
    Rscl = exposure.rescale_intensity(R, in_range=(p2, p98))

    G = scaled_data[:,:,35]/N.amax(scaled_data[:,:,35])
    p2, p98 = N.nanpercentile(G,(1, 98))
    Gscl = exposure.rescale_intensity(G, in_range=(p2, p98))

    B = scaled_data[:,:,19]/N.amax(scaled_data[:,:,19])
    p2, p98 = N.nanpercentile(B, (1, 98))
    Bscl = exposure.rescale_intensity(B, in_range=(p2, p98))

    rgb = N.dstack((Rscl,Gscl,Bscl))

    ''' Debugging Block
    plt.figure()
    plt.imshow(N.stack((R,G,B), axis = 2))
    plt.show()
    plt.figure()
    plt.imshow(rgb)
    plt.show()
    '''

    # Double check to make sure all values are in range, set any out of range
    #   values to the range limits
    rgb[N.greater(rgb,1)] = 1
    rgb[N.less(rgb,0)] = 0

    # Return
    ret_val = datacube(scaled_data,resp,nav,flags,rgb)
    return ret_val

if __name__ == '__main__':
    out = load_nis('/Users/wrightad/Documents/Data/NEON/Flight_ROIs/20150608/NIS01_20150608_165842_rad')
