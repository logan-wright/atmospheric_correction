#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:46:57 2017

Loads a text file
Can be called in a script by:
    output = load_acd(filepath)

@author: wrightla
"""


import numpy as np

class model(object):
    def __init__(self,directory,name):
        self.irradiance = load_flx(directory+name+'.flx')
        self.transmittance = load_acd(directory+name+'.acd')
        # TO BE IMPLEMENTED LATER
        #self.atmosphere = load_tp6(directory+name+'.tp6')
        #self.input = load_tp5(directory+name+'.tp5')

def load_flx(file):

    fid = open(file)
    hdr = list()
    for i in range(16):
       temp = fid.readline()
       if temp.strip('\n- '):
           hdr.append(temp.strip('\n -'))
    fid.close()
    
    temp = hdr[3].split()
    n = len(temp)
    
    altitudes = [temp[i] for i in np.arange(0,n,2,dtype = int)]
    rawdata = np.genfromtxt(file,skip_header = 16,skip_footer = 2)
    
    direct_downwelling = [rawdata[:,i] for i in np.arange(3,n,3,dtype = int)]
    diffuse_downwelling = [rawdata[:,i] for i in np.arange(2,n,3,dtype = int)]
    upwelling = [rawdata[:,i] for i in np.arange(1,n,3,dtype = int)]

    ret_dict = dict([(
                     (
                     ()))])
    return ret_dict

def load_7sc(filepath):


    return ret_dict

def load_acd(filepath):
    # Call Function and specify number of lines in header to skip
    rawdata = np.loadtxt(filepath,skiprows = 5)
    # Creates a dictionary containing each of the numpy arrays corresponding
    #   to each term in the file
    ret_dict = dict([('freq',rawdata[:,0]),
                     ('los',rawdata[:,1]),
                     ('kint',rawdata[:,2]),
                     ('kweight',rawdata[:,3]),
                     ('ts',rawdata[:,4]),
                     ('Ts', rawdata[:,5]),
                     ('t', rawdata[:,6]),
                     ('T', rawdata[:,7]),
                     ('sph', rawdata[:,8])])

    # Returns the construction dictionary
    return ret_dict

def load_modtran(filepath):

    return modtran_model

if __name__ == '__main__':
    transfer_model = model('/Users/wrightad/Documents/MODTRAN/NEON_ATMCORR/','NEON_20150608_Baseline')
