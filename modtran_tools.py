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

def load_flx(filepath):

    fid = open(filepath,'r')
    hdr = list()
    for i in range(16):
       temp = fid.readline()
       if temp.strip('\n- '):
           hdr.append(temp.strip('\n -'))
    fid.close()

    # Get Altitude Data
    temp = hdr[3].split()
    n = len(temp)
    altitudes = [temp[i] for i in np.arange(1,n,2,dtype = int)]

    # Get Flux Data
    rawdata = np.genfromtxt(filepath,skip_header = 16,skip_footer = 2)
    n = rawdata.shape[1]

    direct_downwelling = [rawdata[:,i] for i in np.arange(3,n,3,dtype = int)]
    diffuse_downwelling = [rawdata[:,i] for i in np.arange(2,n,3,dtype = int)]
    upwelling = [rawdata[:,i] for i in np.arange(1,n,3,dtype = int)]
    wvl = rawdata[:,0]

    # conv = 10000   # Scale Factor to convert [W cm^-2 nm^-1] to [W m^-2 nm^-2], used for .flx MODTRAN output file
    downwelling = list()
    for i in range(len(altitudes)):
        downwelling.append(direct_downwelling[i] + diffuse_downwelling[i])

    ret_dict = dict([('header',hdr),
                     ('wvl',wvl),
                     ('altitude',altitudes),
                     ('downwelling_direct', direct_downwelling),
                     ('downwelling_diffuse',diffuse_downwelling),
                     ('upwelling',upwelling),
                     ('downwelling',downwelling)])
    return ret_dict

def load_7sc(filepath):
    data = np.genfromtxt(filepath,dtype = float, skip_header = 11, skip_footer = 1)
    # wvl_7sc = data_7sc.data[:,0]
    # r0 = data_7sc.data[:,9]
    # Iup = data_7sc.data[:,10]
    ret_dict = dict([('wvl',data[:,0]),
                     ('path_radiance',data[:,8]),
                     ('path_irradiance',data[:,9])])
    return ret_dict

def load_acd(filepath):
    # Call Function and specify number of lines in header to skip
    rawdata = np.loadtxt(filepath,skiprows = 5)
    # Creates a dictionary containing each of the numpy arrays corresponding
    #   to each term in the file
    ret_dict = dict([('freq',rawdata[:,0]),     # Frequency (cm-1)
                     ('wvl',1e7/rawdata[:,0]),  # Converts Frequency (cm-1) to wavelength [nm]
                     ('los',rawdata[:,1]),      # ?
                     ('kint',rawdata[:,2]),     # ?
                     ('kweight',rawdata[:,3]),  # ?
                     ('ts',rawdata[:,4]),       # Sun to Ground Diffuse Transmittance
                     ('Tso', rawdata[:,5]),     # Sun to Ground to Observer Direct Transmittance
                     ('t', rawdata[:,6]),       # Observer to Ground Embedded Diffuse Transmittance
                     ('T', rawdata[:,7]),       # Observer to Ground Direct Transmittance
                     ('sph', rawdata[:,8])])    # Spherical Albedo from Ground

    # Derive Ts and resolve any divide by 0, or NaN ambiguities
    ret_dict['Ts'] = ret_dict['Tso']/ret_dict['T'] # Sun to Ground Direct Transmittance
    ret_dict['Ts'][np.isnan(ret_dict['Ts'])] = 0
    ret_dict['Ts'][np.isinf(ret_dict['Ts'])] = 0

    # Returns the constructed dictionary
    return ret_dict

def load_modtran(filepath):

    modtran_model = 'Not Yet IMPLEMENTED'

    return modtran_model

if __name__ == '__main__':
    transfer_model = model('/Users/wrightad/Documents/MODTRAN/NEON_ATMCORR/','NEON_20150608_Baseline')
