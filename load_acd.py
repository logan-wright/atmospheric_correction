#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:46:57 2017

Loads a text file
Can be called in a script by:
    output = load_acd(filepath)

@author: wrightla
"""
import numpy as N

def load_acd(filepath):
    # Call Function and specify number of lines in header to skip
    rnyawdata = N.loadtxt(filepath,skiprows = 5) 
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
    
if __name__ == '__main__':
    output = load_acd('H2010159042916_Clearsky.acd')
    