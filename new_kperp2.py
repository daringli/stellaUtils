#!/usr/bin/env python


from netcdf_util import netcdf_file

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from stella_input import Stella_input

from scipy.signal import find_peaks

import sys



def get_kperp2(dirname):
    stella_output = dirname + "/" + 'stella.out.nc'
    with netcdf_file(stella_output,'r',mmap=False) as f:
        # phi_vs_t(t, tube, zed, kx, ky, ri)
        phi = f.variables['phi_vs_t'][()]
        phi_vs_z  = phi[-1,0,:,:,:,:] # last time, first tube, real and imaginary
        phi2_vs_z  = np.sum(phi_vs_z**2,axis=-1) # sum real and imaginary squared
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        z = f.variables['zed'][()]
        kperp2 = f.variables['kperp2'][()]
        kperp2 = kperp2[:,0,0,:]
        gds2 = f.variables['gds2'][()] # nzed, nalpha
        jacob = f.variables['jacob'][()] # nzed, nalpha    
        

        phi2_vs_z = phi2_vs_z[:,0,:]
        kperp2 = ky**2 * gds2[:,0,None] * phi2_vs_z * jacob[:,0, None]/np.sum(phi2_vs_z * jacob,axis=0)
        kperp2 = np.sum(kperp2,axis=0)
        return kperp2

