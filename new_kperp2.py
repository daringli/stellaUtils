#!/usr/bin/env python


from netcdf_util import netcdf_file

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from stella_input import Stella_input

from scipy.signal import find_peaks
from scipy.interpolate import interp1d

import sys



def get_kperp2(dirname, override_ky = None):
    stella_output = dirname + "/" + 'stella.out.nc'
    with netcdf_file(stella_output,'r',mmap=False) as f:
        # phi_vs_t(t, tube, zed, kx, ky, ri)
        phi = f.variables['phi_vs_t'][()]
        phi_vs_z  = phi[-1,0,:,:,:,:] # last time, first tube, real and imaginary
        phi2_vs_z  = np.sum(phi_vs_z**2,axis=-1) # sum real and imaginary squared. zed, kx, ky
        if override_ky is None:
            ky = f.variables['ky'][()]
        else:
            _ky = f.variables['ky'][()]
            _ky = np.concatenate(([0.0], _ky, [_ky[-1] * 2]))
            tmp_zeros = np.zeros_like(phi2_vs_z[:,:,0])
            tmp_zeros = np.array([tmp_zeros])
            #print(tmp_zeros.shape)
            tmp_zeros = np.transpose(tmp_zeros, axes= [1,2,0])
            #print(tmp_zeros.shape)
            #print(phi2_vs_z.shape)
            phi2_vs_z = np.concatenate((tmp_zeros, phi2_vs_z , tmp_zeros), axis=2)
            ky = override_ky
            phi2_vs_z = interp1d(_ky,phi2_vs_z)(ky)
        kx = f.variables['kx'][()]
        I = np.argsort(kx)
        kx = kx[I]
        #print(kx)
        phi2_vs_z = phi2_vs_z[:,I,:]
        z = f.variables['zed'][()]
        #kperp2 = f.variables['kperp2'][()]
        #kperp2 = kperp2[:,0,0,:]
        gds2 = f.variables['gds2'][()] # nzed, nalpha
        gds21 = f.variables['gds21'][()] # nzed, nalpha
        gds22 = f.variables['gds22'][()] # nzed, nalpha
        shat = f.variables['shat'][()] # scalar
        jacob = f.variables['jacob'][()] # nzed, nalpha    
        

        # used to only work for Nkx = 1
        #phi2_vs_z = phi2_vs_z[:,0,:]

        # this is really kperp^2 times phi2 and the jacobian, as needed for the flux surface average
        # '0' is since we only use the first flux tube
        ky = ky[None,None,:]
        kx = kx[None,:, None]
        gds2 = gds2[:,0,None,None]
        gds21 = gds21[:,0,None,None]
        gds22 = gds22[:,0,None,None]
        jacob = jacob[:,0,None,None]
        kperp2 = (ky**2 * gds2 + 2*kx*ky*gds21/shat + kx**2 * gds22/shat**2) * phi2_vs_z * jacob/np.sum(phi2_vs_z * jacob,axis=0)
        kperp2 = np.sum(kperp2,axis=0)
        return kperp2

