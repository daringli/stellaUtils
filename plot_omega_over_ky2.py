#!/usr/bin/env python

from scipy.io import netcdf

import os
import numpy as np
import matplotlib.pyplot as plt
from plot_omega import omega_from_dir

from stella_input import Stella_input

path = os.path.abspath(__file__)
thisfile = path.split('/')[-1]

thisdir = os.getcwd().split('/')[-1]

def get_FSAkperp2(dirname):

    if not os.path.isfile(dirname + '/stella.out.nc'):
        return np.nan, np.nan

    ikx = 0
    with netcdf.netcdf_file(dirname + '/stella.out.nc','r',mmap=False) as f:
        # phi_vs_t(t, tube, zed, kx, ky, ri)
        phi = f.variables['phi_vs_t'][()]
        jacob = f.variables['jacob'][()]
        
        #print(phi.shape)
        phi_vs_z  = phi[-1,0,:,:,:,:] # last time, first tube, real and imaginary
        #print(phi_vs_z.shape)
        phi2_vs_z  = np.sum(phi_vs_z**2,axis=-1) # sum real and imaginary squared
        #print(phi2_vs_z.shape)
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        z = f.variables['zed'][()]
        

        # should we use this instead?
        nfield_periods = Stella_input(dirname).nfield_periods
        z2 = z * nfield_periods

        dz = z[1] - z[0]

        print("These will all be the same for uniform z spacing:")
        print(np.min(np.diff(z)),np.max(np.diff(z)))

        kperp2 = f.variables['kperp2'][()]
        kperp2 = kperp2[:,0,ikx,:] #kperp2(zed, alpha, kx, ky)

        product = phi2_vs_z[:,ikx,:]*kperp2[:,:] # zed, ky

        # sum / integrate over z:
        # just summing kperp and phi_vs_t from output .nc file
        # this is correct if phi or kperp2 do no include the "jacobian" dl/dz 1/|B|
        FSAkperp2 = np.sum(product, axis = 0)/np.sum(phi2_vs_z[:,ikx,:], axis =0)
    
        # weight  the sums with dl/|B| as in a flux-surface average
        # this is correct if phi or kperp2 do no include the "jacobian" dl/dz 1/|B|
        dl_over_B = jacob[0,:] * dz # dz does not matter as it cancels out, for uniform z spacing
        FSA2kperp2 = np.sum(product * dl_over_B, axis = 0)/np.sum(phi2_vs_z[:,ikx,:] * dl_over_B, axis =0)
    
        return FSAkperp2, FSA2kperp2

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']


    
    fig, axes =plt.subplots(4,sharex=True)
    for d in dirs:

        print(d)
        FSAkperp2, FSA2kperp2 = get_FSAkperp2(d)

        if np.any(np.isnan(FSAkperp2)):
            FSAkperp2, FSA2kperp2 = np.nan, np.nan
            ky, omega, gamma = np.nan, np.nan, np.nan
        else:
            ky, omega, gamma = omega_from_dir(d)
            if len(gamma.shape) > 1:
                if gamma.shape[1] == 1:
                    gamma = gamma[:,0]
        
        print(ky)
        print(gamma)
        print(FSAkperp2)
        axes[0].plot(ky, gamma)
        axes[1].plot(ky, gamma/ky**2)
        axes[2].plot(ky, gamma/FSAkperp2)
        axes[3].plot(ky, gamma/FSA2kperp2)
        


    axes[0].set_ylabel(r"$\gamma$")
    axes[1].set_ylabel(r"$\gamma/k_y^2$")
    axes[2].set_ylabel(r"$\gamma/<k_\perp^2>$")
    axes[3].set_ylabel(r"$\gamma/<k_\perp^2>_2$")
    
    axes[-1].set_xlabel(r"$k_y \rho$")
    axes[0].set_ylim(bottom=0.0)
    axes[1].set_ylim(bottom=0.0)
    axes[2].set_ylim(bottom=0.0)
    axes[3].set_ylim(bottom=0.0)
    
    axes[1].legend(dirs)    

    axes[0].set_title("$STELLAUTILS/" + thisfile + ", run from: " + thisdir)
    plt.show()
