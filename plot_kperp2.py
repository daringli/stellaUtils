#!/usr/bin/env python

from scipy.io import netcdf

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from stella_input import Stella_input

from scipy.signal import find_peaks

import sys


with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
    # phi_vs_t(t, tube, zed, kx, ky, ri)
    phi = f.variables['phi_vs_t'][()]
    print(phi.shape)
    phi_vs_z  = phi[-1,0,:,:,:,:] # last time, first tube, real and imaginary
    print(phi_vs_z.shape)
    phi2_vs_z  = np.sum(phi_vs_z**2,axis=-1) # sum real and imaginary squared
    print(phi2_vs_z.shape)
    ky = f.variables['ky'][()]
    kx = f.variables['kx'][()]
    z = f.variables['zed'][()]
    kperp2 = f.variables['kperp2'][()]
    kperp2 = kperp2[:,0,0,:]

nfield_periods = Stella_input('.').nfield_periods
z = z * nfield_periods

    
fig,axes = plt.subplots(3,sharex=True)

if len(kx) == 1:

    Nky = len(ky)
    cmap = cm.get_cmap("rainbow",Nky)
    ikx = 0
    # last time-step
    for iky,aky in enumerate(ky):
        norm = np.max(phi2_vs_z[:,ikx,iky])
        product = phi2_vs_z[:,ikx,iky]*kperp2[:,iky]
        norm2 = np.max(product)
        axes[0].plot(z,phi2_vs_z[:,ikx,iky]/norm,color=cmap(iky))
        axes[2].plot(z,product/norm2,color=cmap(iky))
    kperp2_normalized = kperp2[:,-1]/ky[-1]**2 
    axes[1].plot(z,kperp2_normalized)
    for i in find_peaks(kperp2_normalized)[0]:
        axes[0].axvline(z[i],linestyle='dashed',color='silver')
        axes[1].axvline(z[i],linestyle='dashed',color='silver')
        axes[2].axvline(z[i],linestyle='dashed',color='silver')


    #plt.legend(ky)

    axes[-1].set_xlabel(r"$z$")
    axes[0].set_ylabel(r"$|\phi|_k^2/\rm{max}(|\phi|_k^2)$")
    axes[1].set_ylabel(r"$k_{\perp}^2$")
    axes[2].set_ylabel(r"$|\phi|_k^2 k_{\perp}^2$")
    norm = mpl.colors.Normalize(vmin=ky[0], vmax=ky[-1])
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

    cax.set_title(r'$k_y \rho$')
    
    divider1 = make_axes_locatable(axes[1])
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax1.remove()

    divider2 = make_axes_locatable(axes[2])
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

    cax2.set_title(r'$k_y \rho$')
    
    if len(sys.argv) > 1:
        plt.savefig("kperp2.png")
    else:
        plt.show()
