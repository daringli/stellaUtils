#!/usr/bin/env python

from netcdf_util import netcdf_file


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from stella_input import Stella_input




with netcdf_file('stella.out.nc','r',mmap=False) as f:
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
nfield_periods = Stella_input('.').nfield_periods
z = z * nfield_periods

    
fig,ax = plt.subplots(1)

if len(kx) == 1:

    Nky = len(ky)
    cmap = cm.get_cmap("rainbow",Nky)
    ikx = 0
    # last time-step
    for iky,aky in enumerate(ky):
        norm = np.max(phi2_vs_z[:,ikx,iky])
        ax.plot(z,phi2_vs_z[:,ikx,iky]/norm,color=cmap(iky))
    #plt.legend(ky)

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$|\phi|_k^2/\rm{max}(|\phi|_k^2)$")
    norm = mpl.colors.Normalize(vmin=ky[0], vmax=ky[-1])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

    cax.set_title(r'$k_y \rho$')
    plt.show()
