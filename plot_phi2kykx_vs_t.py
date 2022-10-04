#!/usr/bin/env python

from scipy.io import netcdf

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


import sys

<<<<<<< HEAD
=======
BIGNUM = 1e200
>>>>>>> 3fb6fbc10cc8dc1c280b62e4fdd31cb99437c535

with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
    phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
    ky = f.variables['ky'][()]
    kx = f.variables['kx'][()]
    t = f.variables['t'][()]
<<<<<<< HEAD

=======
    # remove inf parts of phi2
    iinf = np.where(phi2_vs_kxky > BIGNUM)[0]
    if len(iinf) > 0:
        t = t[:iinf[0]]
        phi2_vs_kxky = phi2_vs_kxky[:iinf[0]]
    
>>>>>>> 3fb6fbc10cc8dc1c280b62e4fdd31cb99437c535
lines = []
gamma = np.zeros_like(ky)
def onclick(event):
    if event.inaxes is not None:
        ax = event.inaxes

    if ax == axes[0]:
        global gamma
        gamma = []
        clicked_t = event.xdata
        phi2s = np.transpose(phi2_vs_kxky[:,0,:])
        for phi2 in phi2s:
            i  = (np.abs(t-clicked_t)).argmin()
            y  = np.log(phi2[i:])
            x = t[i:]
            coefs = np.polyfit(x,y,1,full=False)
            gamma.append(coefs[0]/2)
        gamma = np.array(gamma)
        
        line = axes[0].axvline(clicked_t)
        lines.append(line)
        if len(lines) > 1:
            lines[0].remove()
            del lines[0]
        gamma_line.set_ydata(gamma)
        #update limits must be done manually
        axes[1].relim()
        axes[1].autoscale()
        axes[1].set_ylim(bottom=0.0)
        plt.draw()
        np.save('my_growthrate.npy',gamma)
        np.save('my_ky.npy',ky)
        np.save("my_tfit.npy",t[i])
        plt.savefig("phi2kykx_vs_t.png")

fig,axes = plt.subplots(2)

cid = fig.canvas.mpl_connect('button_press_event', onclick)

if len(kx) == 1:

    Nky = len(ky)
    cmap = cm.get_cmap("rainbow",Nky)

    # last time-step
    for iky,aky in enumerate(ky):
        phi2_vs_t  = phi2_vs_kxky[:,0,iky]
        axes[0].semilogy(t,phi2_vs_t,color=cmap(iky))
    #plt.legend(ky)

    axes[0].set_xlabel(r"$t$")
    axes[0].set_ylabel(r"$|\phi|_k^2$")
    norm = mpl.colors.Normalize(vmin=ky[0], vmax=ky[-1])
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

    #cb1.set_label(r'$k_y \rho$', labelpad=-40, y=1.05, rotation=0)
    cax.set_title(r'$k_y \rho$')
    gamma_line = axes[1].plot(ky,gamma)
    gamma_line = gamma_line[0]
    axes[1].set_xlabel(r"$k_y \rho$")
    axes[1].set_ylabel(r"$\gamma$")
    axes[1].set_ylim(bottom=0.0)

    #plt.show()
    plt.show()
