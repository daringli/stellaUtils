#!/usr/bin/env python

from netcdf_util import netcdf_file

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from stella_input import Stella_input

def get_latest_phi2(dirname):
    with netcdf_file(dirname + '/stella.out.nc','r',mmap=False) as f:
        # phi_vs_t(t, tube, zed, kx, ky, ri)
        phi = f.variables['phi_vs_t'][()]
        print(phi.shape)
        phi_vs_z  = phi[-1,0,:,:,:,:] # last time, first tube, real and imaginary
        print(phi_vs_z.shape)
        phi2_vs_z  = np.sum(phi_vs_z**2,axis=-1) # sum real and imaginary squared
        phi2_vs_z  = np.sum(phi2_vs_z,axis=(1,2)) # over ky and kx

        print(phi2_vs_z.shape)

        z = f.variables['zed'][()]
    nfield_periods = Stella_input(dirname).nfield_periods
    return z * nfield_periods, phi2_vs_z

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']

    fig, ax =plt.subplots(1,sharex=True)
    for d in dirs:
        z,phi2 = get_latest_phi2(d)
        

        norm = np.max(phi2)
        ax.plot(z,phi2/norm)
        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"$|\phi|^2/\rm{max}(|\phi|^2)$")
    ax.legend(dirs)
    plt.show()
