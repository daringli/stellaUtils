#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf


def kperp2_from_dir(dirname, usemy = True):
    tmp = os.getcwd()
    os.chdir(dirname)
    with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
        kperp2 = f.variables['kperp2'][()]
        jacob = f.variables['jacob'][()]
        zed = f.variables['zed'][()]
        delzed  =zed[1] - zed[0]
        ky = f.variables['ky'][()]
    kperp2 = kperp2[:,0,0,:] * jacob[0,:] * delzed
    
    os.chdir(tmp)
    return ky, kperp2 
    
if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']


    fig, ax =plt.subplots(1)
    axes = [ax]
    for d in dirs:
        ky, kperp2 = kperp2_from_dir(d)
        kperp2 = np.sum(kperp2,axis=0) # sum over z
        axes[0].plot(ky**2, kperp2)
        axes[0].set_xlabel(r"$k_y^2")
        axes[0].set_ylabel(r"$k_{\perp}^2$")
        
    axes[0].legend(dirs)    
    plt.show()
