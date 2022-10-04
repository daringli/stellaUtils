#!/usr/bin/env python

from scipy.io import netcdf

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import os




def get_gamma_from_dir(dirname,tfit=80.0):

    tmp = os.getcwd()

    os.chdir(dirname)

    with netcdf.netcdf_file('stella.out.nc','r',mmap=False) as f:
        phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        t = f.variables['t'][()]
    
    gamma = []
    phi2s = np.transpose(phi2_vs_kxky[:,0,:])
    for phi2 in phi2s:
        i  = (np.abs(t-tfit)).argmin()
        y  = np.log(phi2[i:])
        x = t[i:]
        coefs = np.polyfit(x,y,1,full=False)
        gamma.append(coefs[0]/2)
    gamma = np.array(gamma)
    np.save('my_growthrate.npy',gamma)
    np.save('my_ky.npy',ky)
    np.save("my_tfit.npy",tfit)

    os.chdir(tmp)
    return ky, gamma

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    argv = sys.argv
    argc = len(argv)
    if argc > 2:
        t = float(argv[1])
        d = argv[2:]

    for dirname in d:
        print(d)
        ky,gamma = get_gamma_from_dir(dirname,tfit=t)
        plt.plot(ky,gamma)


    plt.legend(d)
    plt.show()
