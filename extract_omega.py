#!/usr/bin/env python

from netcdf_util import netcdf_file

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import os


BIGNUM = 1e200

def get_gamma_from_dir(dirname,tfit=80.0):

    filename = 'stella.out.nc'
    
    tmp = os.getcwd()

    os.chdir(dirname)

    if not os.path.isfile(filename):
        print("No output found in '" + dirname + "'.")
        os.chdir(tmp)
        return np.nan, np.nan
    
    with netcdf_file(filename,'r',mmap=False) as f:
        phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        t = f.variables['t'][()]
        
    gamma = []
    phi2s = np.transpose(phi2_vs_kxky[:,0,:])
    for iky, phi2 in enumerate(phi2s):
        i  = (np.abs(t-tfit)).argmin()
        iinf = np.where(phi2 > BIGNUM)[0]
        if len(iinf) > 0:
            y  = np.log(phi2[i:iinf[0]])
            x = t[i:iinf[0]]
        else:
            y  = np.log(phi2[i:])
            x = t[i:]
        N = np.count_nonzero(~np.isnan(y))
        if N >= 3:
            coefs0 = np.polyfit(x,y,1,full=False)[0]
        else:
            # no growth rate defined
            print("Warning: no growth rate for iky = " + str(iky))
            coefs0 = np.nan
        gamma.append(coefs0/2)
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


    exclude_dirs = []
    for dirname in d:
        print(dirname)
        ky,gamma = get_gamma_from_dir(dirname,tfit=t)
        if np.all(np.isnan(ky)):
            exclude_dirs.append(dirname)
        else:
            plt.plot(ky,gamma)

    i = 0
    N = len(d)
    for j in range(N):
        if d[i] in exclude_dirs:
            del d[i]
            i = i - 1
        i = i + 1
    print(d)
    plt.legend(d)
    plt.xlabel(r"$k_y$")
    plt.ylabel(r"$\gamma$")
    plt.show()
