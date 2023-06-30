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
    full_kx = False
    filename = 'stella.out.nc'
    
    tmp = os.getcwd()

    os.chdir(dirname)

    if not os.path.isfile(filename):
        print("No output found in '" + dirname + "'.")
        os.chdir(tmp)
        return np.nan, np.nan, np.nan
    
    with netcdf_file(filename,'r',mmap=False) as f:
        phi2_vs_kxky  = f.variables['phi2_vs_kxky'][()]
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        t = f.variables['t'][()]

    Nkx = len(kx)
    Nky = len(ky)
    phi2s = np.transpose(phi2_vs_kxky,axes=[1,2,0])  # kx, ky, t

    if not full_kx:
        # remove negative kx
        phi2s = phi2s[:(Nkx//2+1)]
        kx = kx[:(Nkx//2+1)]
        Nkx = Nkx//2 + 1
        kx0_index = 0
    else:
        I = np.argsort(kx)
        kx = kx[I]
        phi2s = phi2s[I]
        kx0_index = Nkx//2 + 1

    gamma = np.zeros((Nkx,Nky))
    iinfs = np.full((Nkx,Nky),None)


    for ikx, phi2_x in enumerate(phi2s):
        for iky, phi2 in enumerate(phi2_x):
            imin  = (np.abs(t-tfit)).argmin()
            iinf = np.where(phi2 > BIGNUM)[0]
            if len(iinf) == 0:
                iinf = np.where(np.isnan(phi2))[0]
            if len(iinf) > 0:
                imax = iinf[0]
            else:
                imax = None
            y  = np.log(phi2[imin:imax])
            x = t[imin:imax]
            N = np.count_nonzero(~np.isnan(y))
            if N >= 3:
                coefs0 = np.polyfit(x,y,1,full=False)[0]
            else:
                # no growth rate defined
                print("Warning: no growth rate for iky = " + str(iky))
                coefs0 = np.nan
            gamma[ikx,iky] = coefs0/2

    np.save('my_growthrate.npy',gamma)
    np.save('my_ky.npy',ky)
    np.save('my_kx.npy',kx)
    np.save("my_tfit.npy",tfit)

    os.chdir(tmp)
    return ky, kx, gamma

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
        ky,kx, gamma = get_gamma_from_dir(dirname,tfit=t)
        if np.all(np.isnan(ky)):
            exclude_dirs.append(dirname)
        else:
            plt.plot(ky,gamma[0])

            
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
    #plt.show()
