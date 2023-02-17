#!/usr/bin/env python


import numpy as np
from netcdf_util import netcdf_file
from new_kperp2 import get_kperp2


BIG_NUM = 1.0

def get_gamma(dirname, tfit=20.0):
    with netcdf_file(dirname + '/stella.out.nc','r',mmap=False) as f:
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
    return ky, gamma


def rogerio_objective(dirname):
    try:
        ky, gamma = get_gamma(dirname)
    except FileNotFoundError:
        print("failed to extract gamma from: " + dirname)
        return BIG_NUM
    kperp2 = get_kperp2(dirname)
    Nky = len(kperp2)
    return np.trapz(gamma/kperp2,ky)



if __name__=="__main__":
    # Print the value of the rogerio objective of a stella simulation.
    print("start")
    import sys
    argc = len(sys.argv)
    if argc > 1:
        dirname = sys.argv[1]
    else:
        dirname = '.'
    val = rogerio_objective(dirname)
    print("2222?")
    print(val)
