#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from ql_fluxes import get_latest_omega, read_omega


def omega_from_dir(dirname, usemy = True):
    tmp = os.getcwd()
    os.chdir(dirname)
    
    if os.path.exists("my_growthrate.npy") and usemy:
        Imomega_avg = np.load("my_growthrate.npy")
        ky = np.load("my_ky.npy")
        kx = np.array([0])

        omega_tuple = read_omega()
        t, ky, kx, _Imomega, _Reomega, Reomega_avg, _Imomega_avg = get_latest_omega(omega_tuple)

    else:
        omega_tuple = read_omega()
        t, ky, kx, Imomega, Reomega, Reomega_avg, Imomega_avg = get_latest_omega(omega_tuple)

    #Ky,Kx = np.meshgrid(ky,kx,indexing='ij')
    os.chdir(tmp)
    Imomega_avg = np.fmax(Imomega_avg,0.0)
    return ky, Reomega_avg, Imomega_avg
    
if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']


    fig, axes =plt.subplots(2,sharex=True)
    for d in dirs:
        ky, omega, gamma = omega_from_dir(d)
        axes[0].plot(ky, omega)
        axes[1].plot(ky, gamma)
        axes[1].set_xlabel(r"$k_y \rho$")
        axes[0].set_ylabel(r"$\omega$")
        axes[1].set_ylabel(r"$\gamma$")
    axes[1].set_ylim(bottom=0.0)
    axes[0].legend(dirs)    
    plt.show()
