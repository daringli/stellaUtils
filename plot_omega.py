#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from ql_fluxes import get_latest_omega, read_omega

from itertools import cycle

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])

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
        color = next(colors)
        ky, omega, gamma = omega_from_dir(d)

        print(omega.shape)
        print(gamma.shape)
        print(ky.shape)
        if len(gamma.shape) == 2:
            gamma = gamma[0]
        if len(omega.shape) == 2:
            omega = omega[:,0]

            
        axes[0].plot(ky, omega, color = color)
        axes[1].plot(ky, gamma, color = color)
        axes[1].set_xlabel(r"$k_y \rho$")
        axes[0].set_ylabel(r"$\omega/[v_{Ti}/a]$")
        axes[1].set_ylabel(r"$\gamma/[v_{Ti}/a]$")
        imax = np.argmax(gamma)
        axes[0].plot(ky[imax], omega[imax],marker='x',label='_nolegend_',color='k')
        axes[1].plot(ky[imax], gamma[imax],marker='x',label='_nolegend_',color='k')
        
    axes[1].set_ylim(bottom=0.0)
    axes[0].legend(dirs)    
    plt.show()
