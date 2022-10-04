#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from ql_fluxes import get_latest_omega, read_omega
from plot_omega import omega_from_dir
    
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
        print(ky.shape)
        print(gamma.shape)
        
        axes[0].plot(ky, gamma)
        axes[1].plot(ky, gamma/ky**2)
        axes[-1].set_xlabel(r"$k_y \rho$")
        axes[0].set_ylabel(r"$\gamma$")
        axes[1].set_ylabel(r"$\gamma/k_y^2$")
    axes[1].set_ylim(bottom=0.0)
    axes[0].set_ylim(bottom=0.0)
    axes[0].legend(dirs)    
    plt.show()
