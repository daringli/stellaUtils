#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from ql_fluxes import get_latest_omega, read_omega
from plot_omega import omega_from_dir

from new_kperp2 import get_kperp2

if __name__ == "__main__":
    import sys
    argv = sys.argv
    argc = len(argv)
    if argc > 1:
        dirs = argv[1:]
    else:
        dirs = ['.']

    exclude_dirs = []

    tmp = os.getcwd()
    
    
    fig, axes =plt.subplots(2,sharex=True)
    for d in dirs:
        try:
            ky, omega, gamma = omega_from_dir(d)
        except FileNotFoundError:
            exclude_dirs.append(d)
            print(d)
            os.chdir(tmp)
            continue
        kperp2 = get_kperp2(d)
        if len(gamma.shape) > 1:
            if gamma.shape[1] == 1:
                gamma = gamma[:,0]
        
        axes[0].plot(ky, gamma)
        axes[1].plot(ky, gamma/kperp2)
        axes[-1].set_xlabel(r"$k_y$")
        axes[0].set_ylabel(r"$\gamma$")
        axes[1].set_ylabel(r"$\gamma/<k_{\perp}^2>$")
    axes[1].set_ylim(bottom=0.0)
    axes[0].set_ylim(bottom=0.0)

    i = 0
    N = len(dirs)
    for j in range(N):
        if dirs[i] in exclude_dirs:
            del dirs[i]
            i = i - 1
        i = i + 1
    
    axes[0].legend(dirs)    
    plt.show()
