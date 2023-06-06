#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from simsopt.mhd.vmec_diagnostics import vmec_fieldlines
from simsopt.mhd import Vmec

import sys

dirname = sys.argv[1]
wout_filename  = dirname + "/wout_vmec.nc"
vmec = Vmec(wout_filename, mpi=None)

s = 0.25
alpha = 0.0
npol = 6
ntheta = 250
theta1d =  np.linspace(-npol * np.pi, npol * np.pi, num=ntheta)

f = vmec_fieldlines(vmec, s, alpha, theta1d=theta1d, phi_center=0, plot=True, show=False)

plt.savefig(dirname + "/vmec_fieldlines.png", dpi=300)

modB = f.modB

fig, ax = plt.subplots()
ax.plot(f.phi[0][0],modB[0][0])
plt.show()
