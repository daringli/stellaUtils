#!/usr/bin/env python

diff_method = 'forward'
rel_step = 1e-05
abs_step = 1e-08
BIG_NUM = 100.0

#import pdbl; pdb.set_trace()

wait_time = 30 # seconds


from scipy.optimize import minimize
import numpy as np
#from mpi4py import MPI
#from simsopt.util.mpi import MpiPartition, log
#from simsopt.mhd.vmec import Vmec
from simsopt.geo.surfacerzfourier import SurfaceRZFourier
from simsopt._core.util import finite_difference_steps

import os
from shutil import copy
from time import sleep

from netcdf_util import netcdf_file


from slurm_utils import ErrOut, sbatchInDir, scancel
from stella_input import Stella_input


from plot_omega import omega_from_dir
from new_kperp2 import get_kperp2

#log()


filename = "input.vmec"

surf = SurfaceRZFourier.from_vmec_input(filename)
surf.change_resolution(mpol=6, ntor=6)

# number of dofs for each iteration:
#8, 24, 48, 80
# we go with 24

#mpols = 5
#max_mode = 2

mpols = 3
max_mode = 1

surf.fix_all()
# surf.fixed_range(mmin=0, mmax=max_mode, 
#                  nmin=-max_mode, nmax=max_mode, fixed=False)
# surf.fix("rc(0,0)") # Major radius
surf.unfix("rc(0,1)")

x = surf.x

base_nml = "stripped_input.vmec"
my_iter = 0
N = len(x)




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


def extract_value(dirname):
    try:
        ky, gamma = get_gamma(dirname)
    except FileNotFoundError:
        print("failed to extract gamma from: " + dirname)
        return BIG_NUM
    kperp2 = get_kperp2(dirname)
    Nky = len(kperp2)
    return np.trapz(gamma/kperp2,ky)


    

def dirname_from_indices(i,n):
    dirname = str(i).zfill(5) + "_" + str(n).zfill(3)
    return dirname

def J(x,n = 0):
    global my_iter
    surf.x = x
    dirname = dirname_from_indices(my_iter,n)

    boundary_string = surf.get_nml()
    boundary_string = boundary_string.rsplit('\n',2)[0]
    boundary_string = boundary_string.split('\n',3)[-1] + "\n"
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        os.mkdir(dirname + "/nc")
        copy("stella.in",dirname)
        copy("job",dirname)
        new_vmec_input = dirname + "/input.vmec" 
        copy(base_nml,new_vmec_input)
        
        # modify vmec file
        with open(new_vmec_input,'r') as f:
            lines = f.readlines()
        Nlines= len(lines)
        if lines[Nlines-1].strip() == "/":
            lines.insert(Nlines-1,boundary_string)
        else:
            lines.insert(Nlines-2,boundary_string)

        with open(new_vmec_input,'w') as f:
            f.write("".join(lines))
    else:
        print("dir already exists: " + dirname)
        #raise ValueError("dir already exists: " + dirname)

def Jfull(x0):
    global my_iter
    xsteps = finite_difference_steps(x0,abs_step,rel_step)
    N = len(x0)
    dirnames = [dirname_from_indices(my_iter,idx) for idx in range(N+1)]
    
    
    J(x0)

    x = np.zeros(N)
    for idx in range(N):
        x[:] = x0
        x[idx] = x[idx] + xsteps[idx]
        J(x,idx+1)

    # submit all jobs
    for dirname in dirnames:
        status = ErrOut(dirname).status
        if status is None:
            sbatchInDir(dirname, 'job')

    # check if simulations are done and extract values
    values = np.zeros(N+1)
    indices = list(range(0,N+1))
    while len(indices) > 0:
        to_del = []
        num_del = 0
        for _i,idx in enumerate(indices):
            dirname = dirnames[idx]
            eo = ErrOut(dirname)
            status = eo.status

            print(dirname + ":" + str(status))
            if (status == "DONE"):
                values[idx] = extract_value(dirname)
                to_del.append(_i)
            elif ((status == "OOM") or (status == "TIME") or (status == "VMECFAIL")  or (status == "STELLAFAIL")):
                jobID = eo.jobID
                scancel(j=jobID)
                values[idx] = BIG_NUM
                to_del.append(_i)
        for _i in to_del:
            del indices[_i - num_del]
            num_del = num_del + 1

        if len(indices)  > 0:
            sleep(wait_time)

    if values[0] == BIG_NUM:
        # use gradients to make the optimizer avoid non-convergence
        values[1:] = np.clip(values[1:], None, BIG_NUM*0.75)
    df = values[1:] - values[0]
    dfdx = df/xsteps
    with open("log.txt",'a') as f:
        f.write(", ".join(( str(my_iter), np.array_str(x0, max_line_width=np.inf, precision=None).replace('\n', ','), str(values[0]), np.array_str(dfdx, max_line_width=np.inf, precision=None).replace('\n', ',') + "\n")))
    my_iter = my_iter + 1
    # flip the sign to maximize
    print((values[0],dfdx))
    return (values[0],dfdx)

header = ["iter"]  + surf.dof_names + ["J"] + ["dJ/d(" + name + ")" for name in surf.dof_names]

with open("log.txt",'w') as f:
    f.write(", ".join(header) + "\n")

MAXITER=500
res = minimize(Jfull, x, jac=True, options={'maxiter': MAXITER, 'maxcor': 300})
print("success = " + str(res.success))
print("x =  " + str(res.x))
print("status = " + str(res.message))
print("f = " + str(res.fun))
print("iteration = " + str(res.nit))
