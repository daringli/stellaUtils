#!/usr/bin/env python

from scipy.io import netcdf
import numpy as np


def extract_array_value(line,array_name):
    # unused
    if array_name in l:
        tmp = l.split(array_name + "(",1)[1]
        tmp = tmp.split(')',1)
        mn = tuple([int(x) for x in tmp[0].split(',')])
        val = tmp[1].split('=',1)[1].replace(',',' ').split(None,1)[0]
        return mn, float(val)

    
def get_boundary(vmec_input):
    # reads the boundary Fourier coefficients from a VMEC input file
    rbc = {}
    zbs = {}
    
    with open(vmec_input,'r') as f:
        lines = f.readlines()
    for l in lines:
        if "RBC" in l:
           tmp = l.split("RBC(",1)[1]
           tmp = tmp.split(')',1)
           mn = tuple([int(x) for x in tmp[0].split(',')])
           val = tmp[1].split('=',1)[1].replace(',',' ').split(None,1)[0]
           rbc[mn] = float(val)

        if "ZBS" in l:
           tmp = l.split("ZBS(",1)[1]
           tmp = tmp.split(')',1)
           mn = tuple([int(x) for x in tmp[0].split(',')])
           val = tmp[1].split('=',1)[1].replace(',',' ').split(None,1)[0]
           zbs[mn] = float(val)
    return rbc, zbs

def sort_func(x):
    return np.fabs(x[0])

def sort_boundary(rbc,zbs):
    values_sorted = [(rbc[key], "rbc " + str(key)) for key in rbc] + [(zbs[key], "zbs " + str(key)) for key in zbs]
    values_sorted.sort(reverse=True, key = sort_func)
    return values_sorted
    
if __name__ == "__main__":
    import sys
    
    argv  = sys.argv
    argc = len(argv)
    if argc > 1:
        vmec_input = argv[1]
    elif argc == 1:
        vmec_input = "input.vmec"
    else:
        print("Usage: ./get_relevant_modes.py [vmec_input]")
    rbc, zbs = get_boundary(vmec_input)
    sorted  =sort_boundary(rbc,zbs)
    print(sorted[:10])
