from __future__ import division, print_function

import subprocess
import os
from shutil import copy

import numpy as np
from glob import glob

# stdout file is out.<job ID>
stdout_prefix = "out."
stderr_prefix = "err."
jobfile = "job"
inputfile = "stella.in"
username = "sbul"

def getLatestJobIDInDir(dirname):
    """Get the ID of the latest job that is running or has been run in a directory.
    returns None if no jobs are found."""

    # see if a job is currently running
    maybeJobID = jobQueuedInDir(dirname)
    if maybeJobID is not None:
        return maybeJobID

    # look for past jobs
    #    filenames = sorted(glob(dirname + "/out.*")) # list of .out files.
    tmp = subprocess.run("ls " + dirname + " | grep out.", shell=True, capture_output=True, text=True) # list of .out files.
    filenames = tmp.stdout.split()
    if len(filenames) > 0:
        ids = [int(fn.rsplit(".",1)[-1]) for fn in filenames]
        ids.sort()
        jobID = str(ids[-1])
    else:
        jobID = None
    return jobID

def getQueuedJobDirs(u=username):
    dirs = []
    out = squeue(u=u)
    tmp = out.split('\n')
    jobs = [t.split()[0] for t in tmp[1:-1]]
    for job in jobs:
        command2 = 'scontrol show job ' + job + ' | grep WorkDir'
        result2 = subprocess.run(command2, shell=True, capture_output=True, text=True)
        dirs.append(result2.stdout.split("=",1)[1][:-1])
    return dirs,jobs


def jobQueuedInDir(d,u=username):
    d = os.path.abspath(d)
    dirs,jobs = getQueuedJobDirs(u=u)
    for (dir,job) in zip(dirs,jobs):
        if os.path.samefile(d,dir):
            return job
    else:
        return None
    

class ErrOut(object):
    """An object to parse stdout and stderr produced by a Stella run. The main purpose is to deduce the "status" of that run.

    Usage:
    ErrOut(jobID).status()"""

    donestring = "ELAPSED TIME" #presence in stdout indicates job finished.
    oomstring = "oom-kill"
    maybeoomstring = "killed"
    timestring = "time limit" # presence in stderr indicates time limit killed
    vmecfail = "ARNORM OR AZNORM EQUAL ZERO IN BCOVAR" # vmec fails and stella goes on
    
    def __init__(self,dirname,jobID=None):
        self.dirname = dirname
        if jobID is None:
            self.jobID = getLatestJobIDInDir(dirname)
        else:
            self.jobID = jobID
        
        if self.jobID is not None:
            if not os.path.isfile(self.stdout):
                raise ValueError("Output does not exist")
    
    @property
    def jobID(self):
        return self._jobID

    @jobID.setter
    def jobID(self,jobID):
        self._jobID = jobID
        if jobID is not None:
            self.stdout = self.dirname + "/" + stdout_prefix + str(self.jobID)
            self.stderr = self.dirname + "/" + stderr_prefix +  str(self.jobID)
        else:
            self.stdout = None
            self.stderr = None


    @property
    def status(self):
        if self.stdout is not None and os.path.isfile(self.stdout):
            with open(self.stdout,'r') as f:
                out = f.read()
            with open(self.stderr,'r') as f:
                err = f.read()
            if ErrOut.donestring in out:
                status = "DONE"
            elif ErrOut.vmecfail in out:
                status = "VMECFAIL"
            elif (len(err) == 0) and (len(out) > 0):
                status= "RUNNING"
            elif ErrOut.oomstring in err.lower():
                status = "OOM"
            elif ErrOut.timestring in err.lower():
                status = "TIME"
            elif ErrOut.maybeoomstring in err.lower():
                # sometimes OOM results in a different error message
                status = "OOM"
            elif (len(out) > 0):
                status= "RUNNING"
            else:
                print("somethings odd. status None for stdout:")
                print(self.stdout)
                status = None
        elif jobQueuedInDir(self.dirname) is not None:
            status = "QUEUED"
        else:
            print("somethings odd. no queued job yet no output.")
            status = None
        return status



def squeue(u="sbul",j=None):
    if j is not None:
        tmp = subprocess.run(["squeue","-j",j], capture_output=True, text=True)
    else:
        tmp = subprocess.run(["squeue","-u",u], capture_output=True, text=True)
    return tmp.stdout

def scancel(u="sbul",j=None):
    if j is not None:
        tmp = subprocess.run(["scancel","-j",j], capture_output=True, text=True)
    else:
        tmp = subprocess.run(["scancel","-u",u], capture_output=True, text=True)
    return tmp.stdout
    

def sbatchInDir(dirname, jobfile = jobfile):
    """ Launches a job inside the named directory.
    dirname -- string with path to directory to scan in
    returns: a job object, containing jobid and dirname and functions to check status.
    """
    tmp = os.getcwd()
    os.chdir(dirname)
    # this breaks backwards compatibility, python3.7 and higher only!
    finished = subprocess.run(["sbatch", jobfile], capture_output=True, text=True)
    stdout = finished.stdout
    search = "Submitted batch job "
    i1 = stdout.find(search) + len(search)
    #print(stdout[i1:])
    if i1 == -1:
        raise ValueError("Couldn't submit job. Does the jobfile '" + dirname + "/" + jobfile +"' exist?")
    jobid = stdout[i1:].split(None,1)[0]
    os.chdir(tmp)    
    return jobid


def copyInput(dirname,newdirname):
    if not os.path.exists(newdirname):
        os.makedirs(newdirname)
        copy(dirname + "/" + inputfile, newdirname)
        copy(dirname + "/" + jobfile, newdirname)
