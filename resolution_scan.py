from stella_input import Stella_input
from slurm_utils import sbatchInDir

params = {}
params["nzed"] = [168,208]
params["nfield_periods"] = [8,10]
params["nvgrid"] = [48,60]
params["nmu"] = [32,40]

base = Stella_input(".")

for key in params:
    group = base.get_group(key)
    for value in params[key]:
        new_dir = key + str(value)
        new = base.copy(new_dir)
        new.changevar(group,key,value)
        sbatchInDir(new_dir)
