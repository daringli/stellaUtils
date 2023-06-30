from enum import Enum

# class syntax

class NetCDF(Enum):
    netCDF4 = 0
    scipy = 1

try:
    from netCDF4 import Dataset
except ImportError as e:
    try:
        from scipy.io import netcdf_file as scipy_netcdf_file
    except ImportError as e:
        # this error is fatal
        raise e
    else:
        netcdf = NetCDF.scipy
else:
    netcdf = NetCDF.netCDF4

def netcdf_file(filename, mode='r', mmap=None, format ='NETCDF3_CLASSIC', maskandscale=False):
    if netcdf == NetCDF.scipy:
        if format == 'NETCDF3_CLASSIC':
            version = 1
        elif format == 'NETCDF3_64BIT_OFFSET':
            version = 2
        else:
            raise ValueError('scipy.io netcdf_file only supports "NETCDF3_CLASSIC" or "NETCDF3_64BIT_OFFSET" formats.')
        ret = scipy_netcdf_file(filename, mode, mmap, version, maskandscale)

    elif netcdf == NetCDF.netCDF4:
        ret =  Dataset(filename, mode, format, maskandscale = maskandscale)
        ret.set_always_mask(False)
    return ret
