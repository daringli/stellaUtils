#!/usr/bin/env python

from netcdf_util import netcdf_file
import numpy as np
#from collections import defaultdict

class GeometryFile(object):
    # the names should be the same as in the netcdf file
    extrapol_fudge = 1e-10
    attribs = ['bmag', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']
    
    header_to_name = {
        'bdot_grad_z': 'b_dot_grad_z',
        'qinp': 'q',
        'alpha': 'alpha_grid',
    }

    def get_name(header):
        if header in GeometryFile.header_to_name:
            ret = GeometryFile.header_to_name[header]
        else:
            ret = header
        return ret
    
    def __init__(self, geomfile = 'stella.geometry', geomtype = None):
        self.geomfile = geomfile
        gsplit = self.geomfile.rsplit('.',1)

        if geomtype is None:
            if len(gsplit) > 1:
                file_ending = gsplit[-1]
                if file_ending == 'nc':
                    self.geomtype = 2
                else:
                    self.geomtype = 1
            else:
                self.geomtype = 1
        else:
            self.geomtype = geomtype
            
        self.read_geometry()

        
    def read_geometry(self):
        if self.geomtype == 1:
            self.read_geometry_text()
        elif self.geomtype == 2:
            self.read_geometry_nc()


    def read_geometry_nc(self):
        with netcdf_file(self.geomfile) as f:
            self.zed =    f.variables['zed'][()]
            self.bmag =     f.variables['bmag'][()]

            self.Nalpha = self.bmag.shape[1]
            self.Nz = self.bmag.shape[0]
            
            if 'alpha_grid' in f.variables:
                self.alpha_grid =    f.variables['alpha_grid'][()] # add
                self.alpha0 = self.alpha_grid[0]
            else:
                # alpha = [(alpha0 + ((j - 1) * 2 * pi) / nalpha, j=1, nalpha)]
                self.alpha0 = 0.0
                self.alpha_grid = self.alpha0 + np.arange(self.Nalpha) * 2 * np.pi/self.Nalpha
                
            self.b_dot_grad_z =  f.variables['b_dot_grad_z'][()]
            self.grho  =    f.variables['grho'][()]
            self.gbdrift =  f.variables['gbdrift'][()]
            self.gbdrift0 = f.variables['gbdrift0'][()]
            self.cvdrift =  f.variables['cvdrift'][()]
            self.cvdrift0 = f.variables['cvdrift0'][()]
            self.gds2 =     f.variables['gds2'][()]
            self.gds21 =    f.variables['gds21'][()]
            self.gds22 =    f.variables['gds22'][()]
            
            if 'gds23' in f.variables:
                self.gds23 =    f.variables['gds23'][()]
            else:
                self.gds23 = np.zeros_like(self.bmag) * np.nan

            if 'gds24' in f.variables:
                self.gds24 =    f.variables['gds24'][()]
            else:
                self.gds24 = np.zeros_like(self.bmag) * np.nan

                
            if 'btor' in f.variables:
                self.btor =    f.variables['btor'][()]
            else:
                self.btor =    np.zeros_like(self.bmag)

            if 'bmag_psi0' in f.variables:
                self.bmag_psi0 =    f.variables['bmag_psi0'][()]
            else:
                self.bmag_psi0 =    self.bmag



                
            # scalars
            #self.kxfac =    f.variables['kxfac'][()]
            self.q =        f.variables['q'][()]
            if 'scale' in f.variables:
                self.scale =    f.variables['scale'][()]
            else:
                # this needs to be calculated using weird methods
                self.scale = np.nan #0.5
            #self.Rmaj =     f.variables['Rmaj'][()]
            self.shat =     f.variables['shat'][()]
            self.drhodpsi = f.variables['drhodpsi'][()]
            
            if 'dxdXcoord' in f.variables:
                self.dxdXcoord =    f.variables['dxdXcoord'][()]
            else:
                self.dxdXcoord =    self.drhodpsi

            if 'dydalpha' in f.variables:
                self.dydalpha =    f.variables['dydalpha'][()]
            else:
                self.dydalpha =    1/self.drhodpsi

            if 'rhoc' in f.variables:
                self.rhoc =    f.variables['rhoc'][()]
            else:
                self.rhoc =    self.dydalpha

            if 'rhotor' in f.variables:
                self.rhotor =    f.variables['rhotor'][()]
            else:
                self.rhotor =    self.rhoc
                
            if 'exb_nonlin' in f.variables:
                self.exb_nonlin =    f.variables['exb_nonlin'][()]
            else:
                self.exb_nonlin =    0.5 * self.dydalpha * self.dxdXcoord

            if 'exb_nonlin_p' in f.variables:
                self.exb_nonlin_p =    f.variables['exb_nonlin_p'][()]
            else:
                self.exb_nonlin_p =    0.0

            # metadata
            if 'bref' in f.variables:
                self.bref =    f.variables['bref'][()]
            else:
                self.bref = 1.0

            if 'aref' in f.variables:
                self.aref =    f.variables['aref'][()]
            else:
                self.aref = 1.0

            if 'vmec_file' in f.variables:
                self.vmec_file = f.variables['vmec_file'][()]#.filename # add
            else:
                self.vmec_file = ""

            if 'geo_file' in f.variables:
                self.geo_file = f.variables['geo_file'][()]#.filename # add
            else:
                self.geo_file = ""

            # gds23       gds24     gbdrift     cvdrift    gbdrift0   bmag_psi0        btor
                                

        
        #self.Nalpha = len(self.alpha)
        self.Nz = len(self.zed)
        self.zeta = self.zed/self.scale
    
            
    def read_geometry_text(self):
        # states
        # 0: look for header for scalars
        # 1: look for scalars
        # 2: look for header for arrays
        # 3: look for arrays
        state = 0 
        with open(self.geomfile, 'r') as f:
            lines = f.readlines()
        for l in lines:
            L = l.strip()
            if len(L) == 0:
                continue

            if state == 0:
                if L[0] == '#':
                    if 'rhoc' in L:
                        header = L[1:].split()
                        state = 1
                        Nh = len(header)
                        
            elif state == 1:
                if L[0] == '#':
                    Ls = L[1:].split()
                    if len(Ls) == Nh:
                        for h,val in zip(header,Ls):
                            name = GeometryFile.get_name(h)
                            self.__dict__[name] = float(val)
                        state = 2
                        
            elif state == 2:
                if L[0] == '#':
                    Ls = L.split()
                    if 'alpha' in L:
                        header2 = L[1:].split()
                        Nh2 = len(header2)
                        for h in header2:
                            name = GeometryFile.get_name(h)
                            self.__dict__[name] = []
                        state = 3
                        
            elif state == 3:
                Ls = L.split()
                if len(Ls) == Nh2:
                    for h,val in zip(header2,Ls):
                        name = GeometryFile.get_name(h)
                        self.__dict__[name].append(float(val))

        if state == 3:
            # if we got through the file normally
            # turn the 'header2' quantities into arrays
            
            self.Nalpha = len(set(self.alpha_grid))
            self.Nz = len(self.alpha_grid)//self.Nalpha
            self.scale = 1/np.round(self.zeta[0]/self.zed[0])
            for h in header2:
                name = GeometryFile.get_name(h)
                self.__dict__[name] = np.array(self.__dict__[name]).reshape((self.Nalpha,self.Nz)).transpose()
            self.alpha_grid = self.alpha_grid[0,:]
            self.zed = self.zed[:,0]
            self.zeta = self.zeta[:,0]
            self.grho = np.sqrt(self.gds22/self.shat**2)
            self.cvdrift0 = self.gbdrift0
            self.drhodpsi = self.dxdXcoord
            self.vmec_file = ""
            self.geo_file = ""
                
    @property
    def gradpar(self):
        return np.mean(self.b_dot_grad_z,axis=0)

    # aliases
    @property
    def bdot_grad_z(self):
        return self.b_dot_grad_z

    @bdot_grad_z.setter
    def bdot_grad_z(self,val):
        self.b_dot_grad_z = val

    @property
    def qinp(self):
        return self.q

    @qinp.setter
    def qinp(self,val):
        self.q = val
    
    # gx compatibility
    @property
    def kxfac(self):
        return self.exb_nonlin * 2.0

    @kxfac.setter
    def kxfac(self,val):
        self.exb_nonlin  = val/2.0

    def avg_attrib(self, attrib):
        attr = getattr(self, attrib)
        return np.trapz(attr/self.b_dot_grad_z, self.zed,axis=0)/np.trapz(1/self.b_dot_grad_z,self.zed,axis=0)

    def avg_bmag(self):
        return self.avg_attrib('bmag')

    def avg_gbdrift(self):
        return self.avg_attrib('gbdrift')

    def avg_gbdrift0(self):
        return self.avg_attrib('gbdrift0')
    
    def avg_cvdrift(self):
        return self.avg_attrib('cvdrift')
    
    def avg_cvdrift0(self):
        return self.avg_attrib('cvdrift0')

    def avg_gds2(self):
        return self.avg_attrib('gds2')

    def avg_gds21(self):
        return self.avg_attrib('gds21')

    def avg_gds22(self):
        return self.avg_attrib('gds22')
    
    
    def attrib_interpolator(self, attrib, unscaled_z=True):
        if unscaled_z:
            x = list(self.zed)
        else:
            x = list(self.zeta)
        # to allow for very slight extrapolation
        # so that we can copy from files with the same z range
        x[0] =   x[0] * (1 + GeometryFile.extrapol_fudge)
        x[-1] = x[-1] * (1 + GeometryFile.extrapol_fudge)
        ret = interp1d(x, getattr(self,attrib))
        return ret
        
    def copy_attrib(self, copyFromGeoFile, attrib, unscaled_z=False):
        tmp = copyFromGeoFile.attrib_interpolator(attrib, unscaled_z)
        if unscaled_z:
            setattr(self,attrib,tmp(self.zed))
        else:
            setattr(self,attrib,tmp(self.zeta))

    
    def copy_bmag(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'bmag', scaled_z)

    def copy_gbdrift(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'gbdrift', scaled_z)

    def copy_gbdrift0(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'gbdrift0', scaled_z)

    def copy_cvdrift(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'cvdrift', scaled_z)

    def copy_cvdrift0(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'cvdrift0', scaled_z)

    def copy_gds2(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'gds2', scaled_z)


    def copy_gds21(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'gds21', scaled_z)

    def copy_gds22(self, copyFromGeoFile, scaled_z=True):
        self.copy_attrib(copyFromGeoFile, 'gds22', scaled_z)


        
    def __str__(self):
        # below is the corresponding fortran code we should mimic
        # (except we want more digits)
        # write (geometry_unit,*) is a newline:  https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap02/write-1.html
        # (logic: each write is a newline, so an empty write is an empty newline)
        #      do ia = 1, Nalpha
        #         do iz = -nzgrid, nzgrid
        #            write (geometry_unit, '(15e12.4)') alpha(ia), zed(iz), zeta(ia, iz), bmag(ia, iz), b_dot_grad_z(ia, iz), &
        #               gds2(ia, iz), gds21(ia, iz), gds22(ia, iz), gds23(ia, iz), &
        #               gds24(ia, iz), gbdrift(ia, iz), cvdrift(ia, iz), gbdrift0(ia, iz), &
        #               bmag_psi0(ia, iz), btor(iz)
        #         end do
        #         write (geometry_unit, *) 
        #      end do

        def header_format(s):
            #    0.5000E+00
            return "{:>14.4E}".format(s)

        def table_format(s):
            #  5.0000e-01
            return "{:>12.4e}".format(s)
        
        ret = "#          rhoc          qinp          shat        rhotor          aref          bref     dxdXcoord      dydalpha    exb_nonlin  exb_nonlin_p\n"
        ret = ret + "#" + header_format(self.rhoc) + header_format(self.q) + header_format(self.shat) + header_format(self.rhotor) + header_format(self.aref) + header_format(self.bref) + header_format(self.dxdXcoord) + header_format(self.dydalpha) + header_format(self.exb_nonlin) + header_format(self.exb_nonlin_p)
        ret = ret + "\n\n"
        ret = ret + "     # alpha         zed        zeta        bmag bdot_grad_z        gds2       gds21       gds22       gds23       gds24     gbdrift     cvdrift    gbdrift0   bmag_psi0        btor\n"

        for i in range(self.Nalpha):
            for j in range(self.Nz):
                ret = ret + table_format(self.alpha_grid[i]) + table_format(self.zed[j]) + table_format(self.zeta[j]) + table_format(self.bmag[j,i]) + table_format(self.b_dot_grad_z[j,i]) + table_format(self.gds2[j,i]) + table_format(self.gds21[j,i]) + table_format(self.gds22[j,i]) + table_format(self.gds23[j,i]) + table_format(self.gds24[j,i]) + table_format(self.gbdrift[j,i]) + table_format(self.cvdrift[j,i]) + table_format(self.gbdrift0[j,i]) + table_format(self.bmag_psi0[j,i]) + table_format(self.btor[j,i]) + "\n"
                
        return ret


    def _write_wout(self,f,variable):
        f.createVariable(variable, 'i', ())
        var = f.variables[variable]
        #var.filename = self.wout
        
    def _write_scalar(self,f,variable):
        f.createVariable(variable, 'f', ())
        var = f.variables[variable]
        var.assignValue(getattr(self,variable))

    def _write_zed_array(self,f,variable):
        f.createVariable(variable, 'f', ('zed',))
        var = f.variables[variable]
        var[:] = getattr(self,variable)

    def _write_alpha_array(self,f,variable):
        f.createVariable(variable, 'f', ('alpha',))
        var = f.variables[variable]
        var[:] = getattr(self,variable)

        
    def _write_array(self,f,variable):
        f.createVariable(variable, 'f', ('zed','alpha',))
        var = f.variables[variable]
        var[:] = getattr(self,variable)

    def _write_string(self,f,variable):
        f.createVariable(variable, 'S1', ('char200',))
        var = f.variables[variable]
        var[:] = getattr(self,variable)

        
    def write_netcdf_output(self, output_filename):
        with netcdf_file(output_filename,'w', maskandscale=False) as f:
            f.createDimension('zed', self.Nz)
            f.createDimension('alpha', self.Nalpha)
            f.createDimension('tube', 1) # not sure about this
            f.createDimension('char200', 200)
            
            self._write_alpha_array(f,'alpha_grid')
            self._write_zed_array(f,'zed')
            self._write_zed_array(f,'zeta')
            self._write_array(f,'bmag')
            self._write_array(f,'b_dot_grad_z')
            self._write_array(f,'grho')
            self._write_array(f,'gbdrift')
            self._write_array(f,'gbdrift0')
            self._write_array(f,'cvdrift')
            self._write_array(f,'cvdrift0')
            self._write_array(f,'gds2')
            self._write_array(f,'gds21')
            self._write_array(f,'gds22')
            self._write_array(f,'gds23')
            self._write_array(f,'gds24')
            self._write_array(f,'btor')
            self._write_array(f,'bmag_psi0')

            self._write_scalar(f,'q')
            self._write_scalar(f,'scale')
            self._write_scalar(f,'shat')
            self._write_scalar(f,'drhodpsi')
            self._write_scalar(f,'dxdXcoord')
            self._write_scalar(f,'dydalpha')
            self._write_scalar(f,'rhoc')
            self._write_scalar(f,'rhotor')
            self._write_scalar(f,'exb_nonlin')
            self._write_scalar(f,'exb_nonlin_p')
            self._write_scalar(f,'aref')
            self._write_scalar(f,'bref')
            
            self._write_string(f,'vmec_file')
            self._write_string(f,'geo_file')

            
            
    def write_output(self, output_filename, netcdf = None):
        if netcdf is None:
            if self.geomtype == 2:
                netcdf = True
            else:
                netcdf = False
        if netcdf:
            self.write_netcdf_output(output_filename)
        else:
            with open(output_filename,'w') as f:
                f.write(str(self))
                    
        

        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("geofiles",nargs="+", metavar='geofile')

    parser.add_argument("-o", "--output", action="store", help="Output files to write to. If none, will print to stdout.", dest='outputs', nargs="+", default =[])
    
    parser.add_argument("-z","--zero",action='store', nargs="+", choices=['bmag','gbdrift','gbdrift0','cvdrift','cvdrift0','gds2','gds21','gds22'], help="List of geometric coefficients to set to zero.", metavar='geometric coefficient',dest='zeros', default=[])

    parser.add_argument("-a","--avg",action='store', nargs="+", choices=['bmag','gbdrift','gbdrift0','cvdrift','cvdrift0','gds2','gds21','gds22'], help="List of geometric coefficients to set to its average.", metavar='geometric coefficient',dest='avgs', default=[])
    
    parser.add_argument("-t", "--text", action="store_true", default=False, help="Whether to write output as plain text file, as opposed to netcdf. Only used if --output is specified. Default False", dest='text')

    
    parser.add_argument("--use-unscaled", "--use-unscaled-z", action='store_true', help="Whether to use unscaled zed for interpolation when copying quantities from a different geo file.", dest='unscaled', default=False)


    for attrib in GeometryFile.attribs:
        parser.add_argument("--" + attrib, action='store', help="An geofile to copy " + attrib +" from.", metavar='geofile',dest=attrib)
    
    args = parser.parse_args()

    nfilenames = len(args.geofiles)
    noutputs = len(args.outputs)
    nzeros = len(args.zeros)
    fatal_exit = False
    
    if (noutputs > 0) and (noutputs != nfilenames):
        print("Output must be unspecified or match the number of inputs.")
        exit(1)

    override_geofiles = {}
    for attrib in GeometryFile.attribs:
        inzero = attrib in args.zeros
        inavg = attrib in args.avgs
        copy_geofile = getattr(args,attrib)
        override = copy_geofile is not None
        if inzero and inavg and override:
            print("Cannot both set " + attrib + " to zero AND its average AND copy it from geofile '" + copy_geofile + "'.")
            exit(1)
        elif inzero and inavg:
            print("Cannot both set " + attrib + " to zero AND its average.")
            exit(1)
        elif  inzero and override:
            print("Cannot both copy " + attrib + " from geofile '" + copy_geofile + "' and set it to zero.")
            exit(1)
        elif  inavg and override:
            print("Cannot both copy " + attrib + " from geofile '" + copy_geofile + "' and set it to its average.")
            exit(1)

        if override:
            override_geofiles[attrib] = GeometryFile(copy_geofile)
            
            

    for i,filename in enumerate(args.geofiles):
        geo = GeometryFile(filename)
        nz = len(geo.zed)
        na = len(geo.alpha_grid)
        z = np.zeros((nz,na))
        o = np.ones((nz,na))

        for attrib_to_zero in args.zeros:
            setattr(geo, attrib_to_zero, z)

        for attrib_to_avg in args.avgs:
            setattr(geo, attrib_to_avg, geo.avg_attrib(attrib_to_avg) * o)

        for attrib_to_override in override_geofiles:
            geo.copy_attrib(override_geofiles[attrib_to_override], attrib_to_override, args.unscaled)
            
        if noutputs == 0:
            print(geo)
        else:
            geo.write_output(args.outputs[i], netcdf = not args.text)
    #geo.write_output('geo.nc')
