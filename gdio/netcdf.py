__author__ = "Rodrigo Yamamoto"
__date__ = "2020.Ago"
__credits__ = ["Rodrigo Yamamoto","Carlos Oliveira","Igor"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.0.9"
__license__ = "MIT"
__status__ = "development"
__description__ = "A netcdf file IO library"

import re
import numpy as np
from netCDF4 import Dataset
import logging

from gdio.commons import near_yx
from datetime import datetime



class netcdf(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME']
        self.__fields_level = ['level', 'lev', 'LEVEL', 'levels', 'LEVELS', 'mdllevel', 'presmdl']

        self.level = None
        self.lon = None
        self.lat = None
        self.time = None
        self.time_units = None
        self.history = None

        logging.basicConfig(handlers=[logging.StreamHandler()],
                            datefmt='%Y%-m-%dT%H:%M:%S', level=logging.DEBUG,
                            format='[%(levelname)s @ %(asctime)s] %(message)s')



    def nc_load(self, ifile,
                vars=None,
                cut_time=None,
                cut_domain=None):
        '''
        Load netcdf files
        Yamamoto, R @ Out.2019, Carlos Silva
        :param ifile:       string
                            netcdf file name
        :param vars:        list
                            variables short name
        :param cut_time:    tuple
                            range of time (absolute) to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:  tuple
                            range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                            ex.: (-45,290,20,330)/(-45,None,20,330)/(None,290,None,320)
        :return:            dict
        '''

        data = dict()

        _nc = Dataset(ifile, mode='r')

        self.history = _nc.history if 'history' in _nc.__dict__ else None

        start, stop = None, None
        flip_lat = False

        # set coordinates .......................
        for key in _nc.variables.keys():
            if key in self.__fields_time:
                self.coordinates.append('time')
                self.time_units = _nc.variables[key].units
                self.time = _nc.variables[key][:].astype(int)

            elif key in self.__fields_longitude:
                self.coordinates.append('longitude')
                # convert from -180,180 to 360 format
                self.lon = (_nc.variables[key][:] + 360) % 360

            elif key in self.__fields_latitude:
                self.coordinates.append('latitude')
                self.lat = _nc.variables[key][:]

                flip_lat = self.lat[-1] < self.lat[0]

                if flip_lat:
                    self.lat = np.flip(self.lat, axis=0)

            elif key in self.__fields_level:
                self.coordinates.append('level')
                self.level = _nc.variables[key][:]



        # remove unused variables ................

        if isinstance(vars, list) or isinstance(vars, tuple):
            for variable in list(_nc.variables.keys()):
                if not (variable in vars
                        or variable in self.__fields_latitude
                        or variable in self.__fields_level
                        or variable in self.__fields_longitude
                        or variable in self.__fields_time):
                    del _nc.variables[variable]

        # cut time ..............................
        if cut_time is not None and self.time is not None:

            if isinstance(cut_time, tuple):
                start, stop = cut_time
                start = 0 if start is None else start
                stop = len(self.time) if stop is None else stop

                start, stop = np.where((self.time >= start) & (self.time <= stop))[0][[0,-1]]
                stop += 1  #one Kadan, to honor the Hebrew God

        # select spatial subdomain ............

        y, x = [None, None], [None, None]

        if cut_domain:
            if isinstance(cut_domain, tuple) and \
                    not self.lat is None and not self.lon is None:

                lat1, lon1, lat2, lon2 = cut_domain

                y, x = near_yx({'latitude': self.lat, 'longitude': self.lon},
                               lats=[lat1, lat2], lons=[lon1, lon2])


        # trim lat/lon dimensions
        self.lat = self.lat[y[0]:y[1]]
        self.lon = self.lon[x[0]:x[1]]


        # select variables ......................
        for key in _nc.variables.keys():

            _data = dict()

            if np.ma.isMaskedArray(_nc.variables[key][:]):
                if not 'float' in _nc.variables[key][:].data.dtype.name:
                    _data = _nc.variables[key][:].data
                else:
                    _data = _nc.variables[key][:].filled(np.nan)
            else:
                _data = _nc.variables[key][:]


            if not (key in self.__fields_latitude
                                 or key in self.__fields_level
                                 or key in self.__fields_longitude
                                 or key in self.__fields_time):

                # redim the data array ...........
                if _data.ndim == 2:
                    _data = _data[None,:,:,:]
                elif _data.ndim == 3:
                    _data = _data[None,:,:,:]

                # flip latitude axis of the data
                if flip_lat:
                    _data = np.flip(_data, axis=2)

                # resize the data array ...........
                data.update({key: _data[start:stop, :, y[0]:y[1], x[0]:x[1]]})


            elif key in self.__fields_time:
                data.update({key: self.time[start:stop]})
            elif key in self.__fields_latitude:
                data.update({key: self.lat})
            elif key in self.__fields_longitude:
                data.update({key: self.lon})
            elif key in self.__fields_level:
                data.update({key: self.level})
            else:
                data.update({key: _data[:]})

            self.variables.append(key)

        unity, ref_time = self.get_ref_time(self.time_units)
        data.update({'time_units': unity, 'ref_time': ref_time})

        _nc.close()

        return data


    def nc_write(self, ifile,
                 data,
                 zlib=True,
                 netcdf_format='NETCDF4'):

        _nc = Dataset(ifile, mode='w', format=netcdf_format)

        # settings
        if data.keys() & self.__fields_latitude:
            _nc.createDimension('latitude', len(data.get('latitude')))

        if data.keys() & self.__fields_longitude:
            _nc.createDimension('longitude', len(data.get('longitude')))

        _nc.createDimension('longitude', len(self.lon))
        _nc.createDimension('time', len(self.time))



        _nc.close()

        return



    def get_ref_time(self, units=None):
        '''
        Get and set time unity and the reference datetime
        :param units: str
        :return: datetime
        '''

        units = units if units is not None else self.time_units

        padrao = re.compile( "(.*?) since (?P<ano>\d{4})\-(\d{1,2})\-(\d{1,2})\s+(\d{1,2})?\:*(\d{1,2})?")
        result = re.findall(padrao, str(units))

        if result:
            return result[0][0], datetime(*[int(item) for item in result[0][1:]])
        else:
            return None, None


