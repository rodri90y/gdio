__author__ = "Rodrigo Yamamoto"
__date__ = "2020.Set"
__credits__ = ["Rodrigo Yamamoto","Carlos Oliveira","Igor"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.1.2"
__license__ = "MIT"
__status__ = "development"
__description__ = "A netcdf file IO library"

import re
import numpy as np
from netCDF4 import Dataset
import logging

from gdio.commons import near_yx, objectify, dict_get
from datetime import datetime


class netcdf(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME']
        self.__fields_level = {
                                'isobaricInhPa': ['level','levels','lev','presmdl'],
                                'millibars': ['level','levels','lev','presmdl'],
                                'hybrid': ['mdllevel','level_hybrid'],
                                'Eta_level': ['eta'],
                                'Sigma_level': ['sigma']
                              }
        self.__fields_order = ['time','latitude','longitude']
        self.__fields_order.insert(1, sum(self.__fields_level.values(),[])[0])

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

        data = objectify()

        _nc = Dataset(ifile, mode='r')

        self.history = _nc.history if 'history' in _nc.__dict__ else None

        start, stop = None, None
        flip_lat = False
        cut_domain_roll = 0
        typLev = None

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

            elif key.lower() in sum(self.__fields_level.values(),[]):
                self.coordinates.append(key)
                levels = list(_nc.variables[key][:].astype(int))
                typLev = _nc[key].units

        # remove unused variables ................

        if isinstance(vars, list) or isinstance(vars, tuple):
            for variable in list(_nc.variables.keys()):
                if not (variable in vars
                        or variable in self.__fields_latitude
                        or variable in sum(self.__fields_level.values(),[])
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
                while True:
                    y, x = near_yx({'latitude': self.lat, 'longitude': self.lon},
                                   lats=[lat1, lat2], lons=[lon1, lon2])

                    # if x0>x1 the longitude is rolled of x0 elements
                    # in order to avoid discontinuity 360-0 of the longitude
                    if x[0] > x[1]:
                        cut_domain_roll = -x[0]
                        self.lon = np.roll(self.lon, cut_domain_roll, axis=0)
                    else:
                        break


        # trim lat/lon dimensions
        self.lat = self.lat[y[0]:y[1]]
        self.lon = self.lon[x[0]:x[1]]

        unity, ref_time = self.get_ref_time(self.time_units)
        data.update({'ref_time': ref_time, 'time_units': unity})

        # select variables ......................
        for key in _nc.variables.keys():
            _data = objectify()

            if np.ma.isMaskedArray(_nc.variables[key][:]):
                if not 'float' in _nc.variables[key][:].data.dtype.name:
                    _data = _nc.variables[key][:].data
                else:
                    _data = _nc.variables[key][:].filled(np.nan)
            else:
                _data = _nc.variables[key][:]


            # if necessary roll longitude due discontinuity 360-0 of the longitude
            _data = np.roll(_data, cut_domain_roll, axis=-1)

            if not (key in self.__fields_latitude
                 or key.lower() in sum(self.__fields_level.values(),[])
                 or key in self.__fields_longitude
                 or key in self.__fields_time):

                # redim the data array ...........
                if _data.ndim == 2:
                    _data = _data[None,None,:,:]
                elif _data.ndim == 3:
                    _data = _data[:,None,:,:]

                # flip latitude axis of the data
                if flip_lat:
                    _data = np.flip(_data, axis=2)

                # resize the data array and consolidate ...........
                data[key] = {'value': _data[start:stop, :, y[0]:y[1], x[0]:x[1]]}
                data[key].update({
                                    'param_id': None,
                                    'type_level': typLev,
                                    'level': levels,
                                    'parameter_units': _nc[key].units
                                }
                                )

            elif key in self.__fields_time:
                data.update({key: self.time[start:stop]})
            elif key in self.__fields_latitude:
                data.update({key: self.lat})
            elif key in self.__fields_longitude:
                data.update({key: self.lon})

            self.variables.append(key)

        _nc.close()

        return data


    def nc_write(self, ifile,
                 data,
                 zlib=True,
                 netcdf_format='NETCDF4'):
        '''
        Write netcdf file

        :param ifile:           string
                                file path
        :param data:            dict
                                dataset
        :param zlib:            bool
                                enable compression
        :param netcdf_format:   string
                                netcdf format: NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT
        :return:
        '''

        _nc = Dataset(ifile, mode='w', format=netcdf_format)

        # settings
        if self.history:
            _nc.history = self.history
        else:
            _nc.history = 'Created by gdio @ {date:%Y%m%d%H}'.format(date=datetime.now())

        data = data if isinstance(data, objectify) else objectify(data)
        dims = list()


        for key, val in data.items():

            if self.verbose:
                logging.debug('''writing {var}'''.format(var=key))

            if key in self.__fields_time:
                _nc.createDimension('time', len(val))
                time = _nc.createVariable('time', 'f8', ('time',), zlib=zlib)
                time.standard_name = 'time'
                time.units = "{0} since {1}".format(data.get('time_units'), data.get('ref_time'))
                time.calendar = 'standard'
                time.axis = 'T'

                time[:] = val
                dims.append(key)

            elif key in self.__fields_latitude:
                _nc.createDimension('latitude', len(data.get('latitude')))
                lat = _nc.createVariable('latitude', 'f8', ('latitude',), zlib=zlib)
                lat.standard_name = 'latitude'
                lat.long_name = 'latitude'
                lat.units = 'degrees_north'
                lat.grads_dim = 'Y'

                lat[:] = val
                dims.append(key)

            elif key in self.__fields_longitude:
                _nc.createDimension('longitude', len(val))
                lon = _nc.createVariable('longitude', 'f8', ('longitude',), zlib=zlib)
                lon.standard_name = 'longitude'
                lon.long_name = 'longitude'
                lon.units = 'degrees_east'
                lon.grads_dim = 'x'

                # convert from 360 to -180,180 format
                lon[:] = ((val[:]-180) % 360) - 180
                dims.append(key)

            else:

                if isinstance(val, dict):

                    level_id = None

                    # create multiples z axis dimensions
                    if 'type_level' in val.keys():

                        level_id = self.__fields_level.get(val.type_level, ['level'])[0]

                        if not level_id in _nc.dimensions:
                            _nc.createDimension(level_id, len(val.level))
                            level = _nc.createVariable(level_id, 'f8', (level_id,), zlib=zlib)
                            level[:] = val.level
                            level.units = val.get('type_level')
                            level.axis = 'Z'
                            level.filling = 'off'
                            dims.append(level_id)

                    # add variables
                    if 'value' in val.keys():
                        ncvar = _nc.createVariable(key, "f8",
                                           sorted(dims, key=lambda d: self.__fields_order.index(d)),
                                           zlib=True)

                        if isinstance(val.value, np.ndarray):
                            ncvar[:] = val.value
                            ncvar.missing_value = 9.999e20
                            ncvar.units = val.get('parameter_units')


        _nc.close()




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


