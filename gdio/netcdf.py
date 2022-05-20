__author__ = "Rodrigo Yamamoto"
__date__ = "2022.Mai"
__credits__ = ["Rodrigo Yamamoto", "Carlos Oliveira", "Igor"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.2.9"
__license__ = "MIT"
__status__ = "development"
__description__ = "A netcdf file IO library"

import logging
import re
from datetime import datetime

import numpy as np
from netCDF4 import Dataset

from gdio.commons import near_yx2, objectify, near_yx


class netcdf(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'XLAT_M', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'XLONG_M', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME']
        self.__fields_members = ['member']
        self.__fields_level = {
            'isobaricInhPa': ['level', 'levels', 'lev', 'presmdl'],
            'hybrid': ['mdllevel', 'level_hybrid'],
            'etaLevel': ['eta'],
            'eta': ['eta'],
            'sigmaLevel': ['sigma'],
            'sigma': ['sigma'],
            'theta': ['theta'],
            'isentropic': ['theta'],

            'surface': [],
            'meanSea': [],
            'depthBelowLandLayer': [],
            'heightAboveGround': [],
            'highCloudTop': [],
            'middleCloudTop': [],
            'nominalTop': [],
            'lowCloudTop': [],
            'middleCloudLayer': [],
            'lowCloudLayer': [],
            'atmosphereSingleLayer': [],
            'entireAtmosphere': [],
            'entireOcean': [],
            'convectiveCloudLayer': [],
            'boundaryLayerCloudLayer': [],
            'highCloudBottom': [],
            'lowCloudBottom': [],
            'convectiveCloudBottom': [],
            'convectiveCloudTop': [],
            'middleCloudBottom': [],
            'pressureFromGroundLayer': [],
            'maxWind': [],
            'tropopause': [],
            'isothermZero': [],
            'potentialVorticity': []

        }
        self.__fields_3dlevel = ['isobaricInhPa', 'hybrid', 'sigma', 'eta']
        self.__fields_order = ['ensemble', 'time', 'latitude', 'longitude']
        self.__fields_order = self.__fields_order[:2] + sum(self.__fields_level.values(), []) + self.__fields_order[2:]

        self.lon = None
        self.lat = None
        self.time = None
        self.time_units = None
        self.history = None
        self.centre = 0
        self.grid_description = None
        self.levels = dict()

        logging.basicConfig(datefmt='%Y%-m-%dT%H:%M:%S', level=logging.DEBUG,
                            format='[%(levelname)s @ %(asctime)s] %(message)s')

    def nc_load(self, ifile,
                vars=None,
                cut_time=None,
                cut_domain=None,
                level_type=None,
                rename_vars={}):
        '''
        Load netcdf files
        Rodrigo Yamamoto @ Fev.2021, Carlos Silva
        :param ifile:       string
                            netcdf file name
        :param vars:        list
                            variables short name
        :param cut_time:    tuple
                            range of time (absolute) to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:  tuple
                            range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                            ex.: (-45,290,20,330)/(-45,None,20,330)/(None,290,None,320)
        :param level_type:  list
                            type of level (hybrid, isobaricInhPa, surface)
        :param rename_vars: dictonary
                            rename variables names (key) for a new name (value).
                            Eg. {'tmpmdl': 't', 'tmpprs': 't'}
        :return:            dictonary/attributes
                            multiple time data container
        '''

        data = objectify()

        try:
            _nc = Dataset(ifile, mode='r')

            self.history = _nc.history if 'history' in _nc.__dict__ else None

            start, stop = None, None
            flip_lat = False
            cut_domain_roll = 0

            # set coordinates .......................
            self.levels['surface'] = [0]

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

                elif key.lower() in sum(self.__fields_level.values(), []):
                    typeLev = [k for k, v in self.__fields_level.items() if key in v][0]
                    self.levels['surface'] = [0]
                    self.levels[typeLev] = list(_nc.variables[key][:].astype(int))

            # cut time ..............................
            if cut_time is not None and self.time is not None:

                if isinstance(cut_time, tuple):
                    start, stop = cut_time
                    start = 0 if start is None else start
                    stop = len(self.time) if stop is None else stop

                    stop += 1  # one Kadan, to honor the Hebrew God

            # convert lat lon to 2d mesh coordinates
            if self.lat.ndim == 1 and self.lon.ndim == 1:
                dims = (self.lon.size, self.lat.size)
                self.lat = np.tile(self.lat, (dims[0], 1)).T
                self.lon = np.tile(self.lon, (dims[1], 1))


            # select spatial subdomain ............

            y, x = [None, None], [None, None]

            if cut_domain:
                if isinstance(cut_domain, tuple) and \
                        not self.lat is None and not self.lon is None:

                    lat1, lon1, lat2, lon2 = cut_domain

                    while True:
                        y, x = near_yx2({'latitude': self.lat, 'longitude': self.lon},
                                       lats=[lat1, lat2], lons=[lon1, lon2])

                        # if x0>x1 the longitude is rolled of x0 elements
                        # in order to avoid discontinuity 360-0 of the longitude
                        try:
                            if x[0] > x[1]:
                                cut_domain_roll = -x[0]
                                self.lon = np.roll(self.lon, cut_domain_roll, axis=0)
                            else:
                                break
                        except BaseException:
                            break

            # trim lat/lon dimensions
            self.lat = self.lat[y[0]:y[1], x[0]:x[1]]
            self.lon = self.lon[y[0]:y[1], x[0]:x[1]]

            unity, ref_time = self.get_ref_time(self.time_units)
            data.update({'ref_time': ref_time, 'time_units': unity})

            # select variables ......................
            for key, val in _nc.variables.items():

                _data = objectify()

                if not (key in self.__fields_latitude
                        or key.lower() in sum(self.__fields_level.values(), [])
                        or key in self.__fields_longitude
                        or key in self.__fields_time):

                    if vars is None or key in vars:

                        if np.ma.isMaskedArray(_nc.variables[key][:]):
                            if not 'float' in val[:].data.dtype.name:
                                _data = val[:].data
                            else:
                                _data = val[:].filled(np.nan)
                        else:
                            _data = val[:]

                        # rename variables
                        for k, v in rename_vars.items():
                            if key in k:
                                key = v

                        # if necessary roll longitude due discontinuity 360-0 of the longitude
                        _data = np.roll(_data, cut_domain_roll, axis=-1)

                        typLev = self.__gettypLev(val)

                        if (level_type is None or typLev in level_type):

                            # redim the data array ...........
                            if _data.ndim == 2:
                                _data = _data[None, None, None, :, :]
                            elif _data.ndim == 3:
                                if 'time' in val.dimensions:
                                    _data = _data[None, :, None, :, :]
                                else:
                                    _data = _data[None, None, :, :, :]
                            elif _data.ndim == 4:
                                _data = _data[None, :, :, :, :]

                            # flip latitude axis of the data
                            if flip_lat:
                                _data = np.flip(_data, axis=3)

                            # resize the data array and consolidate ...........
                            __tmp = {
                                typLev: {
                                    'value': _data[:, start:stop, :, y[0]:y[1], x[0]:x[1]],
                                    'level': self.levels[typLev],
                                    'members': list(range(0, _data.shape[0]))
                                },
                                'param_id': None,
                                'long_name': self.get_attr(val, 'long_name'),
                                'parameter_units': self.get_attr(val, 'units'),
                                'latitude': self.lat,
                                'longitude': self.lon
                            }

                            if key in data.keys():
                                data[key].update(__tmp)
                                data[key].level_type.append(typLev)
                            else:
                                data[key] = __tmp
                                data[key].level_type = [typLev]
                                self.variables.append(key)


                elif key in self.__fields_time:
                    data.update({key: self.time[start:stop]})

            _nc.close()

        except Exception as e:
            logging.error('''gdio.nc_load > {0}'''.format(e))

        return data

    def __gettypLev(self, data, typeLev='surface'):
        '''
        Get type of level of dimemsion attribute
        :param data:        obj
                            netcdf object
        :param typeLev:     string
                            default level in case of level type dimension absence
        :return:            string
                            level type
        '''

        for dim in data.dimensions:
            if dim not in self.__fields_latitude + \
                    self.__fields_longitude +  \
                    self.__fields_time:
                # convert level name to grib standard name type
                typeLev = [k for k, v in self.__fields_level.items() if dim in v][0]

        return typeLev

    def nc_write(self, ifile,
                 data,
                 zlib=True,
                 netcdf_format='NETCDF4',
                 complevel=4,
                 least_significant_digit=None
                 ):
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

        :param complevel:       int (default 4)
                                compression level

        :param least_significant_digit: int (default None)
                                        specify the power of ten of the smallest decimal place in the data that is a
                                        reliable value that dramatically improve zlib compression by quantizing
                                        (or truncating) the data
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
                time = _nc.createVariable('time', 'f4', ('time',), zlib=zlib)
                time.standard_name = 'time'
                time.units = "{0} since {1}".format(data.get('time_units'), data.get('ref_time'))
                time.calendar = 'standard'
                time.axis = 'T'

                time[:] = val
                dims.append(key)

            else:

                if isinstance(val, dict):
                    z_dims = list()
                    for typLev in val.level_type:

                        if typLev in self.__fields_level:
                            level_id = self.__fields_level.get(typLev)[:1]
                        else:
                            level_id = self.__fields_level.get('surface')[:1]

                        z_dims = level_id if not z_dims == level_id else z_dims

                        # create multiples z axis dimensions
                        if 'level_type' in val.keys():
                            if level_id and level_id[0] not in _nc.dimensions:
                                _nc.createDimension(level_id[0], len(val[typLev].level))
                                level = _nc.createVariable(level_id[0], 'f4', (level_id[0],),
                                                                 zlib=zlib,
                                                                 complevel=complevel)
                                level[:] = val[typLev].level
                                level.units = typLev
                                level.axis = 'z' if level_id[0] in ['level'] else 'e'
                                level.filling = 'off'

                        # add latitude dimension
                        if any(k in val.keys() for k in self.__fields_latitude):
                            if not 'latitude' in _nc.dimensions:

                                if val.latitude.ndim > 1:
                                    _latitude = val.latitude[:,0]

                                _nc.createDimension('latitude', len(_latitude))
                                lat = _nc.createVariable('latitude', 'f4', ('latitude',),
                                                                 zlib=zlib,
                                                                 complevel=complevel)
                                lat.standard_name = 'latitude'
                                lat.long_name = 'latitude'
                                lat.units = 'degrees_north'
                                lat.grads_dim = 'Y'

                                lat[:] = _latitude
                                dims.append('latitude')

                        # add longitude dimension
                        if any(k in val.keys() for k in self.__fields_longitude):
                            if not 'longitude' in _nc.dimensions:

                                if val.longitude.ndim > 1:
                                    _longitude = val.longitude[0,:]

                                _nc.createDimension('longitude', len(_longitude))
                                lon = _nc.createVariable('longitude', 'f4', ('longitude',),
                                                                 zlib=zlib,
                                                                 complevel=complevel)
                                lon.standard_name = 'longitude'
                                lon.long_name = 'longitude'
                                lon.units = 'degrees_east'
                                lon.grads_dim = 'x'

                                # convert from 360 to -180,180 format
                                lon[:] = ((_longitude - 180) % 360) - 180

                                dims.append('longitude')

                        # add variables
                        if 'value' in val[typLev].keys():

                            # set member dimension
                            if isinstance(val[typLev].value, np.ndarray):

                                if val[typLev].value.ndim > 4 and val[typLev].value.shape[0] > 1:
                                    if not 'ensemble' in _nc.dimensions:
                                        _nc.createDimension('ensemble', val[typLev].value.shape[0])
                                        ens = _nc.createVariable('ensemble', 'f4', ('ensemble',),
                                                                 zlib=zlib,
                                                                 complevel=complevel)
                                        ens.standard_name = 'ensemble'
                                        ens.long_name = 'Ensemble member'
                                        ens.axis = 'e'
                                        dims.append('ensemble')

                                # variable setup
                                ncvar = _nc.createVariable(key, "f4",
                                                           sorted(dims + z_dims,
                                                                  key=lambda d: self.__fields_order.index(d)),
                                                           zlib=True,
                                                           complevel=complevel,
                                                           least_significant_digit=least_significant_digit)

                                if 'ensemble' in _nc.dimensions:
                                    if len(dims + z_dims) == 5:
                                        ncvar[:] = val[typLev].value[:, :, :, :, :]  # level data
                                    elif len(dims + z_dims) == 4:
                                        ncvar[:] = val[typLev].value[:, :, 0, :, :]  # surface data
                                else:
                                    if len(dims + z_dims) == 4:
                                        ncvar[:] = val[typLev].value[0, :, :, :, :]  # level data
                                    elif len(dims + z_dims) == 3:
                                        ncvar[:] = val[typLev].value[0, :, 0, :, :]  # surface data

                                ncvar.missing_value = 999999
                                ncvar.units = str(val.get('parameter_units'))

        _nc.close()

    def get_attr(self, nc, attr, default=None):
        '''
        Get netcdf attribute
        :param nc:      object
                        netcdf object
        :param attr:    string
                        attribute name
        :return:        data
        '''
        return nc.getncattr(attr) if attr in nc.ncattrs() else default

    def get_ref_time(self, units=None):
        '''
        Get and set time unity and the reference datetime
        :param units: str
        :return: datetime
        '''

        units = units if units is not None else self.time_units

        padrao = re.compile("(.*?) since (?P<ano>\\d{4})\\-(\\d{1,2})\\-(\\d{1,2})\\s+(\\d{1,2})?\\:*(\\d{1,2})?")
        result = re.findall(padrao, str(units))

        if result:
            return result[0][0], datetime(*[int(item) for item in result[0][1:]])
        else:
            return None, None
