__author__ = "Rodrigo Yamamoto"
__date__ = "2024.Dez"
__credits__ = ["Rodrigo Yamamoto"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.0.2"
__license__ = "MIT"
__status__ = "development"
__description__ = "A netcdf file IO library"

import logging
import re
from datetime import datetime

import numpy as np
import h5py

from gdio.commons import near_yx2, objectify, near_yx


class hdf(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.__fields_index = ['x', 'y', 'ix', 'iy']
        self.__fields_latitude = [
                                    'latitude', 'lat',
                                    'xlat',
                                    'XLAT_M'
                                ]
        self.__fields_longitude = [
                                    'longitude', 'lon',
                                    'xlon',
                                    'XLONG_M'
                                    ]
        self.__fields_time = ['time']
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
        self.__fields_order = ['ensemble', 'time', 'latitude', 'longitude', 'y', 'x']
        self.__fields_order = self.__fields_order[:2] + sum(self.__fields_level.values(), []) + self.__fields_order[2:]
        self.__fields_ignore = []
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

    def h5py_dataset_iterator(self, node):

        for key, item in node.items():
            if isinstance(item, h5py.Dataset):  # test for dataset
                yield (key, item)
            elif isinstance(item, h5py.Group):  # test for group (go down)
                yield from self.h5py_dataset_iterator(item)

    def hdf_load(self, ifile,
                vars=None,
                cut_time=None,
                cut_domain=None,
                level_type=None,
                rename_vars={}):
        '''
        Load netcdf files
        Rodrigo Yamamoto @ Set.2022
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
            _hf = h5py.File(ifile, mode='r')

            self.history = _hf.attrs.get('history') if 'history' in _hf.attrs else None

            start, stop = None, None
            flip_lat = False
            cut_domain_roll = 0

            # set coordinates .......................
            self.levels['surface'] = [0]


            for key, val in self.h5py_dataset_iterator(_hf):

                if key.lower() in self.__fields_time:
                    self.coordinates.append('time')
                    self.time_units = self.get_attr(val, 'units')
                    self.time = val[:].astype(int)

                elif key.lower() in self.__fields_longitude:
                    self.coordinates.append('longitude')
                    # convert from -180,180 to 360 format
                    self.lon = (val[:] + 360) % 360

                elif key.lower() in self.__fields_latitude:
                    self.coordinates.append('latitude')
                    self.lat = val[:]

                    flip_lat = self.lat[-1] < self.lat[0]
                    flip_lat = flip_lat if isinstance(flip_lat, np.bool_) else all(flip_lat)

                    if flip_lat:
                        self.lat = np.flip(self.lat, axis=0)

                elif key.lower() in sum(self.__fields_level.values(), []):
                    typeLev = [k for k, v in self.__fields_level.items() if key in v][0]
                    self.levels['surface'] = [0]
                    self.levels[typeLev] = list(val[:].astype(int))


            # cut time ..............................
            if cut_time is not None and self.time is not None:

                if isinstance(cut_time, tuple):
                    start, stop = cut_time
                    start = 0 if start is None else start
                    stop = len(self.time) if stop is None else stop

                    stop += 1  # one Kadan, to honor the Hebrew God

            # convert lat lon to 2d mesh coordinates
            try:
                if self.lat.ndim == 1 and self.lon.ndim == 1:
                    dims = (self.lon.size, self.lat.size)
                    self.lat = np.tile(self.lat, (dims[0], 1)).T
                    self.lon = np.tile(self.lon, (dims[1], 1))
            except BaseException:
                pass

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
                                self.lon = np.roll(self.lon, cut_domain_roll, axis=1)
                            else:
                                break
                        except BaseException:
                            break

            # trim lat/lon dimensions
            try:
                self.lat = self.lat[y[0]:y[1], x[0]:x[1]]
                self.lon = self.lon[y[0]:y[1], x[0]:x[1]]
            except BaseException:
                pass

            unity, ref_time = self.get_ref_time(self.time_units)
            data.update({'ref_time': ref_time, 'time_units': unity})

            # select variables ......................
            for key, val in self.h5py_dataset_iterator(_hf):
                _data = objectify()

                if not (key.lower() in self.__fields_latitude
                        or key.lower() in sum(self.__fields_level.values(), [])
                        or key.lower() in self.__fields_longitude
                        or key.lower() in self.__fields_time
                        or key.lower() in self.__fields_index
                        or key.lower() in self.__fields_ignore):

                    if vars is None or key in vars:
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
                            if _data.ndim == 1:
                                _data = _data[None, None, None, None, :]
                            elif _data.ndim == 2:
                                _data = _data[None, None, None, :, :]
                            elif _data.ndim == 3:
                                if 'time' in self.get_attr(val, 'coordinates').lower():
                                    _data = _data[None, :, None, :, :]
                                else:
                                    _data = _data[None, None, :, :, :]
                            elif _data.ndim == 4:
                                _data = _data[None, :, :, :, :]

                            # axis translation (lon,lat > lat,lon)
                            __coord = self.get_attr(val, 'coordinates', '').split()

                            if 'lat' in __coord and 'lon' in __coord:
                                if __coord.index('lat') > __coord.index('lon'):
                                    _data = np.moveaxis(_data, -1, -2)

                            if 'latitude' in __coord and 'longitude' in __coord:

                                if __coord.index('latitude') > __coord.index('longitude'):
                                    _data = np.moveaxis(_data, -1, -2)


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
                                'long_name': self.get_attr(val, 'LongName'),
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

            _hf.close()

        except Exception as e:
            logging.error('''gdio.hdf_load > {0}'''.format(e))

        return data

    def hdf_write(self, ifile,
                 data,
                 compress_type='gzip',
                 complevel=9,
                 least_significant_digit=None,
                 force_reg_grid=False,
                 global_atrib={}
                 ):
        '''
        Write HDF file

        :param ifile:           string
                                file path
        :param data:            dict
                                dataset
        :param compress_type:   string (default gzip)
                                type of compression: zlib, gzip, lzf
        :param complevel:       int (default 4)
                                compression level

        :param least_significant_digit: int (default None)
                                        specify the power of ten of the smallest decimal place in the data that is a
                                        reliable value that dramatically improve the compression quantizing
                                        (or truncating) the data
        :param force_reg_grid:    bool (default False)
                                  disable automatic detection of non-regular grid projection
        :param global_atrib:      dict
                                  add global metadata
        :return:
        '''

        #TODO: verificar padrão para gravar o time ref, data de referencia não aparece no ncview

        _hf = h5py.File(ifile, mode='w')

        # settings
        _hf.attrs['date_created'] = '{date:%Y%m%d%H}'.format(date=datetime.now())
        _hf.attrs['version'] = f'gdio v{__version__}'

        for k, v in global_atrib.items():
            _hf.attrs[k] = v

        if self.history:
            _hf.attrs['history'] = self.history

        data = data if isinstance(data, objectify) else objectify(data)
        dims = list()

        data = data.sort(["ref_time","time"])  #sort dictionary keys to keep the time at first position

        for key, val in data.items():

            if key in self.__fields_time:

                time = _hf.create_dataset("time", val.shape, dtype='f4',
                                          compression=compress_type, compression_opts=complevel)
                time.attrs['standard_name'] = 'time'
                time.attrs['long_name'] = 'time'
                time.attrs['units'] = "{0} since {1}".format(data.get('time_units'), data.get('ref_time'))
                time.attrs['calendar'] = 'standard'
                time.attrs['axis'] = 'T'
                time.make_scale()
                time[:] = val
                dims.append(key)
                #TODO como escrever ref_time propriamente em HDF para ler no ncview
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

                            if level_id and level_id[0] not in _hf.keys():
                                level = _hf.create_dataset(level_id[0], (len(val[typLev].level),), dtype='f4',
                                                          compression=compress_type, compression_opts=complevel)
                                level[:] = val[typLev].level
                                level.attrs['units'] = typLev
                                level.attrs['axis'] = 'z' if level_id[0] in ['level'] else 'e'
                                level.attrs['filling'] = 'off'
                                level.make_scale()

                        # non-regular grid detection
                        d2lat = np.diff(val.latitude, n=2, axis=0)
                        d2lon = np.diff(val.longitude, n=2, axis=1)

                        # nonregular indexing dimension
                        # add latitude dimension
                        if not (np.isclose(d2lon, 0.0).all() and np.isclose(d2lat, 0.0).all())\
                            and not force_reg_grid:

                            if not 'y' in _hf.keys():

                                yi = _hf.create_dataset('y', data=range(val.latitude.shape[0]),
                                                                dtype='i4',
                                                                compression=compress_type,
                                                                compression_opts=complevel)

                                yi.attrs['standard_name'] = 'projection_y_coordinate'
                                yi.attrs['long_name'] = 'Y Coordinate Of Projection'
                                yi.attrs['units'] = 'm'
                                yi.attrs['axis'] = 'Y'
                                yi.make_scale()
                                dims.append('y')

                        # add longitude dimension
                        if not (np.isclose(d2lon, 0.0).all() and np.isclose(d2lat, 0.0).all())\
                            and not force_reg_grid:

                            if not 'x' in _hf.keys():
                                xi = _hf.create_dataset('x', data=range(val.latitude.shape[1]),
                                                              dtype='i4',
                                                              compression=compress_type,
                                                              compression_opts=complevel)
                                xi.attrs['standard_name'] = 'projection_x_coordinate'
                                xi.attrs['long_name'] = 'X Coordinate Of Projection'
                                xi.attrs['units'] = 'm'
                                xi.attrs['axis'] = 'X'
                                xi.make_scale()
                                dims.append('x')

                        # add latitude dimension
                        if any(k in val.keys() for k in self.__fields_latitude):
                            if not 'lat' in _hf.keys() and not 'latitude' in _hf.keys():

                                _latitude = val.latitude

                                if 'y' in _hf.keys():
                                    lat = _hf.create_dataset('latitude', data=_latitude, dtype='f4',
                                                             compression=compress_type, compression_opts=complevel)
                                else:
                                    _latitude = _latitude[:, 0]
                                    lat = _hf.create_dataset('latitude', data=_latitude, dtype='f4',
                                                             compression=compress_type, compression_opts=complevel)
                                    lat.make_scale()
                                    dims.append('latitude')

                                lat.attrs['standard_name'] = 'latitude'
                                lat.attrs['long_name'] = 'latitude'
                                lat.attrs['units'] = 'degrees_north'
                                lat.attrs['grads_dim'] = 'Y'

                        # add longitude dimension
                        if any(k in val.keys() for k in self.__fields_longitude):
                            if not 'lon' in _hf.keys() and not 'longitude' in _hf.keys():

                                _longitude = ((val.longitude - 180) % 360) - 180

                                if 'x' in _hf.keys():
                                    lon = _hf.create_dataset('longitude', data=_longitude,
                                                             dtype='f4',
                                                             compression=compress_type, compression_opts=complevel)
                                else:
                                    _longitude = _longitude[0, :]
                                    lon = _hf.create_dataset('longitude', data=_longitude,
                                                             dtype='f4',
                                                             compression=compress_type, compression_opts=complevel)
                                    lon.make_scale()
                                    dims.append('longitude')

                                lon.attrs['standard_name'] = 'longitude'
                                lon.attrs['long_name'] = 'longitude'
                                lon.attrs['units'] = 'degrees_east'
                                lon.attrs['grads_dim'] = 'x'

                        # add variables
                        if 'value' in val[typLev].keys():

                            # set member dimension
                            if isinstance(val[typLev].value, np.ndarray):

                                if val[typLev].value.ndim > 4 and val[typLev].value.shape[0] > 1:
                                    if not 'ensemble' in _hf.keys():
                                        ens = _hf.create_dataset('ensemble', (val[typLev].value.shape[0],),
                                                         dtype='f4',
                                                        compression=compress_type, compression_opts=complevel)
                                        ens.attrs['standard_name'] = 'ensemble'
                                        ens.attrs['long_name'] = 'Ensemble member'
                                        ens.attrs['axis'] = 'e'
                                        ens.make_scale()
                                        dims.append('ensemble')

                                # variable setup
                                if 'ensemble' in _hf.keys():
                                    if len(dims + z_dims) == 5:
                                        __dat = val[typLev].value[:, :, :, :, :]  # level data
                                    elif len(dims + z_dims) == 4:
                                        __dat = val[typLev].value[:, :, 0, :, :]  # surface data
                                else:
                                    if len(dims + z_dims) == 4:
                                        __dat = val[typLev].value[0, :, :, :, :]  # level data
                                    elif len(dims + z_dims) == 3:
                                        __dat = val[typLev].value[0, :, 0, :, :]  # surface data

                                if least_significant_digit is not None:
                                    __dat = __dat.round(decimals=least_significant_digit)

                                hdfvar = _hf.create_dataset(key,
                                                           data=__dat,
                                                           compression=compress_type, compression_opts=complevel)
                                hdfvar.attrs['missing_value'] = 999999
                                hdfvar.attrs['units'] = str(val.get('parameter_units'))

                                coords = sorted(dims + z_dims, key=lambda d: self.__fields_order.index(d))

                                for i, c in enumerate(coords):
                                    hdfvar.dims[i].attach_scale(_hf[c])

                                hdfvar.attrs['coordinates'] = " ".join(coords)

        _hf.close()

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

        try:
            for dim in self.get_attr(data, 'coordinates').split():
                if dim.lower() not in self.__fields_latitude + \
                        self.__fields_longitude + \
                        self.__fields_index + \
                        self.__fields_time:
                    # convert level name to grib standard name type
                    typeLev = [k for k, v in self.__fields_level.items() if dim in v][0]
        except:
            pass

        return typeLev

    def get_attr(self, nc, attr, default=None):
        '''
        Get netcdf attribute
        :param nc:      object
                        netcdf object
        :param attr:    string
                        attribute name
        :return:        data
        '''

        __attr = nc.attrs.get(attr) if attr in nc.attrs.keys() else default
        if isinstance(__attr, bytes):
            return __attr.decode("utf-8")
        else:
            return __attr

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



    # @staticmethod
    def is_hdf(self, ifile):
        '''
        Check if is hdf file
        from Rodrigo@Set.2022
        :rtype: bool
        :return:
        '''
        if isinstance(ifile, str):
            if h5py.is_hdf5(ifile):
                return True

        return False