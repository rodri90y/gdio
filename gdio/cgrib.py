import numpy as np
import pyproj
from eccodes import *


class cgrib():

    def __init__(self, gid, gribkeys=None):

        self.perturbationNumber = 0

        _gribkeys = [
            'centre',
            'dataDate', 'dataTime', 'unitOfTimeRange',
            'year', 'month', 'day', 'hour', 'minute',
            'validityDate', 'validityTime',
            'step',
            'stepUnits',
            'paramId',
            'typeOfLevel',
            'level',
            'shortName',
            'name',
            'parameterUnits'
        ]

        self.gridkeys = [
            'gridType',
            'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees',
            'latitudeOfLastGridPointInDegrees', 'longitudeOfLastGridPointInDegrees',
            'Ni', 'Nj',
            'iDirectionIncrementInDegrees', 'jDirectionIncrementInDegrees',
            'truncateDegrees', 'grib2divider',
            'angleOfRotationInDegrees', 'latitudeOfSouthernPoleInDegrees', 'longitudeOfSouthernPoleInDegrees',
            'latitudeWhereDxAndDyAreSpecifiedInDegrees', 'projectionCentreFlag', 'orientationOfTheGridInDegrees',
            'LoVInDegrees', 'LaDInDegrees', 'Latin1InDegrees', 'Latin2InDegrees',
            'DxInMetres', 'DyInMetres', 'iScansPositively', 'jScansPositively',
            'shapeOfTheEarth', 'radius', 'scaledValueOfEarthMajorAxis', 'scaledValueOfEarthMinorAxis',
            'LaD', 'Latin',
            'perturbationNumber', 'scaleFactorOfRadiusOfSphericalEarth'
        ]


        _keys = gribkeys if gribkeys else _gribkeys+self.gridkeys

        for k in _keys:

            try:
                _type = eccodes.codes_get_native_type(gid, k)

                if _type is str:
                    setattr(self, k, eccodes.codes_get_string(gid, k))
                elif _type is int:
                    try:
                        setattr(self, k, int(eccodes.codes_get_double(gid, k)))
                    except ArrayTooSmallError as e:
                        setattr(self, k, eccodes.codes_get_array(gid, k))
                elif _type is float:
                    try:
                        setattr(self, k, eccodes.codes_get_double(gid, k))
                    except ArrayTooSmallError as e:
                        setattr(self, k, eccodes.codes_get_array(gid, k))
                else:
                    setattr(self, k, eccodes.codes_get(gid, k))
            except KeyValueNotFoundError as e:
                pass

        # value set
        self.values = codes_get_values(gid)
        # change the flatten data array to actual data dimension
        if isinstance(self.values, np.ndarray):
            self.values = self.values.reshape(self.Nj, self.Ni)

        self.projparams = self._set_projparams()


    def __new__(cls, *args, **kwargs):
        instance = object.__new__(cls)
        instance.__dict__ = dict()
        return instance

    def __setattr__(self, key, value):
        object.__setattr__(self, key, value)

    def __getitem__(self, x):
        return getattr(self, x)

    def keys(self):
        return self.__dict__.keys()

    def items(self):
        return self.__dict__.items()

    def get(self, key, default=None):
        return self[key] if key in self.keys() else default

    def latlons(self):
        '''alias for L{grid}'''
        return self.grid()

    def _set_projparams(self):

        projparams = {}

        # set Earth a,b axis
        # check for radius key, if it exists just use it
        # and don't bother with shapeOfTheEarth
        if self['radius']:
            projparams['a'] = self['radius']
            projparams['b'] = self['radius']
        else:
            if self['shapeOfTheEarth'] == 6:
                projparams['a'] = 6371229.0
                projparams['b'] = 6371229.0
            elif self['shapeOfTheEarth'] in [3, 7]:
                    #simplify version, without considerate scaleFactorOfRadiusOfSphericalEarth key
                    projparams['a'] = self['scaledValueOfEarthMajorAxis']
                    projparams['b'] = self['scaledValueOfEarthMinorAxis']
            elif self['shapeOfTheEarth'] == 4:
                projparams['a'] = 6378137.0
                projparams['b'] = 6356752.314
            elif self['shapeOfTheEarth'] == 2:
                projparams['a'] = 6378160.0
                projparams['b'] = 6356775.0
            elif self['shapeOfTheEarth'] == 1:
                # simplify version, without considerate scaleFactorOfRadiusOfSphericalEarth key
                projparams['a'] = self['scaledValueOfRadiusOfSphericalEarth']
                projparams['b'] = self['scaledValueOfRadiusOfSphericalEarth']
            elif self['shapeOfTheEarth'] == 0:
                projparams['a'] = 6367470.0
                projparams['b'] = 6367470.0
            elif self['shapeOfTheEarth'] == 5:  # WGS84
                projparams['a'] = 6378137.0
                projparams['b'] = 6356752.3142
            elif self['shapeOfTheEarth'] == 8:
                projparams['a'] = 6371200.0
                projparams['b'] = 6371200.0
            else:
                raise ValueError('unknown shape of the earth flag')


        #projections parameters
        if self['gridType'] in ['reduced_gg', 'reduced_ll', 'regular_gg', 'regular_ll']:  # regular lat/lon grid
            projparams['proj'] = 'cyl'
        elif self['gridType'] == 'polar_stereographic':
            projparams['proj'] = 'stere'
            projparams['lat_ts'] = self['latitudeWhereDxAndDyAreSpecifiedInDegrees']
            # simplyfy version only
            if self['projectionCentreFlag'] == 0:
                projparams['lat_0'] = 90.
            else:
                projparams['lat_0'] = -90.
            projparams['lon_0'] = self['orientationOfTheGridInDegrees']
        elif self['gridType'] == 'lambert':
            projparams['proj'] = 'lcc'
            projparams['lon_0'] = self['LoVInDegrees']
            projparams['lat_0'] = self['LaDInDegrees']
            projparams['lat_1'] = self['Latin1InDegrees']
            projparams['lat_2'] = self['Latin2InDegrees']
        elif self['gridType'] == 'mercator':
            #simplyfy version fix scale without grib2divider,

            lon1 = self['longitudeOfFirstGridPointInDegrees']
            lon2 = self['longitudeOfLastGridPointInDegrees']

            if self['truncateDegrees']:
                lon1 = int(lon1)
                lon2 = int(lon2)

            if self.get('LaD'):
                projparams['lat_ts'] = self['LaD'] / 1000.
            else:
                projparams['lat_ts'] = self['Latin'] / 1000.

            if lon2 < lon1: lon2 += 360.  # domain crosses Greenwich
            projparams['lon_0'] = 0.5 * (lon1 + lon2)
            projparams['proj'] = 'merc'

        elif self['gridType'] in ['rotated_ll', 'rotated_gg']:

            rot_angle = self['angleOfRotationInDegrees']
            pole_lat = self['latitudeOfSouthernPoleInDegrees']
            pole_lon = self['longitudeOfSouthernPoleInDegrees']
            projparams['o_proj'] = 'longlat'
            projparams['proj'] = 'ob_tran'
            projparams['o_lat_p'] = -pole_lat
            projparams['o_lon_p'] = rot_angle
            projparams['lon_0'] = pole_lon

        else:  # unsupported grid type.
            projparams = None

        return projparams

    def grid(self):
        '''
        Generate the lat/lon grid from the grib projection parameters
        Based on the pygrib lib
        '''

        if self['gridType'] in ['regular_gg', 'regular_ll']:

            lat1 = self['latitudeOfFirstGridPointInDegrees']
            lat2 = self['latitudeOfLastGridPointInDegrees']
            lon1 = self['longitudeOfFirstGridPointInDegrees']
            lon2 = self['longitudeOfLastGridPointInDegrees']

            delon = self['iDirectionIncrementInDegrees']
            delat = self['jDirectionIncrementInDegrees']

            delat = delat if lat1<lat2 else delat*-1
            lats = np.arange(lat1, lat2 + delat, delat)
            delon = delon if lon1 < lon2 else delon * -1
            lons = np.arange(lon1, lon2 + delon, delon)

            lons, lats = np.meshgrid(lons, lats)

        elif self['gridType'] == 'mercator':

            lon1 = self['longitudeOfFirstGridPointInDegrees']
            lon2 = self['longitudeOfLastGridPointInDegrees']

            if self['truncateDegrees']:
                lon1 = int(lon1)
                lon2 = int(lon2)

            lat1 = self['latitudeOfFirstGridPointInDegrees']
            lat2 = self['latitudeOfLastGridPointInDegrees']

            if self['truncateDegrees']:
                lat1 = int(lat1)
                lat2 = int(lat2)

            pj = pyproj.Proj(self.projparams)
            llcrnrx, llcrnry = pj(lon1, lat1)
            urcrnrx, urcrnry = pj(lon2, lat2)

            nx = self['Ni']
            ny = self['Nj']
            dx = (urcrnrx - llcrnrx) / (nx - 1)
            dy = (urcrnry - llcrnry) / (ny - 1)
            x = llcrnrx + dx * np.arange(nx)
            y = llcrnry + dy * np.arange(ny)
            x, y = np.meshgrid(x, y)
            lons, lats = pj(x, y, inverse=True)

        elif self['gridType'] == 'lambert':
            lat1 = self['latitudeOfFirstGridPointInDegrees']
            lon1 = self['longitudeOfFirstGridPointInDegrees']

            nx = self['Ni']
            ny = self['Nj']
            dx = self['DxInMetres']
            dy = self['DyInMetres']
            pj = pyproj.Proj(self.projparams)
            llcrnrx, llcrnry = pj(lon1, lat1)
            # Set increment direction here for the grid.
            # NOTE: some GRIB files are arranged with first gridpoint
            # in top left, or top right corner for example...
            if self['iScansPositively'] == 0 and dx > 0: dx = -dx
            if self['jScansPositively'] == 0 and dy > 0: dy = -dy
            x = llcrnrx + dx * np.arange(nx)
            y = llcrnry + dy * np.arange(ny)
            x, y = np.meshgrid(x, y)
            lons, lats = pj(x, y, inverse=True)

        elif self['gridType'] in ['rotated_ll', 'rotated_gg']:
            # rotatedlats = self['distinctLatitudes']
            # rotatedlons = self['distinctLongitudes']
            nx = self['Ni']
            ny = self['Nj']
            # d2r = np.pi/180.
            # lonsr, latsr = np.meshgrid(rotatedlons*d2r, rotatedlats*d2r)
            pj = pyproj.Proj(self.projparams)
            # lons,lats = pj(lonsr,latsr,inverse=True)
            # lons = self['longitudes'].reshape(ny, nx)
            # lats = self['latitudes'].reshape(ny, nx)
            raise ValueError('unsupported grid {0}'.format(self['gridType']))
        elif self['gridType'] == 'polar_stereographic':

            lat1 = self['latitudeOfFirstGridPointInDegrees']
            lon1 = self['longitudeOfFirstGridPointInDegrees']
            nx = self['Ni']
            ny = self['Nj']
            dx = self['DxInMetres']
            dy = self['DyInMetres']

            pj = pyproj.Proj(self.projparams)
            llcrnrx, llcrnry = pj(lon1, lat1)
            x = llcrnrx + dx * np.arange(nx)
            y = llcrnry + dy * np.arange(ny)
            x, y = np.meshgrid(x, y)
            lons, lats = pj(x, y, inverse=True)

        else:
            raise ValueError('unsupported grid {0}'.format(self['gridType']))

        return lats, lons

class fopen():

    def __init__(self, filename):
        self.f = open(filename, 'rb')

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __next__(self):

        gid = eccodes.codes_grib_new_from_file(self.f)

        if gid is None:
            raise StopIteration()

        self.data = cgrib(gid)
        eccodes.codes_release(gid)
        return self.data


    def __exit__(self, exc_type, exc_value, traceback):
        self.f.close()



