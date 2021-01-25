__author__ = "Rodrigo Yamamoto"
__date__ = "2021.Jan"
__credits__ = ["Rodrigo Yamamoto"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.1.8.5"
__license__ = "MIT"
__status__ = "development"
__description__ = "A simple and concise gridded data IO library for read multiples grib and netcdf files"

import copy
import logging
import multiprocessing
import os
import warnings
from datetime import datetime, timedelta
from functools import partial

import numpy as np
import numpy.ma as ma

from gdio.commons import near_yx, objectify
from gdio.grib import grib as gblib
from gdio.netcdf import netcdf as nclib

warnings.filterwarnings("ignore")


class gdio(object):

    def __init__(self,
                 verbose=False,
                 remap_n_processes=2):

        self.verbose = verbose
        self.remap_n_processes = remap_n_processes
        self.dataset = list()

        self.coordinates = list()
        self.variables = list()

        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME']
        self.__fields_level = ['level', 'lev', 'LEVEL', 'levels', 'LEVELS']
        self.fields_ensemble = 'perturbationNumber'
        self.fields_ensemble_exception = [0]

        self.lon = None
        self.lat = None
        self.time = None
        self.time_units = None
        self.history = None

        logging.basicConfig(datefmt='%Y%-m-%dT%H:%M:%S', level=logging.DEBUG,
                            format='[%(levelname)s @ %(asctime)s] %(message)s')



    def thread(self, ifile, vars=None, cut_time=None, cut_domain=None, level_type=None, filter_by={}):
        '''
        Load and cutting function
        :param ifile:               string
                                    filename
        :param vars:                list
                                    lista de variaveis
        :param cut_time:            tuple
                                    range of time to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:          tuple
                                    range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                                    ex.: (-45,-90,20,-30)/(-45,None,20,-30)/(None,-90,None,-20)
        :level_type:                list
                                    type of level (hybrid, isobaricInhPa, surface)
        :param filter_by:           dictonary
                                    dict with grib parameters at form of pair key:values (list or single values)
                                    eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]}
                                    or filter_by={'gridType': 'regular_ll'}
        :return:                    dictionary
        '''

        if os.path.isfile(ifile):

            logging.info('''[PID:{0}] io.thread > opening file: {1}'''.format(os.getpid(), ifile))

            _data = None

            gb = gblib(verbose=self.verbose)
            gb.fields_ensemble = self.fields_ensemble
            gb.fields_ensemble_exception = self.fields_ensemble_exception

            nc = nclib(verbose=self.verbose)

            if gb.is_grib(ifile):
                return gb.gb_load(ifile, vars=vars,
                                     cut_time=cut_time,
                                     cut_domain=cut_domain,
                                     level_type=level_type,
                                     filter_by=filter_by)
            else:
                return nc.nc_load(ifile, vars=vars,
                                     cut_time=cut_time,
                                     cut_domain=cut_domain,
                                     level_type=level_type)

        else:
            logging.warning('''[PID:{0}] io.thread > missing file: {1}'''.format(os.getpid(),ifile))
            return None


    def mload(self,
              files,
              vars=None,
              merge_files=True,
              cut_time=None,
              cut_domain=None,
              level_type=None,
              filter_by={},
              uniformize_grid=True,
              inplace=False):
        '''
        Load multiple grib/netcdf files
        :param files:               list
                                    files names
        :param uniformize_grid:     boolean
                                    interpolate all ncs to first nc grid specification
        :param vars:                list
                                    variables names
        :param merge_files:         boolean
                                    merge files
        :param cut_time:            tuple
                                    range of time to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:          tuple
                                    range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                                    ex.: (-45,-90,20,-30)/(-45,None,20,-30)/(None,-90,None,-20)
        :param level_type:          list
                                    type of level (hybrid, isobaricInhPa, surface)
        :param filter_by:           dictonary
                                    dict with grib parameters at form of pair key:values (list or single values)
                                    eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]}
                                    or filter_by={'gridType': 'regular_ll'}
        :param rename_vars:         dictonary
                                    rename variables names (key) for a new name (value).
                                    Eg. {'tmpmdl': 't', 'tmpprs': 't'}
        :return:                    list of dictionaries
        '''

        data = objectify()
        griddes = ()

        ref_time = None
        t_units = 1

        # convert timestep index to timeserie ......
        def dtp(t, unity=1):
            return timedelta(days=float(t * unity))

        vf = np.vectorize(dtp)
        # ..........................................

        pool = multiprocessing.Pool(processes=self.remap_n_processes)

        if isinstance(files, str):
            files = [files]


        for _dat in pool.map(
                partial(self.thread, vars=vars,
                        cut_time=cut_time,
                        cut_domain=cut_domain,
                        level_type=level_type,
                        filter_by=filter_by),
                files):

            if _dat:

                ref_time = _dat.get('ref_time')

                # setting the standard projection/dimensions
                if not griddes:
                    lons_n, lats_n = self.__get_dims(_dat, vars)
                    griddes = lats_n.shape + lons_n.shape


                for key, val in _dat.items():

                    # convert to day unity
                    if _dat.get('time_units').lower() in ['hour', 'hours', 'hrs']:
                        t_units = 1 / 24
                    else:
                        t_units = 1

                    if (vars is None or key in vars) \
                            and not key in ['latitude', 'longitude', 'ref_time', 'time', 'time_units']:

                        for typLev in val.level_type:

                            # uniformize all grids ...........
                            if uniformize_grid:

                                # grid resample, if spatial dimensions are different of first grid(z,lat,lon)
                                if not val[typLev].value.shape[3:] == griddes:

                                    logging.info('''gdio.mload > auto remapping grid @ {0}'''.format(key))

                                    # interpolate through z dimension .........
                                    _tmp = np.ones(val[typLev].value.shape[:3]+griddes) * np.nan

                                    # WARNING: If the z dimension of the source data is different
                                    # from that of the ref data, the interpolated level data may not
                                    # represent the same z as the ref data
                                    for m in range(_tmp.shape[0]):
                                        for z in range(_tmp.shape[2]):
                                            try:
                                                _tmp[m,:,z,:,:] = self.remapbil(val[typLev].value[m,:,z,:,:],
                                                                                 val.longitude, val.latitude,
                                                                                 lons_n, lats_n, order=1, masked=True)
                                            except Exception as e:
                                                logging.error('''gdio.mload > auto remapping grid error {0}'''.format(e))

                                    val[typLev].value = _tmp

                                    del val.longitude, val.latitude
                                    del _tmp

                            # update the lat/lon dimensions
                            data['longitude'], data['latitude'] = lons_n, lats_n

                            # merge files ........................
                            if merge_files:

                                if key in data.keys() and typLev in data[key].keys():

                                    if not (key in self.__fields_latitude
                                            or key in self.__fields_longitude
                                            or key in self.__fields_time
                                            or key in self.__fields_level
                                            or key in ['ref_time', 'time_units']):                            # merge variable field

                                        try:
                                            data[key][typLev].value = np.concatenate((data[key][typLev].value,
                                                                                      val[typLev].value), axis=1)
                                        except Exception as e:
                                            logging.error('''gdio.mload > error @ {0} - {1}'''.format(key, e))
                                else:
                                    if key in data.keys():  # in case of multiples level per variable
                                        data[key].update(val)
                                    else:
                                        data[key] = val


                    else:  # all parameters except variables

                        # set datetime field
                        if key in self.__fields_time:
                            if key in data.keys():  # merge datetime field
                                try:
                                    _time = ref_time + vf(_dat.get('time'), t_units)
                                    data[key] = np.concatenate((data[key], _time))
                                except Exception as e:
                                    logging.error('''gdio.mload > error @ {0} - {1}'''.format(key, e))
                            else:
                                data['time'] = ref_time + vf(_dat.get('time'), t_units)
                                data['ref_time'] = _dat.get('ref_time')
                        else:
                            if key not in data.keys():
                                data[key] = val

                # do not merge files option ..............
                if not merge_files:
                    data.update({'time': ref_time + vf(_dat.get('time'), t_units)})
                    data.update({'ref_time': [_dat.get('ref_time')]})

                    self.dataset.append(data)
                    data = objectify()

            else:
                # in case of missing file ................
                for key in data.keys():

                    if not (key in self.__fields_latitude
                            or key in self.__fields_longitude
                            or key in self.__fields_time
                            or key in self.__fields_level
                            or key in ['ref_time', 'time_units']):
                        for typLev in val.level_type:
                            data[key][typLev].value = np.concatenate((data[key][typLev].value,
                                                            np.ones((1,1)+data[key][typLev].value.shape[2:]) * np.nan),
                                                            axis=1)
                    elif key in self.__fields_time:
                        data[key] = np.concatenate((data[key], [data[key][-1] + timedelta(days=t_units)]))
                    elif key in ['ref_time']:
                        ref_time += timedelta(days=t_units)
                        data[key] = np.concatenate((data[key], [ref_time]))

                logging.warning('''io.load_nc > missing file applying null grid''')

        self.variables = list(data.keys())
        self.coordinates.append('latitude')
        self.coordinates.append('longitude')
        self.coordinates.append('level')

        if inplace:
            if data:
                self.dataset.append(data)
        else:
            return data


    def sel(self,
            __data=None,
            latitude=None,
            longitude=None,
            dates=None,
            level=None,
            member=None,
            date_format='%Y-%m-%d %H:%M'):
        '''
        Select data by coordinates (date, latitude, longitude, levels and members)

        :param __data:       list of dictionaries
                             raw dataset
        :param latitude:     list of floats
                             latitudes
                             range of latitudes to select: [lat1, lat2]
                             especific latitudes (1 or >2) [lat1, lat2, lat3, ...]
        :param longitude:    list of floats
                             range of longitudes to select: [lon1, lon2]
                             especific longitudes (1 or >2) [lon1, lon2, lon3, ...]
        :param dates:        list of datetime/string
                             datetime/string date
                             range of dates to select: [date1, date2]
                             especific dates (1 or >2) [date1, date2, date3, ...]
        :param level:        list of int
                             range of levels to select: [level1, level2]
                             especific levels (1 or >2) [level1, level2, level3, ...]
        :param member:       list of int
                             range of levels to select: [member, member]
                             especific levels (1 or >2) [level1, level2, level3, ...]
        return               dict
        '''

        if __data is None:
            __data = copy.deepcopy(self.dataset)

        for _dat in __data:

            t = dates
            x = None
            y = None
            z = level

            # select time
            if dates:

                for i, dt in enumerate(dates):
                    if isinstance(dt, str):
                        dates[i] = datetime.strptime(dt, date_format)

                if len(dates) == 2:
                    t = (_dat.get('time') >= dates[0]) & (_dat.get('time') <= dates[1])
                elif len(dates) > 0:
                    t = np.isin(_dat.get('time'), dates)


            # select spatial subdomain
            if longitude or latitude:
                y, x = near_yx(_dat, lats=latitude, lons=longitude)


            for k, v in _dat.items():

                # cutting data array
                if isinstance(v, dict):

                    for typLev in v.level_type:

                        # cut data in longitude dimension
                        if x:
                            if len(x) == 2:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, :, x[0]:x[1]]
                            elif len(x) == 1:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, :, x[0]]
                            else:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, :, x]

                        # cut data in latitude dimension
                        if y:
                            if len(y) == 2:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, y[0]:y[1]]
                            elif len(y) == 1:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, y[0]]
                            else:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, :, y]

                        # cut data in levels dimension
                        if z:
                            if len(z) == 2:
                                _dat[k][typLev].value = _dat[k][typLev].value[:, :, z[0]:z[1]]
                                _dat[k][typLev].level = _dat[k][typLev].level[z[0]:z[1]]
                            else:
                                try:
                                    _dat[k][typLev].value = _dat[k][typLev].value[:, :, z]
                                    _dat[k][typLev].level = list(map(_dat[k][typLev].level.__getitem__, z))
                                except:
                                    _dat[k][typLev].value = _dat[k][typLev].value[:, :, -1]
                                    _dat[k][typLev].level = _dat[k][typLev].level[-1]


                        # cut data in time dimension
                        if dates:
                            _dat[k][typLev].value = _dat[k][typLev].value[:, t]

                        # cut data in member dimension
                        if member:
                            _dat[k][typLev].value = _dat[k][typLev].value[member]



                # select cordinates attributes
                else:

                    if k in ['latitude']:  # latitude coordinate

                        if y:
                            if len(y) == 2:
                                _dat[k] = _dat[k][y[0]:y[1]]
                            else:
                                _dat[k] = _dat[k][y]

                    elif k in ['longitude']:  # longitude coordinate

                        if x:
                            if len(x) == 2:
                                _dat[k] = _dat[k][x[0]:x[1]]
                            else:
                                _dat[k] = _dat[k][x]

                    elif k in ['time']:  # time coordinate

                        if dates:
                            _dat[k] = _dat[k][t]

                    else:
                        _dat.update({k: v})

        return __data


    def remapbil(self, data, lon, lat, lon_new, lat_new, order=1, masked=False):
        '''
        Interpolate data to new domain resolution
        :param data:    array
                        3D data (time,lon,lat)
        :param lon:     array
        :param lat:     array
        :param lon_new: array
                        new grid logitudes
        :param lat_new: array
                        new grid latitudes
        :param order:   int
                        0- nearest-neighbor, 1 - bilinear, 2 - cubic spline
        :param masked:  boolean
                        If True, points outside the range of xin and yin
                        are masked (in a masked array). If masked is set to a number

        :return:        3D array
        '''

        _lon_new, _lat_new = np.meshgrid(lon_new, lat_new)

        cpu_num = multiprocessing.cpu_count()

        n_processes = cpu_num if self.remap_n_processes > cpu_num else self.remap_n_processes

        pool = multiprocessing.Pool(processes=n_processes)

        # here we parallelise in each step of time, a kind of magic
        return np.array(pool.map(
                                partial(self.interp, xin=lon[np.argsort(lon)],
                                        yin=lat[np.argsort(lat)], xout=_lon_new, yout=_lat_new,
                                        order=order, masked=masked),
                                data)
                        )


    def interp(self, datain, xin, yin, xout, yout, checkbounds=False, masked=False, order=1):
        """
        From basemap lib
        Interpolate data (``datain``) on a rectilinear grid (with x = ``xin``
        y = ``yin``) to a grid with x = ``xout``, y= ``yout``.
        .. tabularcolumns:: |l|L|
        ==============   ====================================================
        Arguments        Description
        ==============   ====================================================
        datain           a rank-2 array with 1st dimension corresponding to
                        y, 2nd dimension x.
        xin, yin         rank-1 arrays containing x and y of
                        datain grid in increasing order.
        xout, yout       rank-2 arrays containing x and y of desired output grid.
        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        checkbounds      If True, values of xout and yout are checked to see
                        that they lie within the range specified by xin
                        and xin.
                        If False, and xout,yout are outside xin,yin,
                        interpolated values will be clipped to values on
                        boundary of input grid (xin,yin)
                        Default is False.
        masked           If True, points outside the range of xin and yin
                        are masked (in a masked array).
                        If masked is set to a number, then
                        points outside the range of xin and yin will be
                        set to that number. Default False.
        order            0 for nearest-neighbor interpolation, 1 for
                        bilinear interpolation, 3 for cublic spline
                        (default 1). order=3 requires scipy.ndimage.
        ==============   ====================================================
        .. note::
        If datain is a masked array and order=1 (bilinear interpolation) is
        used, elements of dataout will be masked if any of the four surrounding
        points in datain are masked.  To avoid this, do the interpolation in two
        passes, first with order=1 (producing dataout1), then with order=0
        (producing dataout2).  Then replace all the masked values in dataout1
        with the corresponding elements in dataout2 (using numpy.where).
        This effectively uses nearest neighbor interpolation if any of the
        four surrounding points in datain are masked, and bilinear interpolation
        otherwise.
        Returns ``dataout``, the interpolated data on the grid ``xout, yout``.
        """
        # xin and yin must be monotonically increasing.
        if xin[-1] - xin[0] < 0 or yin[-1] - yin[0] < 0:
            raise ValueError('xin and yin must be increasing!')
        if xout.shape != yout.shape:
            raise ValueError('xout and yout must have same shape!')

        # check that xout,yout are
        # within region defined by xin,yin.
        if checkbounds:
            if xout.min() < xin.min() or \
                    xout.max() > xin.max() or \
                    yout.min() < yin.min() or \
                    yout.max() > yin.max():
                raise ValueError('yout or xout outside range of yin or xin')

        # compute grid coordinates of output grid.
        delx = xin[1:] - xin[0:-1]
        dely = yin[1:] - yin[0:-1]

        if max(delx) - min(delx) < 1.e-4 and max(dely) - min(dely) < 1.e-4:
            # regular input grid.
            xcoords = (len(xin) - 1) * (xout - xin[0]) / (xin[-1] - xin[0])
            ycoords = (len(yin) - 1) * (yout - yin[0]) / (yin[-1] - yin[0])
        else:
            # irregular (but still rectilinear) input grid.
            xoutflat = xout.flatten()
            youtflat = yout.flatten()
            ix = (np.searchsorted(xin, xoutflat) - 1).tolist()
            iy = (np.searchsorted(yin, youtflat) - 1).tolist()
            xoutflat = xoutflat.tolist()
            xin = xin.tolist()
            youtflat = youtflat.tolist()
            yin = yin.tolist()
            xcoords = []
            ycoords = []

            for n, i in enumerate(ix):
                if i < 0:
                    xcoords.append(-1)  # outside of range on xin (lower end)
                elif i >= len(xin) - 1:
                    xcoords.append(len(xin))  # outside range on upper end.
                else:
                    xcoords.append(
                        float(i) + (xoutflat[n] - xin[i]) / (xin[i + 1] - xin[i]))

            for m, j in enumerate(iy):
                if j < 0:
                    # outside of range of yin (on lower end)
                    ycoords.append(-1)
                elif j >= len(yin) - 1:
                    ycoords.append(len(yin))  # outside range on upper end
                else:
                    ycoords.append(
                        float(j) + (youtflat[m] - yin[j]) / (yin[j + 1] - yin[j]))

            xcoords = np.reshape(xcoords, xout.shape)
            ycoords = np.reshape(ycoords, yout.shape)

        # data outside range xin,yin will be clipped to
        # values on boundary.
        if masked:
            xmask = np.logical_or(np.less(xcoords, 0),
                                  np.greater(xcoords, len(xin) - 1))
            ymask = np.logical_or(np.less(ycoords, 0),
                                  np.greater(ycoords, len(yin) - 1))
            xymask = np.logical_or(xmask, ymask)

        xcoords = np.clip(xcoords, 0, len(xin) - 1)
        ycoords = np.clip(ycoords, 0, len(yin) - 1)

        # interpolate to output grid using bilinear interpolation.
        if order == 1:
            xi = xcoords.astype(np.int32)
            yi = ycoords.astype(np.int32)
            xip1 = xi + 1
            yip1 = yi + 1
            xip1 = np.clip(xip1, 0, len(xin) - 1)
            yip1 = np.clip(yip1, 0, len(yin) - 1)
            delx = xcoords - xi.astype(np.float32)
            dely = ycoords - yi.astype(np.float32)
            dataout = (1. - delx) * (1. - dely) * datain[yi, xi] + \
                      delx * dely * datain[yip1, xip1] + \
                      (1. - delx) * dely * datain[yip1, xi] + \
                      delx * (1. - dely) * datain[yi, xip1]
        elif order == 0:
            xcoordsi = np.around(xcoords).astype(np.int32)
            ycoordsi = np.around(ycoords).astype(np.int32)
            dataout = datain[ycoordsi, xcoordsi]
        elif order == 3:
            try:
                from scipy.ndimage import map_coordinates
            except ImportError:
                raise ValueError('scipy.ndimage must be installed if order=3')
            coords = [ycoords, xcoords]
            dataout = map_coordinates(datain, coords, order=3, mode='nearest')
        else:
            raise ValueError('order keyword must be 0, 1 or 3')

        if masked:
            newmask = ma.mask_or(ma.getmask(dataout), xymask)
            dataout = ma.masked_array(dataout, mask=newmask)
            if not isinstance(masked, bool):
                dataout = dataout.filled(masked)

        return dataout


    def __get_dims(self, data, var=None):
        '''
        Get grid data dimension
        :param data:    dictionary
                        data
        :param var:     string
                        grid reference variable name
        :return:        tuple
                        grid dimension, lat/lon of grid reference
        '''

        for key, val in data.items():
            if (var is None or key == var[0]) \
                    and not key in ['latitude', 'longitude', 'ref_time', 'time', 'time_units']:
                return val.longitude, val.latitude


