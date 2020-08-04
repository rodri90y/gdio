__author__ = "Rodrigo Yamamoto"
__date__ = "2020.Ago"
__credits__ = ["Rodrigo Yamamoto","Carlos Oliveira","Igor"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.0.8"
__license__ = "MIT"
__status__ = "development"
__description__ = "A simple and concise gridded data IO library for read multiples grib and netcdf files"

import os, copy

import numpy as np
import numpy.ma as ma
import multiprocessing
import logging

from functools import partial
from datetime import datetime, timedelta

from gdio.grib import grib as gblib
from gdio.netcdf import netcdf as nclib
from gdio.commons import near_yx

import warnings

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
        self.__fields_3dlevel = ['depthBelowLandLayer', 'isobaricInhPa', 'hybrid', 'sigma', 'eta']
        self.fields_ensemble = 'perturbationNumber'
        self.fields_ensemble_exception = [0]
        self.level = None
        self.lon = None
        self.lat = None
        self.time = None
        self.time_units = None
        self.history = None

        logging.basicConfig(handlers=[logging.StreamHandler()],
                            datefmt='%Y%-m-%dT%H:%M:%S', level=logging.DEBUG,
                            format='[%(levelname)s @ %(asctime)s] %(message)s')



    def thread(self, ifile, vars=None, cut_time=None, cut_domain=None, level_type=None):
        '''
        Load and cutting function
        :param conf:                list
                                    filename e index
        :param vars:                list
                                    lista de variaveis
        :param cut_time:            tuple
                                    range of time to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:          tuple
                                    range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                                    ex.: (-45,-90,20,-30)/(-45,None,20,-30)/(None,-90,None,-20)
        :level_type:                list
                                    type of level (hybrid, isobaricInhPa, surface)
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
                _data = gb.gb_load(ifile, vars=vars,
                                     cut_time=cut_time,
                                     cut_domain=cut_domain,
                                     level_type=level_type)
            else:

                _data = nc.nc_load(ifile, vars=vars,
                                     cut_time=cut_time,
                                     cut_domain=cut_domain)

            if any(x in ['lon', 'lat'] for x in _data.keys()):
                _data['longitude'] = _data.pop('lon')
                _data['latitude'] = _data.pop('lat')

            return _data

        else:
            logging.warning('''[PID:{0}] io.thread > missing file: {1}'''.format(os.getpid(),ifile))
            return None


    def mload(self,
              files,
              merge_files=True,
              cut_time=None,
              cut_domain=None,
              level_type=None,
              uniformize_grid=True,
              vars=None,
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
        :level_type:                list
                                    type of level (hybrid, isobaricInhPa, surface)
        :return:                    list of dictionaries
        '''

        data = {}
        griddes = ()

        ref_time = None
        t_unit = 1

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
                        level_type=level_type),
                files):


            if not _dat is None:

                ref_time = _dat['ref_time']

                for k in _dat.keys():

                    if (vars is None or k in vars) \
                            and not k in ['latitude', 'longitude', 'lat', 'lon',
                                          'level', 'ref_time', 'time', 'time_unit']:

                        # redim the data array ...........
                        if _dat[k].ndim == 2:
                            _dat[k] = _dat[k][None, None, :, :]
                        elif _dat[k].ndim == 3:
                            _dat[k] = _dat[k][:, None, :, :]


                        if not griddes:
                            griddes = _dat[k].shape
                            lons_n, lats_n = _dat['longitude'], _dat['latitude']


                        # uniformize all grids ...........
                        if uniformize_grid:

                                # grid resample
                            if _dat[k].ndim > 2 and \
                                    not _dat[k].shape[1:] == griddes[1:]:

                                logging.info('''gdio.mload > auto remapping grid @ {0}'''.format(k))

                                # interpolate through z dimension .........
                                _tmp = np.ones(_dat[k].shape[:1]+griddes[1:]) * np.nan

                                # WARNING: If the z dimension of the source data is different
                                # from that of the ref data, the interpolated level data may not
                                # represent the same z as the ref data
                                for z in range(_tmp.shape[1]):
                                    try:

                                        _tmp[:,z,:,:] = self.remapbil(_dat[k][:,z,:,:],
                                                                         _dat['longitude'], _dat['latitude'],
                                                                         lons_n, lats_n, order=1, masked=True)

                                    except Exception as e:
                                        logging.error('''gdio.mload > auto remapping grid error {0}'''.format(k), e)

                                _dat[k] = _tmp

                                del _tmp

                        # update the lat/lon dimensions
                        data['longitude'], data['latitude'] = lons_n, lats_n

                    # convert to day unity
                    if _dat['time_unit'] in ['hours', 'hrs']:
                        t_unit = 1 / 24
                    else:
                        t_unit = 1

                    # merge files ........................
                    if merge_files:

                        if k in data.keys():

                            if not (k in self.__fields_latitude \
                                    or k in self.__fields_longitude \
                                    or k in self.__fields_time \
                                    or k in self.__fields_level \
                                    or k in ['ref_time', 'time_unit']):                            # merge variable field
                                try:
                                    data[k] = np.concatenate((data[k], _dat[k]))
                                except Exception as e:
                                    logging.error('''gdio.mload > error @ {0} - {1}'''.format(k, e))

                            elif k in self.__fields_time:                                          # merge datetime field
                                _time = ref_time + vf(_dat['time'], t_unit)
                                try:
                                    data[k] = np.concatenate((data[k], _time))
                                except Exception as e:
                                    logging.error('''gdio.mload > error @ {0} - {1}'''.format(k, e))
                            elif k in self.__fields_level:                                         # merge level field
                                try:
                                    data[k] = np.unique(np.concatenate((data[k], _dat[k])))
                                except Exception as e:
                                    logging.error('''gdio.mload > error @ {0} - {1}'''.format(k, e))
                        else:                                                                      # set new field
                            data.update({k: _dat[k]})

                            # set datetime field
                            if k in self.__fields_time:
                                data['time'] = ref_time + vf(_dat['time'], t_unit)
                                data['ref_time'] = [_dat['ref_time']]

                    else:
                        data.update({k: _dat[k]})

                        # set datetime field
                        if k in self.__fields_time:
                            data['time'] = ref_time + vf(_dat['time'], t_unit)
                            data['ref_time'] = [_dat['ref_time']]



                # do not merge files option ..............
                if not merge_files:

                    data.update({'time': ref_time + vf(_dat['time'], t_unit)})
                    data.update({'ref_time': [_dat['ref_time']]})

                    self.dataset.append(data)
                    data = dict()

            else:

                # in case of missing file ................
                for k in data.keys():

                    if not (k in self.__fields_latitude \
                            or k in self.__fields_longitude \
                            or k in self.__fields_time \
                            or k in self.__fields_level \
                            or k in ['ref_time', 'time_unit']):
                        data[k] = np.concatenate((data[k], [data[k][-1] * np.nan]))
                    elif k in self.__fields_time:
                        data[k] = np.concatenate((data[k], [data[k][-1] + timedelta(days=t_unit)]))
                    elif k in ['ref_time']:
                        ref_time += timedelta(days=t_unit)
                        data[k] = np.concatenate((data[k], [ref_time]))


                logging.warning('''io.load_nc > missing file applying null grid''')


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
            date_format='%Y-%m-%d %H:%M'):
        '''
        Select data by coordinates (date, latitude, longitude and levels)

        :param __data:       list of dictionaries
                             raw dataset
        :param latitude:     list of floats
                             latitudes
                             range of latitudes to select: [lat1, lat2]
                             especific latitudes (1 or >2) [lat1, lat2, lat3, ...]
        :param longitude:    list of floats
                             range of longitudes to select: [lon1, lon2]
                             especific longitudes (1 or >2) [lon1, lon2, lon3, ...]
        :param dates:        list of floats
                             datetime/string date
                             range of dates to select: [date1, date2]
                             especific dates (1 or >2) [date1, date2, date3, ...]
        :param level:        list of floats
                             range of levels to select: [level1, level2]
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
                    t = np.where((_dat['time'] >= dates[0]) & (_dat['time'] <= dates[1]))
                elif len(dates) > 0:
                    t = np.isin(_dat['time'], dates)


            # select spatial subdomain
            if longitude or latitude:
                y, x = near_yx(_dat, lats=latitude, lons=longitude)


            for k, v in _dat.items():

                if isinstance(v, np.ndarray):

                    if k in ['latitude']:                       # latitude coordinate

                        if y:
                            if len(y) == 2:
                                _dat[k] = _dat[k][y[0]:y[1]]
                            else:
                                _dat[k] = _dat[k][y]

                    elif k in ['longitude']:                    # longitude coordinate

                        if x:
                            if len(x) == 2:
                                _dat[k] = _dat[k][x[0]:x[1]]
                            else:
                                _dat[k] = _dat[k][x]

                    elif k in ['time']:                         # time coordinate

                        if dates:
                            _dat[k] = _dat[k][t]

                    elif k in ['level']:                        # level coordinate
                        if z:
                            if len(z) == 2:
                                _dat[k] = _dat[k][z[0]:z[1]]
                            else:
                                _dat[k] = _dat[k][z]

                    elif not k in ['ref_time']:                 # variables (non coordinates)

                        # cut longitude
                        if x:
                            if len(x) == 2:
                                _dat[k] = _dat[k][:, :, :, x[0]:x[1]]
                            elif len(x) == 1:
                                _dat[k] = _dat[k][:, :, :, x[0]]
                            else:
                                _dat[k] = _dat[k][:, :, :, x]

                        # cut latitude
                        if y:
                            if len(y) == 2:
                                _dat[k] = _dat[k][:, :, y[0]:y[1]]
                            elif len(y) == 1:
                                _dat[k] = _dat[k][:, :, y[0]]
                            else:
                                _dat[k] = _dat[k][:, :, y]

                        # cut levels
                        if z:
                            if len(z) == 2:
                                _dat[k] = _dat[k][:, z[0]:z[1]]
                            else:
                                try:
                                    _dat[k] = _dat[k][:, z]
                                except:
                                    _dat[k] = _dat[k][:, -1]
                                    # cut time
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





