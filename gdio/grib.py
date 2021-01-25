__author__ = "Rodrigo Yamamoto"
__date__ = "2020.Dez"
__credits__ = ["Rodrigo Yamamoto", "Igor Santos"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.1.8.5"
__license__ = "MIT"
__status__ = "development"
__description__ = "A grib file IO library"

import logging
from datetime import datetime, timedelta

import numpy as np

from gdio import cgrib
from gdio.commons import near_yx, objectify


class grib(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.unitOfTimeRange = {
                                    None: None,
                                    0: 'minutes',
                                    1: 'hours',
                                    2: 'days',
                                    3: 'months',
                                    4: 'years',
                                    5: 'decades',
                                    6: 'normal',
                                    7: 'century',
                                    10: '3 hours',
                                    11: '6 hours',
                                    12: '12 hours',
                                    13: 'seconds',
                                    255: 'missing'
                                }
        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME', 'ref_time', 'time_units']
        self.__fields_3dlevel = ['isobaricInhPa', 'hybrid', 'sigma', 'eta',
                                 'heightAboveGround','depthBelowLandLayer','pressureFromGroundLayer']
        self.fields_ensemble = 'perturbationNumber'
        self.fields_ensemble_exception = [0]

        self.centre = 0
        self.lon = None
        self.lat = None
        self.time = None
        self.grid_description = None
        self.time_units = None
        self.history = None

        logging.basicConfig(datefmt='%Y%-m-%dT%H:%M:%S', level=logging.DEBUG,
                            format='[%(levelname)s @ %(asctime)s] %(message)s')



    def gb_load(self, ifile,
                vars=None,
                level_type=None,
                cut_time=None,
                cut_domain=None,
                filter_by={},
                rename_vars={}):
        '''
        Load grib file
        Yamamoto, R @ Out.2019
        :param ifile:       string
                            grib 1 or 2 file name
        :param vars:        list
                            variables short name or id parameter number
        :param cut_time:    tuple
                            range of time to cut ex.: (0,10)/(0,None)/(None,10)
        :param cut_domain:  tuple
                            range of latitudes and longitudes to cut: (lat1, lon1, lat2, lon2)
                            ex.: (-45,290,20,330)/(-45,None,20,330)/(None,290,None,320)
        :param level_type:  list
                            type of level (hybrid, isobaricInhPa, surface)
        :param filter_by:   dictonary
                            dict with grib parameters at form of pair key:values (list or single values)
                            eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]}
                            or filter_by={'gridType': 'regular_ll'}
        :param rename_vars: dictonary
                            rename variables names (key) for a new name (value).
                            Eg. {'tmpmdl': 't', 'tmpprs': 't'}
        :return:            dictonary/attributes
                            multiple time data container
        '''
        _data = objectify()
        data = objectify()

        try:

            _gb = cgrib.fopen(ifile)

            msg = [g for g in _gb]
            msg.sort(key=lambda x: (x.validityDate, x.validityTime, x.paramId, x.typeOfLevel, x.level, x.perturbationNumber))

            forecastDate = None
            fcst_time = 0
            concat_time = False
            ref_time = None

            for n, gr in enumerate(msg):

                start = 0
                stop = len(msg)
                cut_domain_roll = 0

                #filter by grib parameter
                if all([gr[k] in (v if isinstance(v, list) else [v]) for k,v in filter_by.items() if k in gr.keys()]):

                    # initialize time
                    if forecastDate is None:
                        self.history = "Created by gdio @ {date:%Y%m%d%H}".format(date=datetime.now())
                        ref_time = datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)
                        self.grid_description = {k: v for k,v in gr.items() if k in gr.gridkeys}
                        self.centre = gr.centre

                    # set time coordinate ....................
                    if not forecastDate == self.fcstTime(gr):
                        concat_time = True
                        forecastDate = self.fcstTime(gr)
                        fcst_time = int((forecastDate-ref_time).total_seconds()/(self.__unity(gr)*3600))
                        member_num = 0

                    # set temporal subdomain .......
                    if isinstance(cut_time, tuple):
                        start, stop = cut_time
                        start = 0 if start is None else start
                        stop = len(msg) if stop is None else stop


                    if self.verbose:
                        logging.debug('''forecastDate: {0} / cut_time_range: {1}-{2}'''.format(forecastDate,start,stop))


                    # cut time between start and stop time
                    if (not cut_time or (fcst_time >= start and fcst_time <= stop)):

                        typLev = gr.typeOfLevel

                        if (level_type is None or typLev in level_type):

                            if self.verbose:
                                logging.debug('''{0}, {1}, {2}, {3}, {4}, {5} '''.format(gr.dataDate,
                                                             gr.shortName,
                                                             gr.paramId,
                                                             gr.name,
                                                             gr.typeOfLevel,
                                                             gr.level,
                                                             gr.values.shape)
                                              )

                            # handle the variable id
                            idVar = gr.shortName if not gr.shortName in ['', 'unknown'] else str(gr.paramId)

                            # rename variables
                            for k, v in rename_vars.items():
                                if idVar in k:
                                    idVar = v

                            # concatenate variables .......................................
                            if vars is None or gr.shortName in vars or gr.paramId in vars:

                                # setup time ref/unity ...............
                                if not ('ref_time' in data.keys() or 'time_units' in data.keys()):
                                    data.update({'ref_time': ref_time})
                                    unit_time_range = gr.get('unitOfTimeRange', 255)
                                    data.update({'time_units': self.unitOfTimeRange[unit_time_range]})
                                    self.time_units = '{0} since {1}'.format(self.unitOfTimeRange[unit_time_range], ref_time)

                                # merge time ...............
                                if concat_time:
                                    data = self.__concat_time(_data, data, fcst_time)
                                    concat_time = False
                                    _data = objectify()


                                # set spatial coordinates ......
                                self.lat, self.lon = gr.latlons()

                                # convert from -180,180 to 360 format
                                self.lon = (self.lon + 360) % 360

                                flip_lat = self.lat[-1, 0] < self.lat[0, 0]

                                if flip_lat:                    #error with lat/lon 2 dims arrays
                                    self.lat = np.flip(self.lat, axis=0)


                                # select spatial subdomain .......
                                y, x = [None, None], [None, None]

                                if cut_domain:
                                    if isinstance(cut_domain, tuple):
                                        lat1, lon1, lat2, lon2 = cut_domain
                                        while True: # necessary 2 pass to fix 360 - 0 descontinuity
                                            y, x = near_yx({'latitude': self.lat[:, 0], 'longitude': self.lon[0, :]},
                                                                lats=[lat1, lat2], lons=[lon1, lon2])

                                            #if x0>x1 the longitude is rolled of x0 elements
                                            #in order to avoid discontinuity 360-0 of the longitude
                                            try:
                                                if x[0]>x[1]:
                                                    cut_domain_roll = -x[0]
                                                    self.lon = np.roll(self.lon, cut_domain_roll, axis=1)
                                                else:
                                                    break
                                            except:
                                                break


                                # trim lat/lon dimensions ......... Warning except lambert proj
                                self.lat = self.lat[y[0]:y[1], 0]
                                self.lon = self.lon[0, x[0]:x[1]]


                                # if necessary roll longitude due discontinuity 360-0 of the longitude
                                gr.values = np.roll(gr.values, cut_domain_roll, axis=-1)

                                # get data ........................
                                # grab data and flip the latitude axis if necessary
                                if flip_lat:
                                    _tmp = np.flip(gr.values, axis=0)[None, None, None, y[0]:y[1], x[0]:x[1]]
                                else:
                                    _tmp = gr.values[None, None, None, y[0]:y[1], x[0]:x[1]]


                                if idVar in _data.keys() and typLev in _data[idVar].keys():
                                    # concatenate levels
                                    if typLev in self.__fields_3dlevel:
                                        _tmp = np.concatenate((_data[idVar][typLev].value, _tmp), axis=2)
                                        _data[idVar][typLev].level.append(gr.level)

                                    # concatenate members
                                    if gr.get('perturbationNumber',0)>member_num:
                                        member_num = gr.get('perturbationNumber')
                                        _data[idVar][typLev].value = np.concatenate((_data[idVar][typLev].value, _tmp),
                                                                                    axis=0)
                                    else:
                                        _data[idVar][typLev].value = _tmp

                                else:
                                    member_num = 0

                                    __tmp = {
                                                typLev: {'value': _tmp, 'level': [gr.level]},
                                                'param_id': gr.paramId,
                                                'long_name': gr.name,
                                                'parameter_units': gr.parameterUnits,
                                                'latitude': self.lat,
                                                'longitude': self.lon
                                    }

                                    if idVar in _data.keys():
                                        _data[idVar].update(__tmp)
                                        _data[idVar].level_type.append(typLev)
                                    else:
                                        _data[idVar] = __tmp
                                        _data[idVar].level_type = [typLev]

                # consolidate data for last time block  ................
                if n + 1 == len(msg):
                    data = self.__concat_time(_data, data)

                self.variables = list(data.keys())
                self.coordinates.append('latitude')
                self.coordinates.append('longitude')
                self.coordinates.append('level')

        except Exception as e:
            logging.error('''gdio.gb_load > {0}'''.format(e))

        return data


    def gb_write(self, ifile,
                 data,
                 gri_format='GRIB1'):
        '''
        Write grib file

        :param ifile:           string
                                file path
        :param data:            dict
                                dataset
        :param gri_format:   string
                                netcdf format: GRIB1 or GRIB2
        :return:
        '''

        # print(data)
        print(data.keys())
        #ensemble loop
        # time loop
        for t, timestep in enumerate(data.time):

            print(t, timestep)
            # variable loop

            for idVar in self.__get_vars(data):
                print(idVar)

                # if self.verbose:
                #     logging.debug('''writing {var}'''.format(var=key))
                # levels lopp

        return


    def __get_vars(self, data):
        '''
        Extract variables from data keys
        :param data:    dict
                        data dictionary
        :return:        list
                        list of data variables
        '''
        b = self.__fields_latitude + self.__fields_longitude + self.__fields_time
        return list(set(data.keys())-set(b))


    def fcstTime(self, gr):
        '''
        Convert ref time + time step to forecast time
        :param gr:       object
                         grib message object
        '''

        if gr.step > 0:
            return datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute) + timedelta(hours=gr.step*self.__unity(gr))
        else:
            return datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)


    def __unity(self, gr):
        '''
        Scale fator to time unity transformation
        :param gr:  object
                    pygrig object
        :return:    float
        '''

        scale = 1.0

        if gr.stepUnits == 0:  # minute
            scale = 1/60
        elif gr.stepUnits == 1:  # hour
            scale = 1
        elif gr.stepUnits == 2:  # day
            scale = 24
        elif gr.stepUnits == 3:  # month
            scale = 24 * 30
        elif gr.stepUnits == 4:  # year
            scale = 24 * 365
        elif gr.stepUnits == 5:  # decade
            scale = 10 * 24 * 365

        return scale


    def __concat_time(self, _data, __data, fcst_time=None):
        '''
        Concatenate time dimension
        Yamamoto, Rodrigo @ Out.2019
        :param _data:       dict
                            single time data source
        :param __data:      dict
                            multiple time data container
        :param fcst_time:   int
                            timestep
        :return:            dict
                            multiple time data container
        '''

        for k,v in _data.items():
            if k in __data.keys():
                for l in v.keys():
                    if not l in ['level_type','param_id','long_name',
                                 'parameter_units','latitude','longitude']:
                        __data[k][l].value = np.concatenate((__data[k][l].value, _data[k][l].value), axis=1)
            else:
                __data[k] = _data[k]


        # add time information
        if not fcst_time is None:
            if 'time' in __data.keys():
                __data['time'] = np.unique(np.concatenate((__data['time'], np.array([fcst_time]))))
            else:
                __data['time'] = np.array([fcst_time])

            self.time = __data['time']

        return __data


    @staticmethod
    def is_grib(ifile):
        '''
        Check if is grib file
        from Igor@Out.2019
        :rtype: bool
        :return:
        '''

        if isinstance(ifile, str):
            with open(ifile, 'rb') as f:
                header = str(f.readline()[:20])
                f.close()
                if 'GRIB' in header:
                    return True

        return False


