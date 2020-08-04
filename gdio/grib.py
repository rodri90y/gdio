__author__ = "Rodrigo Yamamoto"
__date__ = "2020.Ago"
__credits__ = ["Rodrigo Yamamoto","Carlos Oliveira","Igor"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.0.8"
__license__ = "MIT"
__status__ = "development"
__description__ = "A grib file IO library"

import numpy as np
import pygrib
import logging

from gdio.commons import near_yx
from texttable import Texttable
from datetime import datetime, timedelta


class grib(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.__fields_3dlevel = ['isobaricInhPa', 'hybrid', 'sigma', 'eta',
                                 'heightAboveGround','depthBelowLandLayer','pressureFromGroundLayer']
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



    def gb_info(self, ifile):
        '''
        List grib info
        Yamamoto, R @ Out.2020
        :param ifile:       string
                            grib 1 or 2 file name

        :return:            dict
        '''
        try:
            _gb = pygrib.open(ifile)

            msg = [g for g in _gb]
            msg.sort(key=lambda x: (x.validDate, x.typeOfLevel, x.paramId, x.level))

            table = Texttable()
            table.set_deco(Texttable.HEADER)
            table.set_cols_dtype(['t', 'i', 't', 'i', 'f', 'f', 'a'])
            table.set_cols_align(["l", "c", "l", "l", "l", "l", "l"])
            table.set_cols_width([8, 4, 10, 7, 35, 20, 7])
            table.add_rows([["Time", "N", "shortName", "paramId", "longName", "typeOfLevel", "level"]])

            for k, gr in enumerate(msg):
                table.add_row([gr.dataDate, gr.messagenumber, gr.shortName, gr.paramId,
                               gr.name, gr.typeOfLevel, gr.level])

            return table.draw()

        except Exception as e:
            logging.error('''gdio.gb_load > '''.format(e))


    def gb_load(self, ifile,
                vars=None,
                level_type=None,
                cut_time=None,
                cut_domain=None):
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
        :level_type:        list
                            type of level (hybrid, isobaricInhPa, surface)
        :return:            dict
                            multiple time data container
        '''
        _data = dict()
        data = dict()


        # try:
        _gb = pygrib.open(ifile)

        msg = [g for g in _gb]
        msg.sort(key=lambda x: (x.validDate, x.step, x.paramId, x.typeOfLevel, x.level))

        forecastDate = None
        fcst_time = 0
        concat_time = False

        for k, gr in enumerate(msg):
            start = 0
            stop = len(msg)

            # initialize time
            if forecastDate is None:
                self.history = None
                ref_time = datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)

                data.update({'ref_time': ref_time})
                data.update({'time_unit': gr.fcstimeunits})
                self.time_units = '{0} since {1}'.format(gr.fcstimeunits, ref_time)


            # set time coordinate ....................
            if not forecastDate == self.fcstTime(gr):
                concat_time = True
                forecastDate = self.fcstTime(gr)
                fcst_time = int((forecastDate-ref_time).total_seconds()/(self.__unity(gr)*3600))


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

                    if not typLev in ['surface', 'isobaricInhPa']:
                        idVar = f'{idVar}_{typLev}'.replace(' ', '_')


                    if self.fields_ensemble in gr.keys():
                        if not gr[self.fields_ensemble] in self.fields_ensemble_exception:
                            idVar += '_m{0}'.format(gr[self.fields_ensemble]).replace(' ', '_')


                    # concatenate variables .......................................
                    if vars is None or gr.shortName in vars or gr.paramId in vars:

                        # merge time ...............
                        if concat_time:
                            data = self.__concat_time(_data, data, fcst_time)
                            concat_time = False
                            _data = dict()


                        # set spatial coordinates ......
                        self.lat, self.lon = gr.latlons()

                        # convert from -180,180 to 360 format
                        self.lon = (self.lon + 360) % 360

                        flip_lat = self.lat[-1, 0] < self.lat[0, 0]

                        if flip_lat:
                            self.lat = np.flip(self.lat, axis=0)


                        # select spatial subdomain .......

                        y, x = [None, None], [None, None]

                        if cut_domain:
                            if isinstance(cut_domain, tuple):
                                lat1, lon1, lat2, lon2 = cut_domain

                                y, x = near_yx({'latitude': self.lat[:, 0], 'longitude': self.lon[0, :]},
                                                    lats=[lat1, lat2], lons=[lon1, lon2])


                        # trim lat/lon dimensions .........
                        self.lat = self.lat[y[0]:y[1], 0]
                        self.lon = self.lon[0, x[0]:x[1]]

                        data.update({'latitude': self.lat})
                        data.update({'longitude': self.lon})


                        # get data ........................
                        if idVar in _data.keys():
                            try:
                                # concatenate levels
                                if typLev in self.__fields_3dlevel:
                                    if flip_lat:
                                        _data[idVar] = np.concatenate((_data[idVar],
                                                                       np.flip(gr.values, axis=0)[None, None,
                                                                       y[0]:y[1], x[0]:x[1]]),
                                                                      axis=1)
                                    else:
                                        _data[idVar] = np.concatenate((_data[idVar],
                                                                       gr.values[None, None, y[0]:y[1],
                                                                       x[0]:x[1]]),
                                                                      axis=1)
                            except Exception as e:
                                print('[E] gdio.gb_load > ', e)
                        else:
                            if flip_lat:
                                _data.update({idVar: np.flip(gr.values, axis=0)[None, None, y[0]:y[1], x[0]:x[1]]})
                            else:
                                _data.update({idVar: gr.values[None, None, y[0]:y[1], x[0]:x[1]]})


            # consolidate data for last time block  ................
            if k + 1 == len(msg):
                data = self.__concat_time(_data, data)

        self.variables = list(data.keys())
        self.coordinates.append('latitude')
        self.coordinates.append('longitude')

        # except Exception as e:
        #     logging.error('''gdio.gb_load > '''.format(e))

        return data


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
        for k in _data.keys():
            if k in __data.keys():
                __data[k] = np.concatenate((__data[k], _data[k]), axis=0)
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
        :param file_or_buffer:
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

