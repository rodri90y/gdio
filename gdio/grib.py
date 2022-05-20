__author__ = "Rodrigo Yamamoto"
__date__ = "2022.Mai"
__credits__ = ["Rodrigo Yamamoto", "Igor Santos"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.2.9"
__license__ = "MIT"
__status__ = "development"
__description__ = "A grib file IO library"

import logging
from datetime import datetime, timedelta
import os
import numpy as np

from gdio import cgrib
from gdio.commons import near_yx2, objectify, dict_get, timestep_to_datetime, datetime_to_timestep
from .definitions.Table_4_4 import UNIT_TIME_RANGE

class grib(object):

    def __init__(self, verbose=False):

        self.verbose = verbose

        self.coordinates = list()
        self.variables = list()

        self.unitOfTimeRange = UNIT_TIME_RANGE

        self.__fields_latitude = ['latitude', 'lat', 'xlat', 'LATITUDE']
        self.__fields_longitude = ['longitude', 'lon', 'xlon', 'LONGITUDE']
        self.__fields_time = ['time', 'TIME', 'ref_time', 'time_units']
        self.__fields_3dlevel = ['isobaricInhPa', 'hybrid', 'sigma', 'eta', 'theta',
                                 'sigmaLevel', 'isentropic']
        self.fields_ensemble = 'perturbationNumber'
        self.__non_data_variables = [
                                         'centre',
                                         'dataType',
                                         'level_type',
                                         'param_id',
                                         'long_name',
                                         'parameter_units',
                                         'latitude',
                                         'longitude',
                                         'grid_type',
                                         'projparams'
                                     ]
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
                rename_vars={},
                sort_before=False):
        '''
        Load grib file
        Yamamoto, R @ Mar.2022
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
        :param sort_before: bool
                            Sort fields before process validityDate, validityTime, paramId, typeOfLevel, perturbationNumber and level
                            Warning high consumption of memory, just use when the grib data structure is not standard
        :return:            dictonary/attributes
                            multiple time data container
        '''

        _data = objectify()
        data = objectify()

        try:

            with cgrib.fopen(ifile) as msg:

                    # sort fields before use, warning high consumption of memory
                    if sort_before:
                        msg = [g for g in msg
                               if (vars is None or g.shortName in vars or g.paramId in vars)
                               and all(
                                   [g[k] in (v if isinstance(v, list) else [v]) for k, v in filter_by.items() if k in g.keys()])
                               ]

                        msg.sort(key=lambda x: (
                            x.dataDate, x.dataTime, x.step, x.paramId, x.typeOfLevel, x.level, x.perturbationNumber))

                    forecastDate = None
                    fcst_time = 0
                    step_time = -1
                    concat_time = False
                    ref_time = None
                    msg_len = len(msg)

                    for n, gr in enumerate(msg):

                        start = 0
                        stop = float('inf')
                        cut_domain_roll = 0

                        # filter by grib parameter
                        if all([gr[k] in (v if isinstance(v, list) else [v]) for k, v in filter_by.items() if k in gr.keys()]):

                            # initialize time
                            if forecastDate is None:
                                self.history = "Created by gdio @ {date:%Y%m%d%H}".format(date=datetime.now())
                                ref_time = datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)
                                self.grid_description = {k: v for k, v in gr.items() if k in gr.gridkeys}
                                self.centre = gr.centre

                            # set time coordinate ....................
                            if not forecastDate == self.fcstTime(gr):
                                concat_time = True
                                forecastDate = self.fcstTime(gr)
                                fcst_time = int((forecastDate - ref_time).total_seconds() / (self.__unity(gr) * 3600))
                                step_time += 1
                                member_num = 0

                            # set temporal subdomain .......
                            if isinstance(cut_time, (tuple, list)):
                                start, stop = cut_time
                                start = 0 if start is None else start
                                stop = float('inf') if stop is None else stop

                            if self.verbose:
                                logging.debug(f'forecastDate: {forecastDate} / cut_time_range: {start}-{stop} / {gr.shortName}')

                            # cut time between start and stop time
                            if (not cut_time or (step_time >= start and (stop is float('inf') or step_time <= stop))):

                                typLev = gr.typeOfLevel

                                if (level_type is None or typLev in level_type):

                                    if self.verbose:
                                        logging.debug(f'''centre: {gr.centre}
                                        dataDate: {gr.dataDate}
                                        dataTime: {gr.dataTime}
                                        step: {gr.step}
                                        shortName: {gr.shortName}
                                        paramId: {gr.paramId}
                                        name: {gr.name}
                                        typeOfLevel: {gr.typeOfLevel}
                                        level: {gr.level}
                                        data shape: {gr.values.shape}
                                        perturbationNumber: {gr.get('perturbationNumber')} 
                                        gridType: {gr.gridType}
                                        projparams: {gr.projparams}
                                        ''')

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
                                            unit_time_range = gr.get('unitOfTimeRange', gr.get('stepUnits', 255))
                                            data.update({'time_units': self.unitOfTimeRange[unit_time_range]})
                                            self.time_units = '{0} since {1}'.format(
                                                self.unitOfTimeRange[unit_time_range], ref_time)

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

                                        if flip_lat:  # error with lat/lon 2 dims arrays
                                            self.lat = np.flip(self.lat, axis=0)

                                        # select spatial subdomain .......
                                        y, x = [None, None], [None, None]

                                        if cut_domain:

                                            if isinstance(cut_domain, (tuple, list)):
                                                lat1, lon1, lat2, lon2 = cut_domain
                                                while True:  # necessary 2 pass to fix 360 - 0 descontinuity
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

                                        # trim lat/lon dimensions .........
                                        self.lat = self.lat[y[0]:y[1], x[0]:x[1]]
                                        self.lon = self.lon[y[0]:y[1], x[0]:x[1]]

                                        # if necessary roll longitude due discontinuity 360-0 of the longitude
                                        gr.values = np.roll(gr.values, cut_domain_roll, axis=-1)

                                        # get data ........................
                                        # grab data and flip the latitude axis if necessary
                                        if flip_lat:
                                            _tmp = np.flip(gr.values, axis=0)[None, None, None, y[0]:y[1], x[0]:x[1]]
                                        else:
                                            _tmp = gr.values[None, None, None, y[0]:y[1], x[0]:x[1]]

                                        if idVar in _data.keys() and typLev in _data[idVar].keys():

                                            member_num = gr.get('perturbationNumber')

                                            # concatenate levels
                                            _data[idVar][typLev].value = np.concatenate((_data[idVar][typLev].value, _tmp),
                                                                                        axis=2)

                                            if gr.level not in _data[idVar][typLev].level:
                                                _data[idVar][typLev].level.append(gr.level)

                                            if member_num not in _data[idVar][typLev].members:
                                                _data[idVar][typLev].members.append(gr.perturbationNumber)

                                        else:
                                            member_num = gr.get('perturbationNumber', 0)

                                            # to add a new variable update the
                                            # "self.__non_data_variables" list

                                            __tmp = {
                                                typLev: {'value': _tmp,
                                                         'level': [gr.level],
                                                         'members': [member_num]},
                                                'centre': gr.centre,
                                                'dataType': gr.dataType,
                                                'param_id': gr.paramId,
                                                'long_name': gr.name,
                                                'parameter_units': gr.parameterUnits,
                                                'latitude': self.lat,
                                                'longitude': self.lon,
                                                'grid_type': gr.gridType,
                                                'projparams': gr.projparams
                                            }


                                            if idVar in _data.keys():
                                                _data[idVar].update(__tmp)
                                                _data[idVar].level_type.append(typLev)
                                            else:
                                                _data[idVar] = __tmp
                                                _data[idVar].level_type = [typLev]

                        # consolidate data for last time block or stop time ................
                        if n + 1 == msg_len:

                            data = self.__concat_time(_data, data)

                            self.variables = list(data.keys())
                            self.coordinates.append('latitude')
                            self.coordinates.append('longitude')
                            self.coordinates.append('level')
                            self.coordinates.append('members')

                    # rearrange data (set member axis, sort levels and member data)
                    data = self.__arrange_data(data)

        except Exception as e:
            logging.exception(f'gdio.gb_load: {e}')

        return data

    def __arrange_data(self, data):
        '''
        rearrange data (set member axis, sort levels and member data)
        :param data:    dict
                        data dictionary
        :return:        dict
                        data dictionary
        '''

        for k, v in data.items():
            if isinstance(v, dict):
                try:
                    for l in v.keys():
                        if not l in self.__non_data_variables:
                            dims = list(data[k][l].value.shape)
                            levels = data[k][l].level
                            members = data[k][l].members
                            dims[0], dims[2] = len(members), len(levels)

                            # set member to dimension 0
                            data[k][l].value = data[k][l].value.reshape(dims)

                            # sort levels
                            data[k][l].value = data[k][l].value[:, :, np.argsort(levels)]
                            data[k][l].level = sorted(levels)

                            # sort levels
                            data[k][l].value = data[k][l].value[np.argsort(members)]
                            data[k][l].members = sorted(members)
                except Exception as e:
                    logging.exception(f'gdio.__arrange_data: {e}')
        return data

    def gb_write(self,
                 ofile,
                 data,
                 packingType='grid_simple',
                 least_significant_digit=3,
                 **kwargs) -> None:
        '''
        Write grib file

        :param ifile:           string
                                file path
        :param data:            dict
                                dataset
        :param packingType:     string
                                packingType	Type of packing:
                                    grid_simple
                                    spectral_simple
                                    grid_simple_matrix
                                    grid_jpeg
                                    grid_png
                                    grid_ieee
                                    grid_simple_log_preprocessing
                                    grid_second_order
        :param least_significant_digit: int (default None)
                                        specify the power of ten of the smallest decimal place in the data that is a
                                        reliable value that dramatically improve the compression by quantizing
                                        (or truncating) the data
        :param kwargs:              key-value parameter
                                    additional grib key: edition, editionNumber, centre, subCentre, discipline,
                                                         dataType, missingValue

        :return:
        '''

        data = data if isinstance(data, objectify) else objectify(data)

        step_type = data.time_units

        if isinstance(data.time_units, str):
            _, step_type = dict_get(UNIT_TIME_RANGE, key=data.time_units)

        # convert timestep to datetime if necessary
        if isinstance(data.time[0], (int, np.int64)):
            time = data.ref_time + timestep_to_datetime(data.time, units=self.__unity(step_type))
        else:
            time = data.time

        # step calculation
        dt = 0
        tshift = int((time[0] - data.ref_time).total_seconds() / (3600 * self.__unity(step_type)))

        if len(time) > 1:
            dt = int((time[1] - time[0]).total_seconds() / (3600 * self.__unity(step_type)))

        grb = cgrib.fwrite(filename=ofile)

        for t, timestep in enumerate(time):

            for idVar in self.__get_vars(data):

                data_type = kwargs.get('dataType', 'fc') if kwargs.get('dataType') else data[idVar].dataType
                timestep = data.ref_time if data_type in ['fc'] else timestep

                # convert lat lon to 2d mesh coordinates
                if data[idVar].latitude.ndim == 1 and data[idVar].longitude.ndim == 1:
                    dims = (data[idVar].longitude.size, data[idVar].latitude.size)

                    data[idVar].latitude = np.tile(data[idVar].latitude,
                                                   (dims[0], 1)
                                                   ).T
                    data[idVar].longitude = np.tile(data[idVar].longitude,
                                                    (dims[1], 1))

                msg = {
                        'edition': kwargs.get('edition', 2),
                        'editionNumber': kwargs.get('editionNumber', 2),
                        'centre': kwargs.get('centre', data[idVar].centre),
                        'subCentre': kwargs.get('subCentre', 0),
                        'discipline': kwargs.get('discipline', 0),
                        'dataType': data_type,

                        # time namespace
                        'dataDate': int(f'{timestep:%Y%m%d}'),
                        'dataTime': int(f'{timestep:%H%M}'),
                        'stepUnits': step_type,
                        'step': t * dt + tshift,

                        # variable namespace
                        'paramId': data[idVar].param_id,
                        'shortName': idVar,
                        'missingValue': kwargs.get('missingValue', 99999),
                        'packingType': packingType,
                        'changeDecimalPrecision': least_significant_digit,

                        # projection namespace
                        'gridType': data[idVar].grid_type,
                        'latitude': data[idVar].latitude,
                        'longitude': data[idVar].longitude
                        }

                for level_type in data[idVar].get('level_type', ['surface']):

                    msg.update({
                                'typeOfLevel': level_type
                                })

                    if not level_type in self.__non_data_variables:

                        dims = data[idVar][level_type].value.shape

                        # level
                        for l in range(dims[2]):

                            msg.update({'level': data[idVar][level_type].level[l]})

                            # ensemble loop
                            for m in range(dims[0]):
                                msg.update({'value': data[idVar][level_type].value[m,t,l].squeeze()})
                                grb.write(message=msg)


        grb.close()



    def __get_vars(self, data):
        '''
        Extract variables from data keys
        :param data:    dict
                        data dictionary
        :return:        list
                        list of data variables
        '''
        b = self.__fields_latitude + self.__fields_longitude + self.__fields_time
        return list(set(data.keys()) - set(b))

    def fcstTime(self, gr):
        '''
        Convert ref time + time step to forecast time
        :param gr:       object
                         grib message object
        '''

        if gr.step > 0:
            return datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute) + timedelta(hours=gr.step * self.__unity(gr.stepUnits))
        else:
            return datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)

    def __unity(self, stepUnits=1):
        '''
        Scale fator to time unity transformation
        :param stepUnits:  int
                           pygrig object
        :return:           float
        '''

        scale = 1.0

        if stepUnits == 0:  # minute
            scale = 1 / 60
        elif stepUnits == 1:  # hour
            scale = 1
        elif stepUnits == 2:  # day
            scale = 24
        elif stepUnits == 3:  # month   # problem
            scale = 24 * 30
        elif stepUnits == 4:  # year    # problem
            scale = 24 * 365
        elif stepUnits == 5:  # decade
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

        for k, v in _data.items():
            if k in __data.keys():
                for l in v.keys():
                    if not l in self.__non_data_variables:
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
