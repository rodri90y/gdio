__author__ = "Rodrigo Yamamoto"
__date__ = "2026.Ago"
__credits__ = ["Rodrigo Yamamoto", "Igor Santos"]
__maintainer__ = "Rodrigo Yamamoto"
__email__ = "codes@rodrigoyamamoto.com"
__version__ = "version 0.3.5"
__license__ = "MIT"
__status__ = "development"
__description__ = "A grib file IO library"

import gc
import logging
from datetime import datetime, timedelta
import numpy as np

from gdio import cgrib
from gdio.commons import near_yx2, objectify, dict_get, timestep_to_datetime
from .definitions.Table_4_4 import UNIT_TIME_RANGE

class grib(object):

    def __init__(self, verbose=False, debug=False):

        self.verbose = verbose
        self.debug = debug

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
                                         'projparams',
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
                filter_by=None,
                rename_vars=None,
                **kwargs):
        '''
        Load grib file
        Yamamoto, R @ Ago.2026
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

        data = objectify()

        filter_by = {} if filter_by is None else filter_by
        rename_vars = {} if rename_vars is None else rename_vars
        vars = vars if vars is None else list(vars)
        cut_time = cut_time if cut_time is None else tuple(cut_time)
        cut_domain = cut_domain if cut_domain is None else tuple(cut_domain)
        level_type = level_type if level_type is None else list(level_type)

        # Get grib file metadata
        meta = self.__get_metadata(ifile, vars, level_type, cut_time, filter_by, rename_vars)

        if not meta:
            return data

        ref_time = meta['ref_time']
        members_set = meta['members_set']
        times_set = meta['times_set']
        levels_by_var = meta['levels_by_var']
        var_meta = meta['var_meta']
        unit_time_range = meta['unit_time_range']

        # Sort the metadata (Members and Date/Time/Step)
        sorted_members = sorted(list(members_set))
        sorted_times = sorted(list(times_set))

        idx_member = {m: i for i, m in enumerate(sorted_members)}
        idx_time = {t: i for i, t in enumerate(sorted_times)}
        idx_level = {}

        gr = next(iter(var_meta.values()))

        # set spatial coordinates ......
        self.lat, self.lon = gr.latlons()
        # convert from -180,180 to 360 format
        self.lon = (self.lon + 360) % 360
        flip_lat = self.lat[-1, 0] < self.lat[0, 0]

        if flip_lat:
            self.lat = np.flip(self.lat, axis=0)

        # select spatial subdomain .......
        y, x = [None, None], [None, None]

        cut_domain_roll = 0

        if cut_domain and isinstance(cut_domain, (tuple, list)):
            lat1, lon1, lat2, lon2 = cut_domain
            while True: # necessary 2 pass to fix 360 - 0 descontinuity
                y, x = near_yx2({'latitude': self.lat, 'longitude': self.lon}, lats=[lat1, lat2],
                                lons=[lon1, lon2])
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
        xul = x[1] if x[1] is None else x[1] + 1
        yul = y[1] if y[1] is None else y[1] + 1
        self.lat = self.lat[y[0]:yul, x[0]:xul]
        self.lon = self.lon[y[0]:yul, x[0]:xul]
        shape_2d = self.lat.shape


        data.update({'ref_time': ref_time})
        data.update({'time_units': self.unitOfTimeRange[unit_time_range]})
        self.time_units = f"{self.unitOfTimeRange[unit_time_range]} since {ref_time}"
        data.update({'time': np.array(sorted_times)})
        self.time = data['time']

        # Pre-allocation of the data matrix .................
        for (idVar, typLev), lev_set in levels_by_var.items():
            s_levels = sorted(list(lev_set))
            idx_level[(idVar, typLev)] = {l: i for i, l in enumerate(s_levels)}

            gr_ref = var_meta[(idVar, typLev)]
            shape_5d = (len(sorted_members), len(sorted_times), len(s_levels), shape_2d[0], shape_2d[1])

            __tmp = {
                typLev: {
                    'value': np.empty(shape_5d, dtype=np.float32),
                    'level': s_levels,
                    'members': sorted_members
                },
                'centre': gr_ref.centre,
                'dataType': gr_ref.dataType,
                'param_id': gr_ref.paramId,
                'long_name': gr_ref.name,
                'parameter_units': gr_ref.parameterUnits,
                'latitude': self.lat,
                'longitude': self.lon,
                'grid_type': gr_ref.gridType,
                'projparams': gr_ref.projparams
            }

            if idVar in data.keys():
                data[idVar].update(__tmp)
                data[idVar].level_type.append(typLev)
            else:
                data[idVar] = __tmp
                data[idVar].level_type = [typLev]

        # Filling data .................................
        try:
            #
            del var_meta, levels_by_var, meta
            gc.collect()

            with cgrib.fopen(ifile) as msg:
                for gr in msg:

                    typLev = gr.typeOfLevel
                    idVar = gr.shortName if gr.shortName not in ['', 'unknown'] else str(gr.paramId)
                    for k, v in rename_vars.items():
                        if idVar in k:
                            idVar = v

                    key = (idVar, typLev)

                    # Check the idx_level key existence
                    if key in idx_level:
                        f_date = self.fcstTime(gr)
                        fcst_time = self.__forecast_step(f_date, ref_time, gr.stepUnits)

                        # Check if the GRIB timestep is in the mapped dictionary
                        if fcst_time in idx_time:
                            member_num = self.__get_member(gr)

                            # skip members not selected
                            if member_num not in idx_member:
                                if hasattr(gr, 'gid'):
                                    cgrib.eccodes.codes_release(gr.gid)
                                continue

                            m_idx = idx_member[member_num]

                            t_idx = idx_time[fcst_time]
                            l_idx = idx_level[key][gr.level]

                            # if necessary roll longitude due discontinuity 360-0 of the longitude
                            index = (m_idx, t_idx, l_idx, y, yul, x, xul)
                            self.__slice_and_assign(data[idVar][typLev].value, gr.values, index, cut_domain_roll, flip_lat)

                            if hasattr(gr, 'gid'):
                                cgrib.eccodes.codes_release(gr.gid)

                self.variables = list(data.keys())
                self.coordinates = ['latitude', 'longitude', 'level', 'members']

        except Exception as e:
            logging.exception(f'gdio.gb_load filling data: {e}')

        return data

    def __slice_and_assign(self, data, gr, index, cut_domain_roll, flip_lat):
        '''
            Roll and flip lat lon grid dimensions
            Yamamoto, R @ Ago.2026
            :param data:            numpy
                                    data array
            :param gr:              object
                                    grib message object
            :param index:           tuple
                                    index parameters
            :param cut_domain_roll: bool
                                    enable roll the longitude
            :param flip_lat:        bool
                                    flip the latitude axis
        '''
        m_idx, t_idx, l_idx, y, yul, x, xul = index

        if cut_domain_roll != 0:
            split_idx = -cut_domain_roll
            vals_right = gr[y[0]:yul, split_idx:]
            vals_left = gr[y[0]:yul, :split_idx]

            if flip_lat:
                data[m_idx, t_idx, l_idx, :, :split_idx] = np.flip(vals_right, axis=0)
                data[m_idx, t_idx, l_idx, :, split_idx:] = np.flip(vals_left, axis=0)
            else:
                data[m_idx, t_idx, l_idx, :, :split_idx] = vals_right
                data[m_idx, t_idx, l_idx, :, split_idx:] = vals_left
        else:
            if flip_lat:
                data[m_idx, t_idx, l_idx, :, :] = np.flip(gr[y[0]:yul, x[0]:xul], axis=0)
            else:
                data[m_idx, t_idx, l_idx, :, :] = gr[y[0]:yul, x[0]:xul]



    def __get_metadata(self, ifile, vars, level_type, cut_time, filter_by, rename_vars):
        '''
                Mapping the structure of the grib file
                Yamamoto, R @ Ago.2026
                :param ifile:       string
                                    grib 1 or 2 file name
                :param vars:        list
                                    variables short name or id parameter number
                :param level_type:  list
                                    type of level (hybrid, isobaricInhPa, surface)
                :param cut_time:    tuple
                                    range of time to cut ex.: (0,10)/(0,None)/(None,10)
                :param filter_by:   dictonary
                                    dict with grib parameters at form of pair key:values (list or single values)
                                    eg: filter_by={'perturbationNumber': [0,10],'level': [1000,500,250]}
                                    or filter_by={'gridType': 'regular_ll'}
                :param rename_vars: dictonary
                                    rename variables names (key) for a new name (value).
                                    Eg. {'tmpmdl': 't', 'tmpprs': 't'}
                :return:            dictonary/attributes
                                    multiple data container
                '''

        ref_time = None
        forecastDate = None
        step_time = -1
        member_num = 0
        times_set = set()
        members_set = set()
        levels_by_var = {}
        var_meta = {}
        unit_time_range = 255


        # set temporal subdomain .......
        start, stop = 0, float('inf')
        if isinstance(cut_time, (tuple, list)) and len(cut_time) == 2:
            start = 0 if cut_time[0] is None else cut_time[0]
            stop = float('inf') if cut_time[1] is None else cut_time[1]

        try:
            with cgrib.fopen(ifile) as msg:

                for gr in msg:
                    if all([gr[k] in (v if isinstance(v, list) else [v]) for k, v in filter_by.items() if k in gr.keys()]):

                        # initialize time
                        if ref_time is None:
                            ref_time = datetime(gr.year, gr.month, gr.day, gr.hour, gr.minute)
                            self.history = f"Created by gdio @ {datetime.now():%Y%m%d%H}"
                            self.grid_description = {k: v for k, v in gr.items() if k in gr.gridkeys}
                            self.centre = gr.centre
                            unit_time_range = gr.get('unitOfTimeRange', gr.get('stepUnits', 255))

                        # set time coordinate ....................
                        if forecastDate != self.fcstTime(gr):
                            forecastDate = self.fcstTime(gr)
                            step_time += 1

                        # set  member coordinate ..................
                        if member_num != self.__get_member(gr):
                            step_time = 0
                            member_num = self.__get_member(gr)

                        # cut time between start and stop time
                        if not cut_time or (start <= step_time <= stop):

                            typLev = gr.typeOfLevel

                            if self.debug:
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

                            if level_type is None or typLev in level_type:
                                idVar = gr.shortName if gr.shortName not in ['', 'unknown'] else str(gr.paramId)

                                # rename variables
                                for k, v in rename_vars.items():
                                    if idVar in k:
                                        idVar = v

                                if vars is None or gr.shortName in vars or gr.paramId in vars:

                                    fcst_time = self.__forecast_step(forecastDate, ref_time, gr.stepUnits)

                                    members_set.add(member_num)
                                    times_set.add(fcst_time)

                                    key = (idVar, typLev)
                                    if key not in levels_by_var:
                                        levels_by_var[key] = set()
                                        var_meta[key] = gr
                                    levels_by_var[key].add(gr.level)
        except Exception as e:
            logging.exception(f'gdio.__scan_metadata get metadata: {e}')

        if not levels_by_var:
            logging.warning(
                f"gdio.__scan_metadata: No data was found on '{ifile}' "
            )
            return None

        return {
            'ref_time': ref_time,
            'members_set': members_set,
            'times_set': times_set,
            'levels_by_var': levels_by_var,
            'var_meta': var_meta,
            'unit_time_range': unit_time_range
        }


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

        dims = self.__get_dims(data)

        # ensemble loop
        for m in range(dims[0]):
            # time loop
            for t, timestep in enumerate(time):
                # variables loop
                for idVar in self.__get_vars(data):
                    data_type = kwargs.get('dataType', 'fc') if kwargs.get('dataType') else data[idVar].dataType

                    # convert lat lon to 2d mesh coordinates
                    #     #TODO: reformular para nao gravar lat/lon como matriz e sim como parametros
                    if data[idVar].latitude.ndim == 1 and data[idVar].longitude.ndim == 1:
                        dims = (data[idVar].longitude.size, data[idVar].latitude.size)

                        data[idVar].latitude = np.tile(data[idVar].latitude,
                                                       (dims[0], 1)
                                                       ).T
                        data[idVar].longitude = np.tile(data[idVar].longitude,
                                                        (dims[1], 1))

                    timestep = data.ref_time if data_type in ['fc'] else timestep

                    for level_type in data[idVar].get('level_type', ['surface']):
                        if not level_type in self.__non_data_variables:

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
                                'latitude': data[idVar].latitude,  # esta gravando 2d em cada mensagem,
                                'longitude': data[idVar].longitude,
                                'typeOfLevel': level_type
                            }


                            # level (actual level)
                            for l in range(dims[2]):

                                msg.update({'level': data[idVar][level_type].level[l]})

                                if dims[0]>1:
                                    msg.update({
                                        # Enable GRIB2 (Template 4.1)
                                        'productDefinitionTemplateNumber': 1,
                                        'typeOfEnsembleForecast': 0 if m == 0 else 1,
                                        'perturbationNumber': m,
                                        'numberOfForecastsInEnsemble': dims[0],
                                        # MARS keys compatibility
                                        'marsType': 'cf' if m == 0 else 'pf',
                                        'number': m,
                                    })
                                else:
                                    msg.update({
                                        'productDefinitionTemplateNumber': 0,
                                    })
                                try:
                                    msg.update({'value': data[idVar][level_type].value[m,t,l].squeeze()})
                                except:
                                    logging.exception(f'gb_write error: {m}, {t}, {l}, {data[idVar][level_type].value.shape}')

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

    def __get_dims(self, data):

        for k, v in data.items():
            if isinstance(v, dict):
                try:
                    for l in v.get('level_type'):
                          return v[l].value.shape
                except (KeyError, AttributeError):
                    pass


    def __get_member(self, gr):
        '''
        Extract member from data keys
        :param gr:       object
                         grib message object
        :return:        int
                        member number
        '''

        return gr.get('perturbationNumber') or gr.get('number') or 0

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

    def __forecast_step(self, forecast_date, ref_time, step_units):
        return int((forecast_date - ref_time).total_seconds() / (self.__unity(step_units) * 3600))


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
