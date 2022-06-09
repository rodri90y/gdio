import difflib
import os
import numpy as np
import itertools

from datetime import datetime, timedelta


class objectify(dict):
    """ Nested Attribute Dictionary
        A class to convert a nested Dictionary into an object with key-values
        accessible using attribute notation (AttrDict.attribute) in addition to
        key notation (Dict["key"]). This class recursively sets Dicts to objects,
        allowing you to recurse into nested dicts (like: AttrDict.attr.attr)
        """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for key, value in self.items():
            if not isinstance(value, (dict, list, tuple, set)):
                continue

            super().__setitem__(key, self.__set(value))

    def __set(self, value):
        if isinstance(value, self.__class__):
            return value
        if isinstance(value, dict):
            return self.__class__(value)
        if isinstance(value, (list, tuple, set)):
            return type(value)(map(self.__set, value))
        return value

    def __delattr__(self, name):
        return self.__delitem__(name)

    def __getattr__(self, name):
        return self.__getitem__(name)

    def __setattr__(self, name, value):
        return self.__setitem__(name, value)

    def __setitem__(self, key, value):
        value = self.__set(value)
        return super().__setitem__(key, value)

    def __getitem__(self, key):
        '''return {} for missing keys'''
        if key in self:
            return super().__getitem__(key)
        # else:
        #     logging.error(f'Missing key: {key}')
        #     return self.__class__({})

    # tell pickle to treat your class just like a normal one
    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)
    # .......................................................

    def copy(self):
        return self.__class__(super().copy())

    def update(self, *args, **kwargs):
        other = self.__class__(*args, **kwargs)
        return super().update(other)


def timestep_to_datetime(ts, units=1):
    '''
    Convert timestep index to timeserie
    :param ts:          int list
                        timestep list
    :param unity:       int
                        unit of time factor at hours
    :return:            datetime array
    '''

    __convert = np.vectorize(
        lambda t, u: timedelta(hours=float(t * u)) if not isinstance(t, datetime) else t
    )

    return __convert(ts, units)

def datetime_to_timestep(ts, ref, units=1):
    '''
    Convert timeserie index to timestep
    :param ts:          datetime array
    :param ref:         reference datetime
    :param unity:       int
                        unit of time factor at hours
    :return:            int list
                        timestep list
    '''

    __convert = np.vectorize(
        lambda t, u: t.total_seconds() / 3600 * u
    )

    return __convert(ts-ref, units)

def near_yx(data, lats=None, lons=None):
    '''
    Find the nearst coordinate i/j given a lat/lon coordinate
    Warning error with lat/lon parameter with dims>1, only works with
    mercator, regular lat-lon
    :param lat:     float
                    latitude
    :param lon:     float
                    longitude
    :return:        float tuple lat, lon
                    nearest grid point y/x
    '''

    lats = lats if isinstance(lats, list) else [lats]
    lons = lons if isinstance(lons, list) else [lons]

    lons = [(_lo + 360) % 360 if not _lo is None else _lo for _lo in lons]

    x = []
    y = []

    if 'lat' in data.keys() or \
            'lon' in data.keys():

        _lat = data['lat']
        _lon = data['lon']

    elif 'latitude' in data.keys() or \
            'longitude' in data.keys():

        _lat = data['latitude']
        _lon = data['longitude']

    # convert -180,180 to 0,360 format
    _lon = (_lon + 360) % 360

    for lat, lon in itertools.zip_longest(lats, lons):

        _y = None
        _x = None

        if not lat is None:
            if np.min(_lat) <= lat and np.max(_lat) >= lat:
                _y = np.nanargmin(np.abs(_lat - lat)) if lat is not None else lat

        if not lon is None:
            if np.min(_lon) <= lon and np.max(_lon) >= lon:
                _x = np.nanargmin(np.abs(_lon - lon)) if lon is not None else lon

        x.append(_x)
        y.append(_y)

    return y, x

def near_yx2(data, lats=None, lons=None):
    '''
    Find the nearst coordinate i/j given a lat/lon coordinate
    Warning error with lat/lon parameter with dims>1, only works with
    mercator, regular lat-lon
    :param data:    dict
                    latitude and longitude mesh data
    :param lat:     float
                    latitude
    :param lon:     float
                    longitude
    :return:        float tuple lat, lon
                    nearest grid point y/x
    '''

    lats = lats if isinstance(lats, list) else [lats]
    lons = lons if isinstance(lons, list) else [lons]

    lons = [(_lo + 360) % 360 if not _lo is None else _lo for _lo in lons]

    x = []
    y = []

    if 'lat' in data.keys() or \
            'lon' in data.keys():

        _lat = data['lat']
        _lon = data['lon']

    elif 'latitude' in data.keys() or \
            'longitude' in data.keys():

        _lat = data['latitude']
        _lon = data['longitude']

    # convert -180,180 to 0,360 format
    _lon = (_lon + 360) % 360


    for lat, lon in itertools.zip_longest(lats, lons):

        _x = np.zeros(_lon.shape)
        _y = np.zeros(_lat.shape)

        if not lat is None:
            if np.min(_lat) <= lat and np.max(_lat) >= lat:
                _y = np.abs(_lat - lat)
            else:
                lat = None

        if not lon is None:
            if np.min(_lon) <= lon and np.max(_lon) >= lon:
                _x = np.abs(_lon - lon)
            else:
                lon = None


        xy_min = np.nanargmin((_y + _x))

        #convert the 1D index to 2D index system
        lat_index = xy_min // _lon.shape[1]
        lon_index = xy_min - (lat_index * _lon.shape[1])

        lat_index = None if lat is None else lat_index
        lon_index = None if lon is None else lon_index

        # print(lat, lon)
        x.append(lon_index)
        y.append(lat_index)

    return y, x




def dict_get(data, key=None, by='values'):
    '''
    Get closest matches of list keys, return
    list of the best "good enough" matches
    :param data:    dict

    :param key:     str
    :param by:      str
    :return:        list
                    lists or keys criteria
    '''
    if by == 'values':
        for k, v in data.items():
            v = v if isinstance(v, list) else [v]
            _k = difflib.get_close_matches(key, v, 1, 0.85)

            if _k:
                return _k, k
    else:
        return data.get(key)


def show_data_structure(data):
    '''
    Describe the data structure
    :param data:    dict
    return:
    '''

    for _dat in data:
        __data_tree(_dat, depth=0)



def __data_tree(data, depth=0):
    '''
    Get data structure tree
    :param data:    dict
    :param depth:   int
    '''
    depth += 1
    for k, v in data.items():

        if isinstance(v, dict):
            print("{s:<{depth}}+-- {key} ".format(s="", depth=depth * 4, key=k)
                  )

            __data_tree(v, depth)
        else:
            if isinstance(v, np.ndarray):

                print("{s:<{depth}}+-- {key}: {type} {shape}".format(s="",
                                                                    depth=depth * 4,
                                                                    key=k,
                                                                    type=type(v),
                                                                    shape=v.shape)
                      )

            else:
                print("{s:<{depth}}+-- {key}: {value}".format(s="",
                                                             depth=depth * 4,
                                                             key=k,
                                                             value=str(v))
                     )

