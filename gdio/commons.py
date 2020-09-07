import numpy as np
import difflib

class objectify(dict):
    """ Nested Attribute Dictionary
        A class to convert a nested Dictionary into an object with key-values
        accessible using attribute notation (AttrDict.attribute) in addition to
        key notation (Dict["key"]). This class recursively sets Dicts to objects,
        allowing you to recurse into nested dicts (like: AttrDict.attr.attr)
        """

    def __init__(self, mapping=None):
        super(objectify, self).__init__()
        if mapping is not None:
            for key, value in mapping.items():
                self.__setitem__(key, value)

    def __setitem__(self, key, value):
        if isinstance(value, dict):
            value = objectify(value)
        super(objectify, self).__setitem__(key, value)
        self.__dict__[key] = value  # for code completion in editors

    def __getattr__(self, item):
        try:
            return self.__getitem__(item)
        except KeyError:
            raise AttributeError(item)

    __setattr__ = __setitem__

def near_yx(data, lats=None, lons=None):
    '''
    Find the nearst coordinate i/j given a lat/lon coordinate
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

    for lat in list(lats):
        _y = None
        if not lat is None:
            if np.min(_lat) <= lat and np.max(_lat) >= lat:
                _y = np.nanargmin(np.abs(_lat - lat)) if lat is not None else lat
        y.append(_y);

    for lon in list(lons):
        _x = None
        if not lon is None:
            if np.min(_lon) <= lon and np.max(_lon) >= lon:
                _x = np.nanargmin(np.abs(_lon - lon)) if lon is not None else lon
        x.append(_x);

    return y, x


def dict_get(data, key=None, by='values'):
    '''
    Get closest matches of list keys
    :param data:    dict

    :param key:     str
    :param by:      str
                    lists or keys criteria
    :return:
    '''
    if by=='values':
        for k, v in data.items():
            v = v if isinstance(v, list) else list(v)
            _k = difflib.get_close_matches(key, v, 1, 0.85)
            if _k:
                return _k, k
    else:
        return data.get(key)
