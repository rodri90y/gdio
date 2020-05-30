import numpy as np


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