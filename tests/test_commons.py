from gdio.commons import near_yx2 as near_yx
import os
import sys
import numpy as np
import pyproj
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

class TestGribFiles(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass


    def setUp(self):

        self.nx = 111
        self.ny = 101

        self.coords =(-40, 270, 20, 330)



    def test_nearyx_regular_single_point(self):

        lons, lats = self.set_projection(proj='cyl')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=0, lons=290),
                         ([67], [37]),
                         'incorrect number of variables')

    def test_nearyx_regular_band(self):

        lons, lats = self.set_projection(proj='cyl')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=0),
                         ([67], [None]),
                         'incorrect results')


    def test_nearyx_regular_list_points(self):

        lons, lats = self.set_projection(proj='cyl')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=[10, -5], lons=[280, None]),
                              ([83, 58], [18, None]),
                              'incorrect results')

    def test_nearyx_mercator_single_point(self):

        lons, lats = self.set_projection(proj='cyl')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=0, lons=290),
                         ([67], [37]),
                         'incorrect results')

    def test_nearyx_mercator_band(self):

        lons, lats = self.set_projection(proj='cyl')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lons=290),
                              ([None], [37]),
                              'incorrect results')

    def test_nearyx_mercator_list_points(self):

        lons, lats = self.set_projection(proj='merc')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=[10, -5], lons=[280]),
                              ([84, 60], [18, None]),
                              'incorrect number of variables')


    def test_nearyx_lambert_single_point(self):

        lons, lats = self.set_projection(proj='lcc')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=0, lons=290),
                              ([93], [36]),
                              'incorrect number of variables')

    def test_nearyx_lambert_list_points(self):

        lons, lats = self.set_projection(proj='lcc')
        self.assertTupleEqual(near_yx({'latitude': lats, 'longitude': lons}, lats=[10, -5], lons=[280]),
                              ([None, 82], [17, None]),
                              'incorrect number of variables')




    def set_projection(self, proj='cyl'):
        '''
        create regular, mercator and lambert sintetic coordinates mesh
        :param proj:    string
                        kind of projection: cly  - regular lat-lon
                                            merc - mercator
                                            lcc - lambert
        :return: longitude and latitude mesh array
        '''

        projparams = {}
        lat1, lon1, lat2, lon2 = self.coords
        projparams['a'] = 6371229.0
        projparams['b'] = 6371229.0

        # regular
        projparams['proj'] = proj

        if projparams['proj'] in ['merc']:

            # mercator
            projparams['lat_ts'] = 0
            projparams['lon_0'] = 0.5 * (lon1 + lon2)
            pj = pyproj.Proj(projparams)
            llcrnrx, llcrnry = pj(lon1, lat1)
            urcrnrx, urcrnry = pj(lon2, lat2)
            dx = (urcrnrx - llcrnrx) / (self.nx - 1)
            dy = (urcrnry - llcrnry) / (self.ny - 1)
            x = llcrnrx + dx * np.arange(self.nx)
            y = llcrnry + dy * np.arange(self.ny)

        elif projparams['proj'] in ['lcc']:

            # lambert
            projparams['lon_0'] = lon1 + (lon2 - lon1) / 2  # LoVInDegrees
            projparams['lat_0'] = lat1 + (lat2 - lat1) / 2  # LaDInDegrees
            projparams['lat_1'] = lat1 + (lat2 - lat1) / 3  # Latin1InDegrees
            projparams['lat_2'] = lat2 - (lat2 - lat1) / 3  # Latin2InDegrees
            dx = 50000  # DxInMetres
            dy = 50000  # DyInMetres
            pj = pyproj.Proj(projparams)
            llcrnrx, llcrnry = pj(lon1, lat1)
            x = llcrnrx + dx * np.arange(self.nx)
            y = llcrnry + dy * np.arange(self.ny)

        elif projparams['proj'] in ['cyl']:
            dx = (lon2 - lon1) / (self.nx - 1)
            dy = (lat2 - lat1) / (self.ny - 1)

            y = np.arange(lat1, lat2, dy)
            x = np.arange(lon1, lon2, dx)

        lons, lats = np.meshgrid(x, y)

        if projparams['proj'] not in ['cyl']:
            lons, lats = pj(lons, lats, inverse=True)

        return lons, lats


if __name__ == '__main__':
    unittest.main()