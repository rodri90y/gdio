import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gdio.netcdf import netcdf
from gdio.commons import near_yx


class TestNcFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        nc = netcdf()

        path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), ".")),
                            'data/era5_20191227_lev.nc')

        self.nc = nc.nc_load(path,
                             cut_domain=(-30, -60, 10, -40),
                             cut_time=(12, 24))

    def setUp(self):

        self.expected_dim = (1, 1, 7, 80, 40)
        self.expected_times = [12]
        self.expected_level_type = ['isobaricInhPa']
        self.expected_units = 'm s**-1'
        self.expected_coordinate = ([13], [27])
        self.expected_levels = [200, 300, 500, 700, 800, 950, 1000]


    def test_open_netcdf(self):
        self.assertTrue(not self.nc is {})

    def test_netcdf_variables_test(self):

        self.assertEqual(list(self.nc.keys()),
                         ['ref_time', 'time_units', 'time', 'r', 't', 'u', 'v'],
                         'incorrect number of variables')

    def test_netcdf_varible_dimension(self):
        self.assertEqual(self.nc.get('u').isobaricInhPa.value.shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_grib_levels(self):
        self.assertEqual(self.nc.get('u').isobaricInhPa.level, self.expected_levels,
                         'levels of u variable incorrect')

    def test_grib_level_type(self):
        self.assertEqual(self.nc.get('u').level_type, self.expected_level_type,
                         'level type of u variable incorrect')

    def test_grib_varible_units(self):

        self.assertEqual(self.nc.get('u').parameter_units, self.expected_units,
                         'units of u variable incorrect')

    def test_grib_cut_time(self):

        self.assertListEqual(list(self.nc.get('time')), self.expected_times,
                         'incorrect time cut')

    def test_grib_cut_space(self):
        self.assertEqual(near_yx({'latitude': self.nc.get('u').latitude,
                                  'longitude': self.nc.get('u').longitude},
                                 lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')


if __name__ == '__main__':
    unittest.main()