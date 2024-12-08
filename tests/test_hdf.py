from gdio.commons import near_yx2
from gdio.hdf import hdf
import os
import sys
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "")))


class TestNcFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), "."))

        hd = hdf()
        self.hd = hd.hdf_load(os.path.join(root, 'data/era5_2019122712_lev.hdf'),
                              cut_domain=(-30, -60, 10, -40),
                              cut_time=(1, 1),
                              rename_vars={'t': 't2m'}
                             )

        # write new hdf
        hd.hdf_write(os.path.join(root, 'tmp.hdf'), self.hd)

        # open new hdf
        self.new_hd = hd.hdf_load(os.path.join(root, 'tmp.hdf'))

        os.remove(os.path.join(root, 'tmp.hdf'))


    def setUp(self):
        self.expected_dim = (1, 1, 7, 80, 40)
        self.expected_times = [12]
        self.expected_level_type = ['isobaricInhPa']
        self.expected_units = 'm s**-1'
        self.expected_coordinate = ([13], [27])
        self.expected_variables = sorted(['ref_time', 'time_units', 'time', 'r', 't2m', 'u', 'v'])
        self.expected_levels = [200, 300, 500, 700, 800, 950, 1000]


    def test_open_hdf5(self):
        self.assertTrue(not self.hd is {})

    def test_variables_test(self):
        self.assertEqual(sorted(list(self.hd.keys())), self.expected_variables,
                         'incorrect number of variables')

    def test_varible_dimension(self):
        self.assertEqual(self.hd.get('u').isobaricInhPa.value.shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_levels(self):
        self.assertEqual(self.hd.get('u').isobaricInhPa.level, self.expected_levels,
                         'levels of u variable incorrect')

    def test_level_type(self):
        self.assertEqual(self.hd.get('u').level_type, self.expected_level_type,
                         'level type of u variable incorrect')

    def test_varible_units(self):
        self.assertEqual(self.hd.get('u').parameter_units, self.expected_units,
                         'units of u variable incorrect')

    def test_cut_time(self):
        self.assertListEqual(list(self.hd.get('time')), self.expected_times,
                             'incorrect time cut')

    def test_cut_space(self):
        self.assertEqual(near_yx2({'latitude': self.hd.get('u').latitude,
                                  'longitude': self.hd.get('u').longitude},
                                 lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')

    def test_write_data(self):
        np.testing.assert_almost_equal(self.hd.get('u').isobaricInhPa.value,
                                       self.new_hd.get('u').isobaricInhPa.value,
                                       decimal=3)

    def test_write_coord(self):

        np.testing.assert_almost_equal(self.hd.get('u').latitude,
                                       self.new_hd.get('u').latitude,
                                       decimal=3)
        np.testing.assert_almost_equal(self.hd.get('r').longitude,
                                       self.new_hd.get('r').longitude,
                                       decimal=3)

    def test_write_variables(self):
        self.assertEqual(sorted(self.new_hd.keys()), self.expected_variables,
                         'incorrect number of variables')


if __name__ == '__main__':
    unittest.main()

