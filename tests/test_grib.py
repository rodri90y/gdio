from gdio.commons import near_yx2
from gdio.grib import grib
import os
import sys
import numpy as np
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


class TestGribFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        root = os.path.abspath(os.path.join(os.path.dirname(__file__), "."))
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addClassCleanup(self.tmpdir.cleanup)
        tmp_grib = os.path.join(self.tmpdir.name, 'tmp.grib')

        # open new new grib
        gr = grib(verbose=False)
        self.gbr = gr.gb_load(os.path.join(root, 'data/era5_2membros_sintetico.grib'),
                              cut_domain=(-30, 300, 10, 320),
                              cut_time=(1, 2),
                              rename_vars={'t': '2t'})

        # filtered member
        self.filtered_gbr = gr.gb_load(os.path.join(root, 'data/era5_2membros_sintetico.grib'),
                                       cut_domain=(-30, 300, 10, 320),
                                       cut_time=(1, 2),
                                       rename_vars={'t': '2t'},
                                       filter_by={'perturbationNumber': [1]}
                                       )


        # write new grib
        gr.gb_write(tmp_grib, self.gbr,
                    least_significant_digit=3,
                    packingType='grid_jpeg')

        # open new new grib
        self.new_gbr = gr.gb_load(tmp_grib,
                               rename_vars={'t': '2t'}
                                  )





    def setUp(self):

        self.expected_dim = (2, 2, 7, 161, 81)
        self.expected_variables = sorted(['ref_time', 'time_units', 'time', 'r', '2t', 'u', 'v'])
        self.expected_level_type = ['isobaricInhPa']
        self.expected_units = 'm s**-1'
        self.expected_times = [12, 24]
        self.expected_coordinate = ([26], [53])
        self.expected_levels = [200, 300, 500, 700, 800, 950, 1000]
        self.filtered_dim = (1, 2, 7, 161, 81)
        self.filtered_members = [1]
        self.expected_members = [0, 1]

    def test_open_grib(self):

        self.assertTrue(self.gbr)

    def test_variables(self):
        self.assertEqual(sorted(self.gbr.keys()), self.expected_variables,
                         'incorrect number of variables')


    def test_varible_dimension(self):
        self.assertEqual(self.gbr.get('u').isobaricInhPa.value.shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_levels(self):
        self.assertEqual(self.gbr.get('u').isobaricInhPa.level, self.expected_levels,
                         'levels of u variable incorrect')

    def test_level_type(self):
        self.assertEqual(self.gbr.get('u').level_type, self.expected_level_type,
                         'level type of u variable incorrect')

    def test_varible_units(self):

        self.assertEqual(self.gbr.get('u').parameter_units, self.expected_units,
                         'units of u variable incorrect')

    def test_cut_time(self):
        self.assertListEqual(list(self.gbr.get('time')), self.expected_times,
                             'incorrect time cut')

    def test_cut_space(self):
        self.assertEqual(near_yx2({'latitude': self.gbr.get('u').latitude,
                                  'longitude': self.gbr.get('u').longitude},
                                 lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')


    def test_filtered_dimension(self):
        self.assertEqual(self.filtered_gbr.get('u').isobaricInhPa.value.shape, self.filtered_dim,
                         'dimension shape of u variable incorrect')

    def test_filtered_members(self):
        self.assertEqual(sorted(self.filtered_gbr.get('u').isobaricInhPa.members), self.filtered_members,
                         'incorrect members')


    def test_write_data(self):
        # print(self.gbr.get('u').isobaricInhPa.value.shape)
        np.testing.assert_almost_equal(self.gbr.get('u').isobaricInhPa.value,
                                       self.new_gbr.get('u').isobaricInhPa.value,
                                       decimal=3)

    def test_write_coord(self):

        np.testing.assert_almost_equal(self.gbr.get('u').latitude,
                                       self.new_gbr.get('u').latitude,
                                       decimal=3)
        np.testing.assert_almost_equal(self.gbr.get('r').longitude,
                                       self.new_gbr.get('r').longitude,
                                       decimal=3)

    def test_write_variables(self):
        self.assertEqual(sorted(self.new_gbr.keys()), self.expected_variables,
                         'incorrect number of variables')

    def test_write_members(self):
        self.assertEqual(sorted(self.new_gbr.get('u').isobaricInhPa.members), self.expected_members,
                         'incorrect number of members')



if __name__ == '__main__':
    unittest.main()
