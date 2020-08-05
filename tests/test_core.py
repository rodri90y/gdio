import os, sys
import unittest
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gdio.core import gdio
from gdio.commons import near_yx
from datetime import datetime

class TestNcFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        self.ds = gdio(verbose=False)

        self.ds.mload(['data/era5_20191226-27_lev.grib','data/era5_20191227_lev.nc'],
                      merge_files=True,
                      uniformize_grid=True,
                      cut_domain=(-30, 300, 10, 320),
                      cut_time=(12, 36),
                      inplace=True
                      )

    def setUp(self):

        self.expected_dim = (4, 7, 160, 80)
        self.expected_times = [datetime(2019,12,26,12,0), datetime(2019,12,27,0,0),
                               datetime(2019,12,27,12,0), datetime(2019,12,27,12,0)]
        self.expected_coordinate = ([26], [53])
        self.expected_variables = ['ref_time', 'time_units', 'time', 'latitude', 'longitude', 't', 'u', 'v', 'r', 'level']
        self.expected_corrcoef = 0.94
        self.expected_sel = (1, 4, 6, 18)


    def test_open_multiples_files(self):
        self.assertTrue(not self.ds.dataset is [])

    def test_grid_variables_test(self):

        self.assertEqual(list(self.ds.dataset[0].keys()), self.expected_variables,
                         'incorrect number of variables')

    def test_grid_varible_dimension(self):

        self.assertEqual(self.ds.dataset[0].get('u').shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_grid_cut_time(self):

        self.assertListEqual(list(self.ds.dataset[0].get('time')), self.expected_times,
                         'incorrect time cut')


    def test_grid_interpolation(self):

        a = self.ds.dataset[0].get('u')[2, 0].flatten()
        b = self.ds.dataset[0].get('u')[3, 0].flatten()

        corrcoef = np.corrcoef(np.nan_to_num(a), np.nan_to_num(b))[0, 1]

        self.assertTrue(corrcoef > self.expected_corrcoef,
                         'incorrect interpolation (corrcoef={0:0.4f})'.format(corrcoef))

    def test_grid_cut_space(self):

        self.assertEqual(near_yx({'latitude': self.ds.dataset[0].get('latitude'),
                                  'longitude': self.ds.dataset[0].get('longitude')},
                                lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')

    def test_grid_select_coords(self):

        sub_set = self.ds.sel(dates=[datetime(2019,12,26,12,0)],
                              latitude=[-23.54,-22], longitude=[-46.64,-42.2], level=[2,6])

        self.assertEqual(sub_set[0]['u'].shape, self.expected_sel,
                         'dimension shape of u variable incorrect')

if __name__ == '__main__':
    unittest.main()



