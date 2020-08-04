import os, sys
import unittest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gdio.grib import grib
from gdio.commons import near_yx


class TestGribFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        gr = grib(verbose=False)

        path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), ".")),
                            'data/era5_20191226-27_lev.grib')


        self.gbr = gr.gb_load(path,
                              cut_domain=(-30, 300, 10, 320),
                              cut_time=(12, 24))


    def setUp(self):

        self.expected_dim = (2, 7, 160, 80)
        self.expected_variables = ['ref_time', 'time_unit', 'time', 'latitude', 'longitude', 't', 'u', 'v', 'r']
        self.expected_times = [12, 24]
        self.expected_coordinate = ([26], [53])


    def test_open_grib(self):

        self.assertTrue(not self.gbr is {})

    def test_grib_variables_test(self):

        self.assertEqual(list(self.gbr.keys()), self.expected_variables,
                         'incorrect number of variables')

    def test_grib_varible_dimension(self):

        self.assertEqual(self.gbr.get('u').shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_grib_cut_time(self):

        self.assertListEqual(list(self.gbr.get('time')), self.expected_times,
                         'incorrect time cut')

    def test_grib_cut_space(self):

        self.assertEqual(near_yx({'latitude': self.gbr.get('latitude'),
                                  'longitude': self.gbr.get('longitude')},
                                lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')



if __name__ == '__main__':
    unittest.main()