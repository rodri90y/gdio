import os, sys
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

        self.expected_dim = (1, 7, 80, 40)
        self.expected_times = [12]
        self.expected_coordinate = ([13], [27])



    def test_open_netcdf(self):
        self.assertTrue(not self.nc is {})

    def test_netcdf_variables_test(self):

        self.assertEqual(list(self.nc.keys()),
                         ['time', 'longitude', 'latitude', 'level', 'r', 't', 'u', 'v', 'time_unit', 'ref_time'],
                         'incorrect number of variables')

    def test_netcdf_varible_dimension(self):

        self.assertEqual(self.nc.get('u').shape, self.expected_dim,
                         'dimension shape of u variable incorrect')

    def test_grib_cut_time(self):

        self.assertListEqual(list(self.nc.get('time')), self.expected_times,
                         'incorrect time cut')


    def test_grib_cut_space(self):

        self.assertEqual(near_yx({'latitude': self.nc.get('latitude'),
                                  'longitude': self.nc.get('longitude')},
                                lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')



if __name__ == '__main__':
    # unittest.main()
    pass



