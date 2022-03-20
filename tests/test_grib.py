from gdio.commons import near_yx
from gdio.grib import grib
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


class TestGribFiles(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        gr = grib(verbose=False)

        path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), ".")),
                            'data/era5_20191226-27_lev.grib')

        self.gbr = gr.gb_load(path,
                              cut_domain=(-30, 300, 10, 320),
                              cut_time=(12, 24),
                              rename_vars={'t': 't2m'})


    def setUp(self):

        self.expected_dim = (1, 2, 7, 160, 80)
        self.expected_variables = ['ref_time', 'time_units', 'time', 'r', 't2m', 'u', 'v']
        self.expected_level_type = ['isobaricInhPa']
        self.expected_units = 'm s**-1'
        self.expected_times = [12, 24]
        self.expected_coordinate = ([26], [53])
        self.expected_levels = [200, 300, 500, 700, 800, 950, 1000]

    def test_open_grib(self):

        self.assertTrue(not self.gbr is {})

    def test_variables_test(self):
        self.assertEqual(list(self.gbr.keys()), self.expected_variables,
                         'incorrect number of variables')

    def test_variables_rename(self):
        self.assertFalse(list(self.gbr.keys()) in self.expected_variables,
                         'variable rename fail')

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
        self.assertEqual(near_yx({'latitude': self.gbr.get('u').latitude,
                                  'longitude': self.gbr.get('u').longitude},
                                 lats=-23.54, lons=-46.64), self.expected_coordinate,
                         'problem with the spatial dimension')


if __name__ == '__main__':
    # unittest.main()
    gr = grib(verbose=False)

    path = 'data/era5_20191226-27_lev.grib' #'/home/rodrigo/data/ecmwf/ens/U_isobaric_ens.grb' #

    ds = gr.gb_load(path)

    # print(ds.time_units)
    # print(ds.keys(),gbr['10u'].surface.value.shape)
    # print(ds.keys(),gbr['t'].isobaricInhPa.value.shape)


    # import numpy as np
    # from datetime import datetime
    #
    # gr = grib(verbose=False)
    #
    # ds = {'ref_time': datetime(2019, 12, 27, 0, 0),
    #       'time_units': 'hours',
    #       'time': np.array(
    #           [0, 12]
    #           # [ datetime(2019, 12, 27, 0, 0), datetime(2019, 12, 27, 12, 0)]
    #       ),
    #       'u': {'isobaricInhPa': {'value': np.random.random((3, 2, 7, 80, 40)),
    #                               'level': [200, 300, 500, 700, 800, 950, 1000],
    #                               'members': [0, 1, 2]
    #                               },
    #             'param_id': None,
    #             'long_name': 'U component of wind',
    #             'level_type': ['isobaricInhPa'],
    #             'parameter_units': 'm s**-1',
    #             'longitude': np.array([300., 300.5, 301., 301.5, 302., 302.5, 303., 303.5,
    #                                    304., 304.5, 305., 305.5, 306., 306.5, 307., 307.5,
    #                                    308., 308.5, 309., 309.5, 310., 310.5, 311., 311.5,
    #                                    312., 312.5, 313., 313.5, 314., 314.5, 315., 315.5,
    #                                    316., 316.5, 317., 317.5, 318., 318.5, 319., 319.5]),
    #             'latitude': np.array([-30., -29.5, -29., -28.5, -28., -27.5, -27., -26.5,
    #                                   -26., -25.5, -25., -24.5, -24., -23.5, -23., -22.5,
    #                                   -22., -21.5, -21., -20.5, -20., -19.5, -19., -18.5,
    #                                   -18., -17.5, -17., -16.5, -16., -15.5, -15., -14.5,
    #                                   -14., -13.5, -13., -12.5, -12., -11.5, -11., -10.5,
    #                                   -10., -9.5, -9., -8.5, -8., -7.5, -7., -6.5,
    #                                   -6., -5.5, -5., -4.5, -4., -3.5, -3., -2.5,
    #                                   -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5,
    #                                   2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5,
    #                                   6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5]),
    #             }
    #       }
    #
    ofile = 'test.grb2'
    gr.gb_write(ofile, ds)
