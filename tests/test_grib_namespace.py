import unittest

from gdio.definitions.grib_namespace import DATA_TIME_KEYS


class TestGribNamespace(unittest.TestCase):

    def test_validity_and_unit_time_keys_are_separate(self):
        self.assertIn('validityTime', DATA_TIME_KEYS)
        self.assertIn('unitOfTimeRange', DATA_TIME_KEYS)
        self.assertNotIn('validityTimeunitOfTimeRange', DATA_TIME_KEYS)


if __name__ == '__main__':
    unittest.main()
