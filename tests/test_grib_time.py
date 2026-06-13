from datetime import datetime
import importlib
import sys
import types
import unittest


class TestGribTime(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._sentinel = object()
        cls._module_names = ['eccodes', 'numpy', 'pyproj', 'gdio.commons', 'gdio.cgrib', 'gdio.grib']
        cls._saved_modules = {
            name: sys.modules.get(name, cls._sentinel)
            for name in cls._module_names
        }

        for name in cls._module_names:
            sys.modules.pop(name, None)

        fake_eccodes = types.ModuleType('eccodes')
        fake_eccodes.CODES_PRODUCT_GRIB = 1
        fake_eccodes.KeyValueNotFoundError = type('KeyValueNotFoundError', (Exception,), {})
        fake_eccodes.ArrayTooSmallError = type('ArrayTooSmallError', (Exception,), {})
        fake_numpy = types.ModuleType('numpy')
        fake_pyproj = types.ModuleType('pyproj')

        for name in [
            'codes_get_native_type',
            'codes_get_string',
            'codes_get_double',
            'codes_get_array',
            'codes_get',
            'codes_get_values',
            'codes_get_message',
            'codes_count_in_file',
            'codes_clone',
            'codes_set',
            'codes_set_values',
            'codes_grib_new_from_file',
            'codes_release',
            'codes_write',
            'codes_new_from_samples',
        ]:
            setattr(fake_eccodes, name, lambda *args, **kwargs: None)

        sys.modules['eccodes'] = fake_eccodes
        sys.modules['numpy'] = fake_numpy
        sys.modules['pyproj'] = fake_pyproj
        cls.grib_module = importlib.import_module('gdio.grib')

    @classmethod
    def tearDownClass(cls):
        for name in cls._module_names:
            sys.modules.pop(name, None)

            original = cls._saved_modules[name]
            if original is not cls._sentinel:
                sys.modules[name] = original

    def test_forecast_step_respects_grib_day_units(self):
        gb = self.grib_module.grib(verbose=False)

        self.assertEqual(
            gb._grib__forecast_step(datetime(2020, 1, 3), datetime(2020, 1, 1), 2),
            2,
        )


if __name__ == '__main__':
    unittest.main()
