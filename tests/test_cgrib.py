import importlib
import sys
import types
import unittest


class TestCgribWriter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._sentinel = object()
        cls._module_names = ['eccodes', 'numpy', 'pyproj', 'gdio.cgrib']
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
        cls.cgrib = importlib.import_module('gdio.cgrib')

    @classmethod
    def tearDownClass(cls):
        for name in cls._module_names:
            sys.modules.pop(name, None)

            original = cls._saved_modules[name]
            if original is not cls._sentinel:
                sys.modules[name] = original

    def test_writer_uses_message_grid_type_for_mercator_projection(self):
        writer = self.cgrib.fwrite.__new__(self.cgrib.fwrite)
        message = {
            'gridType': 'mercator',
            'value': FakeGrid([[0.0, 0.0], [0.0, 0.0]]),
            'longitude': FakeGrid([[10.0, 11.0], [10.0, 11.0]]),
            'latitude': FakeGrid([[0.0, 0.0], [1.0, 1.0]]),
        }

        self.cgrib.fwrite._fwrite__set_grib_proj_keys(writer, message)

        self.assertEqual(message['Ni'], 2)
        self.assertEqual(message['Nj'], 2)
        self.assertIn('LaD', message)


class FakeGrid:

    def __init__(self, values):
        self.values = values
        self.shape = (len(values), len(values[0]))

    def __getitem__(self, index):
        if isinstance(index, tuple):
            row, col = index
            return self.values[row][col]

        return self.values[index]


if __name__ == '__main__':
    unittest.main()
