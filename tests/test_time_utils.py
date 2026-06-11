from datetime import datetime
import importlib
import importlib.util
import sys
import types
import unittest


class TestTimeUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._sentinel = object()
        cls._saved_modules = {
            name: sys.modules.get(name, cls._sentinel)
            for name in ['numpy', 'gdio.commons']
        }

        sys.modules.pop('gdio.commons', None)

        if importlib.util.find_spec('numpy') is None:
            sys.modules['numpy'] = types.ModuleType('numpy')

        cls.commons = importlib.import_module('gdio.commons')

    @classmethod
    def tearDownClass(cls):
        for name in ['gdio.commons', 'numpy']:
            sys.modules.pop(name, None)

            original = cls._saved_modules[name]
            if original is not cls._sentinel:
                sys.modules[name] = original

    def test_time_unit_to_hours(self):
        self.assertEqual(self.commons.time_unit_to_hours('seconds'), 1 / 3600)
        self.assertEqual(self.commons.time_unit_to_hours('minutes'), 1 / 60)
        self.assertEqual(self.commons.time_unit_to_hours('hours'), 1)
        self.assertEqual(self.commons.time_unit_to_hours('days'), 24)
        self.assertEqual(self.commons.time_unit_to_hours('months'), 24 * 30)
        self.assertEqual(self.commons.time_unit_to_hours('years'), 24 * 365)

    def test_parse_time_units(self):
        unit, ref_time = self.commons.parse_time_units('hours since 2019-12-27 12:30')

        self.assertEqual(unit, 'hours')
        self.assertEqual(ref_time, datetime(2019, 12, 27, 12, 30))

    def test_parse_invalid_time_units(self):
        self.assertEqual(self.commons.parse_time_units('not a time unit'), (None, None))


if __name__ == '__main__':
    unittest.main()
