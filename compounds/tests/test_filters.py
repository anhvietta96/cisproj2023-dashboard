
import unittest

from compounds.filters import build_molecules_filter


class FilterFunctionTests(unittest.TestCase):

    def test_filter_range(self):
        attribute_name = "abc"
        low = -9
        high = 6
        actual_filter = build_molecules_filter(attribute_name, (high, low))

        self.assertEquals(actual_filter, {
            "abc__gte": low,
            "abc__lte": high
        })

    def test_filter_contains(self):
        attribute_name = "abc"
        contains = "CoV"
        actual_filter = build_molecules_filter(attribute_name, contains)

        self.assertEquals(actual_filter, {
            "abc__icontains": contains,
        })

    def test_filter_invalid_type(self):
        with self.assertRaises(TypeError):
            build_molecules_filter("abc", 3)
