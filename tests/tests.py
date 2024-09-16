# test_core.py
import unittest
from vira.core import Vira

class TestVira(unittest.TestCase):
    def test_summary(self):
        vira = Vira('data/sample_data.csv')
        self.assertEqual(vira.summary(), "Data summary: count=3, sum=6.0")

if __name__ == "__main__":
    unittest.main()
