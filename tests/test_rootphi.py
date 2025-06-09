import unittest
from rootphi import rootphi

class RootPhiTest(unittest.TestCase):
    def test_zero(self):
        self.assertTrue(rootphi.zero().is_zero())

    def test_add_sub(self):
        a = rootphi([1, 2, 3, 4, 5])
        b = rootphi([5, 4, 3, 2, 1])
        self.assertEqual((a + b) - b, a)

    def test_mul_div(self):
        x = rootphi([3, 1, 2, 0, 1])
        y = rootphi([2, 0, 0, 0, 1])
        z = x * y
        self.assertEqual(z / y, x)

    def test_sign(self):
        self.assertEqual(rootphi([1, 0, 0, 0, 1]).sign(), 1)
        self.assertEqual(rootphi([-1, 0, 0, 0, 1]).sign(), -1)

    def test_format_str(self):
        x = rootphi([1, -1, 0, 0, 1])
        self.assertIn("√φ", str(x))
        self.assertTrue(isinstance(str(x), str))

    def test_sympy_compare(self):
        import sympy
        x = rootphi([1, 0, 1, 0, 1])
        sym = x.to_sympy()
        dec = float(sym.evalf(30))
        self.assertAlmostEqual(x.to_float(), dec, places=12)

    def test_decimal_accuracy(self):
        x = rootphi([198753070236573275664441, 0, 0, 0, 156250000000000000000000])
        expected = 198753070236573275664441 / 156250000000000000000000
        self.assertAlmostEqual(x.to_float(), expected, places=12)

if __name__ == "__main__":
    unittest.main()
