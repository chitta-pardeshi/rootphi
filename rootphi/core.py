"""
rootphi.py - Exact arithmetic with numbers in the field Q(φ, √φ),
where φ = (1 + √5) / 2 is the golden ratio.

Supports rational combinations of 1, φ, √φ, and φ√φ:
    (a + b√φ + cφ + dφ√φ) / e
All operations are exact and symbolic.
"""
class rootphi:

    def __init__(self, arr):
        if isinstance(arr, int):
            arr = [arr, 0, 0, 0, 1]
        if not isinstance(arr, list):
            raise ValueError("expecting an integer or a list of 5 integers.")
        if len(arr) != 5:
            raise ValueError("expecting an integer or a list of 5 integers.")
        if not all(isinstance(x, int) for x in arr):
            raise ValueError("expecting an integer or a list of 5 integers.")
        if arr[4] == 0:
            raise ValueError("expecting an integer or a list of 5 integers.")
        g = arr[4]
        for n in arr:
            while n != 0:
                g, n = n, g % n
        g = ( (arr[4] > 0) - (arr[4] < 0) ) * abs(g)
        self.a = arr[0] // g
        self.b = arr[1] // g
        self.c = arr[2] // g
        self.d = arr[3] // g
        self.e = arr[4] // g

    @classmethod
    def zero(cls): 
        return cls([0, 0, 0, 0, 1])
    @classmethod
    def one(cls): 
        return cls([1, 0, 0, 0, 1])
    @classmethod
    def minusone(cls): 
        return cls([-1, 0, 0, 0, 1])
    @classmethod
    def two(cls): 
        return cls([2, 0, 0, 0, 1])
    @classmethod
    def ten(cls): 
        return cls([10, 0, 0, 0, 1])
    @classmethod
    def phi(cls): 
        return cls([0, 0, 1, 0, 1])
    @classmethod
    def sqrtphi(cls): 
        return cls([0, 1, 0, 0, 1])

    @classmethod
    def from_json(cls, data):
        return cls([data["a"], data["b"], data["c"], data["d"], data["e"]])

    def __str__(self):
        return f"({self.a} + {self.b}√φ + {self.c}φ + {self.d}φ√φ) / {self.e}"

    def __repr__(self):
        return f"rootphi([{self.a}, {self.b}, {self.c}, {self.d}, {self.e}])"

    def __format__(self, format_spec):
        return format(self.to_float(), format_spec)

    def __hash__(self):
        return hash((self.a, self.b, self.c, self.d, self.e))

    def __int__(self):
        return int(self.to_float())

    def __float__(self):
        return self.to_float()


    def to_list(self):
        return [self.a, self.b, self.c, self.d, self.e]
    
    def to_json(self):
        return {
            "a": self.a, "b": self.b, "c": self.c,
            "d": self.d, "e": self.e
        }

    def __json__(self):
        return self.to_json()


    def as_tuple(self):
        return (self.a, self.b, self.c, self.d, self.e)

    def to_string(self, grouped: bool = True) -> str:
        def format_term(coeff, symbol=""):
            if coeff == 0:
                return ""
            sign = "-" if coeff < 0 else "+"
            abs_val = abs(coeff)
            if abs_val == 1 and symbol:
                term = symbol
            else:
                term = f"{abs_val}{symbol}"
            return f" {sign} {term}"

        if grouped:
            inner1 = "".join([
                format_term(self.a),
                format_term(self.c, "φ")
            ]).strip()
            inner2 = "".join([
                format_term(self.b),
                format_term(self.d, "φ")
            ]).strip()
            inner1 = inner1.lstrip("+ ").rstrip() or "0"
            if inner2:
                grouped_expr = f"({inner1} + √φ({inner2}))"
            else:
                grouped_expr = f"({inner1})"
        else:
            parts = [
                format_term(self.a),
                format_term(self.b, "√φ"),
                format_term(self.c, "φ"),
                format_term(self.d, "φ√φ")
            ]
            expr = "".join(parts).strip()
            grouped_expr = f"({expr.lstrip('+ ').rstrip() or '0'})"

        if self.e == 1:
            return grouped_expr
        else:
            return f"{grouped_expr} / {self.e}"


    def to_float(self):
        return float(self.as_decimal_string(precision=40))


    def to_decimal(self, precision=40):
        from decimal import Decimal
        return Decimal(self.as_decimal_string(precision=precision))

    def to_sympy(self):
        from sympy import Rational, sqrt, Symbol
        phi = Rational(1, 2) * (1 + sqrt(5))
        return Rational(self.a, self.e) + Rational(self.b, self.e) * sqrt(phi) + \
            Rational(self.c, self.e) * phi + Rational(self.d, self.e) * phi * sqrt(phi)

    # Properties & Predicates

    def simplify(self):
        return rootphi([self.a, self.b, self.c, self.d, self.e])

    def is_integer(self):
        return self.b == 0 and self.c == 0 and self.d == 0 and self.e == 1

    def is_rational(self):
        return self.b == 0 and self.c == 0 and self.d == 0

    def is_zero(self):
        return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

    def is_positive(self):
        return self.sign() > 0

    def is_negative(self):
        return self.sign() < 0


    def denominator(self):
        return self.e

    def coefficient(self, e: int):
        if e == 0:
            return self.a
        elif e == 1:
            return self.b
        elif e == 2:
            return self.c
        elif e == 3:
            return self.d
        else:
            return NotImplemented


    def digits_base_repr(self, base=10, precision=20):
        """
        Returns [sign, integer_part_digits[], fractional_part_digits[]]
        
        sign: -1, 0, or 1
        base: int or rootphi
        precision: number of digits after decimal
        """
        if isinstance(base, int):
            base = rootphi(base)
        if not isinstance(base, rootphi):
            raise ValueError("base must be int or rootphi")

        if not isinstance(precision, int) or precision < 1:
            raise ValueError("precision must be a positive integer")

        sign = self.sign()
        if sign == 0:
            return [0, [0], [0] * precision]

        # Work with absolute value
        n = abs(self)

        # Compute integer part digits
        int_digits = []
        temp = rootphi(1)
        factor = rootphi(1)
        count = 0

        # Find largest power of base ≤ n
        powers = []
        current = rootphi(1)
        while current <= n:
            powers.append(current)
            current = current * base

        powers.reverse()  # Start from highest base power
        remaining = n
        for power in powers:
            digit = 0
            while (remaining - power).sign() >= 0:
                remaining = remaining - power
                digit += 1
            int_digits.append(digit)

        if not int_digits:
            int_digits = [0]

        # Compute fractional digits
        frac_digits = []
        frac = remaining
        for _ in range(precision):
            frac = frac * base
            digit = 0
            while (frac - rootphi(digit + 1)).sign() >= 0:
                digit += 1
            frac_digits.append(digit)
            frac = frac - rootphi(digit)

        return [sign, int_digits, frac_digits]

    def as_decimal_string(self, precision=20):
        sign, int_part, frac_part = self.digits_base_repr(10, precision)
        sign_str = "-" if sign < 0 else ""
        int_str = ''.join(str(d) for d in int_part)
        frac_str = ''.join(str(d) for d in frac_part)
        return f"{sign_str}{int_str}.{frac_str}"

    def sign(self):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        sign_e = ( (e > 0) - (e < 0) )
        sign_c = ( (c > 0) - (c < 0) )
        sign_d = ( (d > 0) - (d < 0) )
        val_2ac = 2*a + c
        sign_2ac = ( (val_2ac > 0) - (val_2ac < 0))
        val_2bd = 2*b + d 
        sign_2bd = ( (val_2bd > 0) - (val_2bd < 0))
        val_5c2ac = sign_e * ( 5 * sign_c * c**2 + sign_2ac * val_2ac**2)
        sign_5c2ac = ( (val_5c2ac > 0) - (val_5c2ac < 0) )
        val_5d2bd = sign_e * (5 * sign_d * d **2 + sign_2bd * val_2bd ** 2)
        sign_5d2bd = ( (val_5d2bd > 0) - (val_5d2bd < 0) )

        terma = sign_5c2ac*e**2*a**2 + (sign_5c2ac*e**2*c**2 + (2*sign_5d2bd*d*b + sign_5d2bd*d**2)*e**2)
        termb = 2*sign_5c2ac*e**2*c*a + (sign_5c2ac*e**2*c**2 + (sign_5d2bd*b**2 + 2*sign_5d2bd*d*b + 2*sign_5d2bd*d**2)*e**2)
        sign_termb = ( (termb > 0) - (termb < 0) )

        val_2termab = 2*terma + termb
        sign_2termab = ( (val_2termab > 0) - (val_2termab < 0) )

        term = 5 * sign_termb * termb**2 + sign_2termab*val_2termab**2
        sign_term = ( (term > 0) - (term < 0))
        return sign_term

    def __abs__(self):
        return self * self.sign()

    def __pos__(self):
        return rootphi([self.a, self.b, self.c, self.d, self.e])

    def __neg__(self):
        return rootphi([-self.a, -self.b, -self.c, -self.d, self.e])
    
    def __add__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = self.a * other.e + other.a * self.e
        b = self.b * other.e + other.b * self.e
        c = self.c * other.e + other.c * self.e
        d = self.d * other.e + other.d * self.e
        e = self.e * other.e
        return rootphi([a, b, c, d, e])

    def __sub__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = self.a * other.e - other.a * self.e
        b = self.b * other.e - other.b * self.e
        c = self.c * other.e - other.c * self.e
        d = self.d * other.e - other.d * self.e
        e = self.e * other.e
        return rootphi([a, b, c, d, e])

    def __mul__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        a = other.a*self.a + other.c*self.c + other.d*self.b  + other.b*self.d + other.d*self.d 
        b = other.b*self.a + other.d*self.c + other.a*self.b + other.c*self.d 
        c = other.c*self.a + other.a*self.c + other.c*self.c + other.b*self.b + other.d*self.b + other.b*self.d + 2*other.d*self.d 
        d = other.d*self.a + other.b*self.c + other.d*self.c + other.c*self.b + other.a*self.d + other.c*self.d 
        e = self.e*other.e
        return rootphi([a, b, c, d, e])

    def __truediv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        f = other.a
        g = other.b
        h = other.c
        i = other.d 
        j = other.e

        na = j*(a*(f**3 + 2*f**2*h - f*(g**2 + 4*g*i + 3*i**2) + g**2*h + 2*g*h*i - h**3 + 2*h*i**2) - b*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - c*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - d*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3))
        nb = j*(-a*(f**2*g + 2*f*h*(g - i) - g**3 - 3*g**2*i + g*(2*h**2 - i**2) - h**2*i + 2*i**3) + b*(f**3 + 2*f**2*h - f*(g**2 + 4*g*i + 3*i**2) + g**2*h + 2*g*h*i - h**3 + 2*h*i**2) - c*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - d*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2))
        nc = j*(-a*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - b*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3) + c*(f**3 + f**2*h - f*(2*g*i + h**2 + i**2) + g**2*h + h*i**2) - d*(f**2*(g + 2*i) - 2*f*h*(g + i) + g**3 + g**2*i - g*i**2 + h**2*i))
        nd = j*(-a*(f**2*i - 2*f*g*h + g**3 + 2*g**2*i - g*h**2 + h**2*i - i**3) - b*(f**2*h - f*(g**2 + 2*g*i - h**2 + 2*i**2) + 2*g*h*i - h**3 + h*i**2) - c*(f**2*(g + i) - 2*f*h*i - g**2*i + g*(h**2 - i**2) + i**3) + d*(f**3 + f**2*h - f*(2*g*i + h**2 + i**2) + g**2*h + h*i**2))
        ne = -e*(-f**4 - 2*f**3*h + f**2*(g**2 + 6*g*i + h**2 + 4*i**2) - 2*f*h*(2*g**2 + 2*g*i - h**2 + 3*i**2) + g**4 + 2*g**3*i - g**2*(h**2 + i**2) + 2*g*i*(2*h**2 - i**2) - h**4 + h**2*i**2 + i**4)
        return rootphi([na, nb, nc, nd, ne])

    def __floordiv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        
        sign_other = other.sign()
        sign_self = self.sign()

        selfval = abs(self)
        othval = abs(other)
        hi = 2
        while  (selfval-othval*hi).sign() > 0:
            hi *= 2
        lo = 0
        while lo < hi - 1:
            m = (lo+hi)//2
            signm = (selfval-othval*m).sign()
            if signm < 0:
                hi = m
            elif signm > 0:
                lo = m
            else:
                hi = lo = m
        return sign_self * sign_other * lo

    def __mod__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        
        return self - other * (self // other)

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError("Exponent must be an integer.")
        
        if exponent == 0:
            return rootphi.one()

        if exponent == 1:
            return self

        if self.sign() == 0:
            if exponent < 0:
                raise ZeroDivisionError("0 cannot be raised to a negative power.")
            return rootphi.zero()
        
        result = rootphi.one()
        base = self
        abs_exponent = abs(exponent)
        
        while abs_exponent > 0:
            if abs_exponent % 2 == 1:
                result *= base
            base *= base
            abs_exponent //= 2
        
        if exponent < 0:
            result = rootphi.one() / result
        
        return result


    def __iadd__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self + other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __isub__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self - other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __imul__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        
        result = self * other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __itruediv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")

        result = self / other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __ifloordiv__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ZeroDivisionError("Division by zero.")
        
        result = self // other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e

    def __imod__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        if other.sign() == 0:
            raise ValueError("Division by zero.")
        
        result = self % other
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e


    def __ipow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError("Exponent must be an integer.")
        
        result = self ** exponent
        self.a = result.a
        self.b = result.b
        self.c = result.c
        self.d = result.d
        self.e = result.e
        return self

    def __eq__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() == 0

    def __ne__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() != 0

    def __lt__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() < 0

    def __le__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() <= 0

    def __gt__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() > 0

    def __ge__(self, other):
        if isinstance(other, int):
            other = rootphi(other)
        if not isinstance(other, rootphi):
            return NotImplemented
        return (self - other).sign() >= 0

