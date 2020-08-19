from mpmath import mp, mpf
from contextlib import contextmanager

from .util import next_power

class IEEEFormat:
    def __init__(self, name, bits, exp, prec):
        self.bits = bits
        self.exp = exp
        self.prec = prec
        self.bias = 2 ** (self.exp - 1) - 1

    def min_exp(self, denorm=True):
        """
        Return the smallest possible exponent.

        If `denorm` is True, include denormalized values.
        """
        min_exp = -(self.bias - 1)
        if denorm:
            min_exp -= self.prec - 1
        return min_exp

    def max_exp(self):
        """
        Return the largest possible exponent.
        """
        return self.bias

    def min_value(self, denorm=True):
        """
        Return the smallest possible value.

        If `denorm` is True, include denormalized values.
        """
        return mp.power(2, self.min_exp(denorm))

    def max_value(self):
        """
        Return the largest possible value.
        """
        return mp.power(2, self.max_exp())

    def ulp(self, x):
        """
        Return the unit of least precision for a value x.

        This is the floating point value which represents the delta between
        x and the closest representible number.
        """
        from .util import zero, one, two
        x = mpf(x)
        n = 0

        min_exp = abs(self.min_exp(denorm=True))
        with mp.workprec(self.prec):
            ulp = next_power(x)
            while (x + ulp) != x and ulp != zero and n < min_exp:
                ulp /= two
                n += 1

        return ulp * two

    def dist(self, f1, f2):
        """
        Return the distance between two values in ULPs.
        """
        with mp.workprec(self.prec):
            lo = mpf(min(f1, f2))
            hi = mpf(max(f1, f2))
            dist = 0
            while lo != hi:
                lo += self.ulp(lo)
                dist +=1

            return dist


ieee_float =     IEEEFormat('single',   32,  8, 24)
ieee_double =    IEEEFormat('double',   64,  11, 53)
intel_extended = IEEEFormat('extended', 80,  15, 63)
ieee_quad =      IEEEFormat('quad',     128, 15, 113)

ieee_formats = dict((f.prec, f) for f in (
    ieee_float, ieee_double, intel_extended, ieee_quad
))


def ieee_format(prec=None):
    """
    Return the IEEE format for the given precision.

    :raises KeyError: If no standard format exists with `prec` bits of
        precision.
    """
    prec = prec or mp.prec
    return ieee_formats[prec]


def ieee_eval(x):
    """
    Return the string representation of a floating point value in each prec.
    """
    x = mpf(x)
    sreps = list()
    for prec in ieee_formats:
        with mp.workprec(prec):
            sreps.append(str(x))
    return sreps


@contextmanager
def workfloat(bytes_or_bits):
    ok = False
    bytes_or_bits = int(bytes_or_bits)
    for fmt in ieee_formats.values():
        if bytes_or_bits in (fmt.bits, fmt.bits // 8):
            ok = True
            with mp.workprec(fmt.prec):
                yield
    if not ok:
        raise KeyError(bytes_or_bits)
