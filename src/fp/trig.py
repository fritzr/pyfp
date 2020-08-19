from mpmath import mpf, pi, cos, sin, cospi, sinpi, cospi_sinpi

from .util import one
from .limits import find_limit
from .ieee import ieee_format


d180 = mpf("180")
pio180 = pi / d180


# As x: 1 -> 0, find where cos(x) ~= 1
def cos_small(prec=None):
    return find_limit(cos, one, prec=prec)


def sin_approx(x):
    return x


# As x: 1 -> 0, find where sin(x) ~= x
def sin_small(prec=None):
    return find_limit(sin, sin_approx, prec=prec)


def rad2deg_ratio(radfunc):
    def degfunc(x):
        return radfunc(x / d180)

    return degfunc


def cosd(x):
    return cospi(x / d180)


def sind(x):
    return sinpi(x / d180)


def tanpi(x):
    return sinpi(x) / cospi(x)


def tand(x):
    return tanpi(x / d180)


# As x -> 0, find where cosd(x) ~= 1
def cosd_small(prec=None):
    return find_limit(cosd, one, prec=prec)


def sind_approx(x):
    return x * pio180


# As x -> 0, find where sind(x) ~= x.
#
# Search in reverse, since this can be quite small.
def sind_small(prec=None, nulp=0):
    return find_limit(sind, sind_approx, prec=prec, nulp=nulp)


# As x ->, tand(x) ~= x.
def tand_small(prec=None, nulp=0):
    return find_limit(tand, sind_approx, prec=prec, nulp=nulp)
