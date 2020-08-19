from mpmath import mp, mpf, workprec

from operator import eq, ne

from .util import zero, one, two
from .ieee import ieee_format


def find_fp(
    func,
    limit,
    xstart,
    xstop,
    xstep,
    compare,
    max_steps=None,
):
    """
    Find the first value such that |func(x) - limit(x)| <= |eps(x)|.

    The initial value is `init`, and each subsequent value is computed by
    step(x).

    We use "<=" so that when `eps` is zero, the expression is equivalent to
    ``func(x) == limit(x)``.

    If `limit` or `eps` are values rather than callables, the corresponding
    expressions ``limit(x)`` or ``eps(x)`` always evaluate to the original value
    of the corresponding parameter.

    If `eps` is not given, eps(x) === ulp(x) (unit of least precision).

    The `prec` temporarily sets overrides ``mp.prec``.

    :returns: (x, n)

        Where x === 2^-n and |func(x) - limit(x)| <= |eps(x)|.
    """

    if not callable(limit):
        limit_value = mpf(limit)

        def limit(x):
            return limit_value

    nsteps = 0
    x = xstart

    while (
        (xstop is None or x != xstop)
        and (max_steps is None or nsteps != max_steps)
        and not compare(func(x), limit(x))
    ):
        x = xstep(x)
        nsteps += 1

    return nsteps, x

def find_limit(
    func,
    limit,
    xstart=None,
    xstop=None,
    nulp=0,
    prec=None,
):

    fmt = ieee_format(prec)

    if xstart is None:
        xstart = mp.one

    if xstop is None:
        xstop = fmt.min_value()

    def xstep(x):
        return x / two


    if nulp:
        def compare(fx, lx):
            return abs(fx - lx) <= fmt.ulp(fx)

    else:
        compare = eq

    with workprec(prec or mp.prec):
        return find_fp(func, limit, xstart, xstop, xstep, eq,
                max_steps=abs(fmt.min_exp()))


def find_limit_reverse(
    func,
    limit,
    xstart=None,
    xstop=None,
    nulp=0,
    prec=None,
    ):

    fmt = ieee_format(prec)

    if xstart is None:
        xstart = fmt.min_value()

    if xstop is None:
        xstop = one

    def xstep(x):
        return x * two

    if nulp:
        def compare(fx, lx):
            return abs(fx - lx) > fmt.ulp(fx)

    else:
        compare = ne

    with workprec(prec or mp.prec):
        steps, value = find_fp(func, limit, xstart, xstop, xstep, ne,
                max_steps=abs(fmt.min_exp()))
        return (abs(fmt.min_exp()) - steps + 1), (value / two)
