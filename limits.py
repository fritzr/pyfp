from mpmath import mp, mpf, workprec

from operator import eq, ne

from .util import zero, one, two
from .ieee import ieee_format


def find_fp_fast(
    func, limit, xstart, xstop, xstep, compare, max_steps=None,
):
    """
    Find the first power-of-2 such that |func(x) - limit(x)| <= |eps(x)|.

    To find the absolute limit, use bisect_limit on [x, x*2].

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


def find_limit_fast(
    func, limit, xstart=None, xstop=None, nulp=0, prec=None,
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
        return find_fp(
            func,
            limit,
            xstart,
            xstop,
            xstep,
            compare,
            max_steps=abs(fmt.min_exp()),
        )


def find_limit_reverse_fast(
    func, limit, xstart=None, xstop=None, nulp=0, prec=None,
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
        steps, value = find_fp(
            func,
            limit,
            xstart,
            xstop,
            xstep,
            compare,
            max_steps=abs(fmt.min_exp()),
        )
        return (abs(fmt.min_exp()) - steps + 1), (value / two)


def bisect_limit(
    func, limit, lower_limit, upper_limit=None, nulp=0, prec=None,
):

    fmt = ieee_format(prec)

    if upper_limit is None:
        upper_limit = 2 * lower_limit

    if nulp:

        def compare(fx, lx):
            return abs(fx - lx) > fmt.ulp(fx)

    else:
        compare = eq

    with workprec(prec or mp.prec):
        x = lower_limit
        last_x = None
        nsteps = 0
        while x != last_x and lower_limit != upper_limit:
            if compare(func(x), limit(x)):
                lower_limit = x
            else:
                upper_limit = x
            last_x = x
            x = (upper_limit + lower_limit) / 2
            nsteps += 1

        return nsteps, lower_limit


def find_limit(
    func, limit, xstart=None, xstop=None, nulp=0, prec=None,
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
        n_min, x_n = find_fp(
            func,
            limit,
            xstart,
            xstop,
            xstep,
            compare,
            max_steps=abs(fmt.min_exp()),
        )

        steps, limit = bisect_limit(func, limit, x_n, x_n * 2, nulp, prec)

    return limit, n_min, steps


def find_limit_reverse(
    func, limit, xstart=None, xstop=None, nulp=0, prec=None,
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
        n_min, x_n = find_fp(
            func,
            limit,
            xstart,
            xstop,
            xstep,
            compare,
            max_steps=abs(fmt.min_exp()),
        )

        n_min = abs(fmt.min_exp()) - n_min + 1
        steps, limit = bisect_limit(func, limit, x_n / two, x_n, nulp, prec)

    return limit, n_min, steps
