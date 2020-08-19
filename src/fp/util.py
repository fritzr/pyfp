from mpmath import mpf, mp

zero = mp.zero
one = mp.one
two = one + one


def axrange(arg0, *args):
    """axrange([start, ]stop[, step])

    Iterator version of mpmath.arange.
    """
    end = arg0
    start = mpf('0')
    step = mpf('0.1')

    if args:
        start = arg0
        end = args[0]

        if len(args) > 1:
            if len(args) > 2:
                raise TypeError("Too many arguments")
            step = args[1]

    cur = start
    while cur < end:
        yield cur
        cur += step


def next_power(x, n=2):
    """
    Return the value sign(x) * n^k such that n^k is the smallest value > |x|
    for integer k.
    """
    x = abs(x)
    return mp.power(n, mp.floor(mp.log(x, n)) + 1)


def ceil_power(x, n=2):
    """
    Return the value sign(x) * n^k such that n^k is the smallest value >= |x|
    for integer k.
    """
    x = abs(x)
    logn = mp.log(x, n)
    if logn == mp.floor(logn):
        return x
    return mp.power(n, mp.floor(logn) + 1)


def prev_power(x, n=2):
    """
    Return the value sign(x) * n^k such that n^k is the largest value < |x|
    for integer k.
    """
    x = abs(x)
    logn = mp.log(x, n)
    if logn == mp.floor(logn):
        logn -= 1

    return mp.power(n, mp.floor(logn))


def floor_power(x, n=2):
    """
    Return the value sign(x) * n^k such that n^k is the largest value <= |x|
    for integer k.
    """
    x = abs(x)
    return mp.power(n, mp.floor(mp.log(x, n)))
