pyfp
====

Simple package for floating-point precision tests, based on the
[mpmath](https://github.com/fredrik-johansson/mpmath) package.

Installing
----------

I find the easiest way to install is:

```bash
$ git clone https://github.com/fritzr/pyfp
$ pip install ./pyfp
```

fp.ieee
-------

Utilities for working in terms of IEEE-754 floating point, where one thinks
about bits before precision.

```python
>>> from mpmath import mp, mpf
>>> from fp.ieee import ieee_float, ieee_double, intel_extended, workfloat
>>> x = mpf("1.234567890987654321")


# Set the precision:
>>> mp.prec = ieee_float.prec
>>> x
mpf('1.23456789')

# Get the current IEEE format (only works for standardized precision values):
>>> current_format = ieee_format(mp.prec)
>>> current_format.bits
32

# Set the precision temporarily using mp.workrprec:
>>> with mp.workprec(ieee_float.prec):  # IEEE-754 singe precision (32 bits)
...     print(x)
...
1.234567

# Set the precision in terms of bits or bytes using ieee.workfloat:
>>> with workfloat(64):  # IEEE-754 double precision (64 bits)
...     print(x)
...
1.23456789098765

>>> with workfloat(intel_extended.bits):  # Intel extended precision (80 bits)
...     print(x)
...
1.23456789098765429

# Unit in the last place (ULP):
>>> ulp_at_one = current_format.ulp(1)
>>> ulp_at_one
mpf('1.1920929e-7')
>>> print(1 == 1 + ulp_at_one)
False
>>> print(1 == 1 + (ulp_at_one / 2))
True

# Distance (in # ULPs):
>>> ieee_float.dist('1.00001001', '1.00001013')
1
>>> ieee_float.dist('1.00001001', '1.000011') < 10  # approximately equal
True

```

fp.trig and fp.limits
---------------------

Extra multi-precision trigonometric functions and limit checking.

```python
from mpmath import mp, mpf
from fp.trig import *

>>> sin(pi / 2) == sind(90)
True
>>> cos(pi / 4) == cosd(45)
True
```

To optimize an implementation of sin(x) we may want to take advantage of the
property that sin(*x*) ~= *x* for small *x*. Then we might wonder how small
must *x* be for us to be able to use the identity without loss of precision?

To formalize this, we define *R_b*{*x*} as the closest number to x which is
representible in the *b*-bit floating point standard. We would like to find
*x_0* such that for all *x* <= *x_0*, |*R_b*{sin(*x_0*)} - *R_b*{*x_0*}| < *e*
for any desired *e* >= 0.

First let's play with a sample number:

```python
>>> x = mpf('0.01')
>>> sin(x)
mpf('0.00999983307')
>>> abs(sin(x) - x)
mpf('1.66706741e-7')
>>> ieee_float.dist(x, sin(x))
179
```

As you can see, sin(*x*) is relatively far from *x* in the sense that there are
179 distinct values we could represent between the two. How can we find the
maximal value *x_0* for which *x* <= *x_0* implies *R_32*{sin(*x*)} ==
*R_32*{*x*}?  Well, `fp.limits.find_limit` computes this for us:

```python
>>> def sin_approx(x):
...     return x
>>> x_0, n_min, b_steps = find_limit(sin, sin_approx)
>>> x_0
mpf('0.000443632889')
>>> sin(x_0) == x_0
True
>>> x_1 = x_0 + ieee_float.ulp(x_0)
>>> sin(x_1) == x_1
False
```

We can see that *x_0* is the largest value that satisfies the condition,
since it no longer holds for the very next representible value *x_1*.

This function first finds *x_l* = 2^*-n*: the largest power of two for which
sin(*x*) == *x*. Then it bisects the region [*x_l*, 2 \* *x_l*] to find the
upper limit *x_0*. We can see this with `find_limit_fast`:

```python
>>> n_min, x_l = find_limit_fast(sin, sin_approx)
>>> n_min
12
>>> x_l == mpf(2) ** -n_min
True
>>> x_l <= x_0 <= 2 * x_l
True
```

In this result, *x_l* certainly satisfies the limit condition, but is not the
maximum such value. To find the upper limit of *x_0*, one can bisect the
region between *x_l* and *x_l* * 2 on the limit condition. This is done by
feeding the result of `find_limit_fast` into `bisect_limit`:

```python
>>> steps, limit = bisect_limit(sin, sin_approx, x_l, x_l * 2)
>>> steps
25
>>> limit == x_0
True
```

There are a few simple functions which compute these limits:

```python
>>> sin_small()[0]  # x_0 : x <= x_0 implies sin(x) == x
mpf('000443632889')
>>> sind_small()[0]  # x_0 : x <= x_0 implies sind(x) == x * pi / 180
mpf('0.0391264558')
>>> tand_small()[0]  # x_0 : x <= x_0 implies tand(x) == x * pi / 180
mpf('0.0139759379')
```

fp.util
-------

Miscellaneous utilities.

```python
from mpmath import mp, mpf
from fp import util

# Generator like mpmath.arange but doesn't hold all the values in memory.
>>> mp.prec = 53
>>> util.axrange(0, mpf('1e130'), mpf('0.01'))
<generator object axrange at ...>

# Powers.

>>> next_power(mpf('128'), 2)
mpf('256.0')

>>> prev_power(mpf('128'), 2)
mpf('64.0')

>>> floor_power(mpf('127'), 2)
mpf('64.0')
>>> floor_power(mpf('128'), 2)
mpf('128.0')
>>> floor_power(mpf('129'), 2)
mpf('128.0')

>>> ceil_power(mpf('127'), 2)
mpf('128.0')
>>> ceil_power(mpf('128'), 2)
mpf('128.0')
>>> ceil_power(mpf('129'), 2)
mpf('256.0')

```
