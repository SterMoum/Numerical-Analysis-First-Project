import math
import sympy as sp
from scipy.misc import derivative
import numpy
import random


def f():
    f = lambda x: 94 * (math.cos(x) ** 3) - 24 * math.cos(x) + 177 * (math.sin(x) ** 2) - 108 * (
            math.sin(x) ** 4) - 72 * (math.cos(x) ** 3) * (
                          math.sin(x) ** 2) - 65

    return f


def get_prime1():
    f = lambda x: 216 * (math.sin(x) ** 3) * (math.cos(x) ** 2) - 432 * (math.sin(x) ** 3) * math.cos(
        x) - 144 * math.sin(x) * (math.cos(x) ** 4) - 282 * math.sin(x) * (
                          math.cos(x) ** 2) + 354 * math.sin(x) * math.cos(x) + 24 * math.sin(x)

    return f


def get_prime2():
    f = lambda x: (1728 - 1800 * math.cos(x)) * (math.sin(x) ** 4) + (2358 * math.cos(x) - 2004) * (
            math.sin(x) ** 2) - 402 * math.cos(x) + 354

    return f


def newton_method(a, b):
    f_prime1 = get_prime1()
    f_prime2 = get_prime2()
    fx = f()
    x = a
    while fx(x) * f_prime2(x) <= 0:
        x += 0.1
    N = 1
    h = (fx(x) * (2 * (f_prime1(x) ** 2) - fx(x) * f_prime2(x))) / (2 * (f_prime1(x) ** 3))
    while abs(h) >= (1 / 2) * (10 ** (-5)):
        N += 1
        x = x - h
        h = (fx(x) * (2 * (f_prime1(x) ** 2) - fx(x) * f_prime2(x))) / (2 * (f_prime1(x) ** 3))
    print("TIMES:", N)
    return x


def bisection_method(x0, x1):
    fx = f()
    fa = fx(x0)
    fb = fx(x1)
    last_root = 0
    N = 1
    while abs(x1 - x0) > (10 ** (-5)) / 2:
        c = random.uniform(x0, x1)
        fc = fx(c)
        last_root = c
        if fc == 0:
            break
        if fa * fc < 0:
            x1 = c
            fb = fc
        else:
            x0 = c
            fa = fc
        N += 1
    print(N)
    return last_root


def secant_method(x0, x1, x2):
    fx = f()
    x3 = 0
    N = 1
    last = 1000000000
    root = x0
    q = fx(x0) / fx(x1)
    r = fx(x2) / fx(x1)
    s = fx(x2) / fx(x0)
    h = (r * (r - q) * (x2 - x1) + (1 - r) * s * (x2 - x0)) / ((q - 1) * (r - 1) * (s - 1))
    while abs(x3 - last) >= 0.5 * (10 ** (-5)):
        x3 = x2 - h
        last = root
        root = x3
        if fx(last) != 0:
            q = fx(root) / fx(x1)
            r = fx(x2) / fx(x1)
            s = fx(x2) / fx(root)
            h = (r * (r - q) * (x2 - x1) + (1 - r) * s * (x2 - root)) / ((q - 1) * (r - 1) * (s - 1))
        else:
            break
        N += 1
    print("TIMES:", N)
    return last


print("BISECTION:", bisection_method(0, 1))
print("BISECTION:", bisection_method(1, 2))
print("BISECTION:", bisection_method(2, 3))
print("NEWTON-RAPHSON:", newton_method(0, 1))
print("NEWTON-RAPHSON:", newton_method(1, 2))
print("NEWTON-RAPHSON:", newton_method(2, 3))
print("SECANT", secant_method(0.3, 0.5, 0.8))
print("SECANT", secant_method(1, 1.05, 1.20))
print("SECANT", secant_method(2, 2.5, 3))
