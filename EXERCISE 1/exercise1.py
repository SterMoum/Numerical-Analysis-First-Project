import math


def get_prime():
    # Defining the first Derivative
    f = lambda x: 3 * math.e ** (math.sin(x) ** 3) * math.cos(x) * math.sin(x) ** 2 + x ** 2 * (
            6 * (x ** 3) - 8 * x - 3)

    return f


def get_prime2():
    # Defining the second Derivative
    f = lambda x: 3 * (math.e ** (math.sin(x) ** 3) * math.sin(x) * (
            (math.sin(x) ** 2) * (3 * ((math.cos(x)) ** 2) * math.sin(x) - 1) + 2 * math.cos(x) ** 2) + 2 *
                       (x - 1) * x * (5 * x ** 2 + 5 * x + 1))

    return f


def bisectionMethod(F, e, x0, x1):
    # Defining Repeats
    N = (math.log(x1 - x0, math.e) - math.log((10 ** (-5)) / 2, math.e)) / math.log(2, math.e)
    N = round(N)
    print("TIMES:", N)

    # initializing values
    fa = F(x0)
    fb = F(x1)
    last_root = 0
    i = 1
    while i <= N:
        c = (x0 + x1) / 2
        fc = F(c)
        last_root = c
        if fc == 0:
            break
        if fa * fc < 0:
            x1 = c
            fb = fc
        else:
            x0 = c
            fa = fc
        i += 1

    return last_root


def newtonMethod(F, e, x0):
    # initializing Derivatives
    f_prime1 = get_prime()
    f_prime2 = get_prime2()
    # initializing first value
    x = x0
    while F(x) * f_prime2(x) <= 0:
        x += 0.1 # random step
    N = 1
    h = F(x) / f_prime1(x)
    while abs(h) >= 1 / 2 * e:
        h = F(x) / f_prime1(x)
        N += 1
        x = x - h
    print("TIMES: ", N)
    return x


def secantMethod(F, e, x0, x1):
    N = 0
    h = (F(x1) * (x1 - x0)) / (F(x1) - F(x0))
    x = x1 - h
    while abs(h) >= 1 / 2 * e:  # xn+1 - xn < 1/2*10^-5
        x0 = x1
        x1 = x
        h = (F(x1) * (x1 - x0)) / (F(x1) - F(x0))
        x = x - h
        N += 1
    print(N)
    return x


# defining f(x) with lower bound x = -2 and upper bound x = 2:

f = lambda x: math.e ** (math.sin(x) ** 3) + x ** 6 - 2 * x ** 4 - x ** 3 - 1

tol = 10 ** (-5)

print("BISECTION: ", bisectionMethod(f, tol, -2, -1))
print("BISECTION: ", bisectionMethod(f, tol, 1, 2))
print("BISECTION: ", bisectionMethod(f, tol, -0.5, 0.5))

# We need only one initial spot for newton-Raphson
print("NEWTON-RAPHSON: ", newtonMethod(f, tol, -2))
print("NEWTON-RAPHSON: ", newtonMethod(f, tol, 1))
print("NEWTON-RAPHSON: ", newtonMethod(f, tol, 0.5))

print("SECANT: ", secantMethod(f, tol, -2, -1))
print("SECANT: ", secantMethod(f, tol, -0.5, 0.5))
print("SECANT: ", secantMethod(f, tol, 1.4, 2))
