import math

import gmpy2
import numpy as np
from gmpy2 import mpz

from SubDef.SD_Primitive_Root import find_primitive_root


def gcd(a, b):
    """Tìm ước chung lớn nhất của 2 số a và b"""
    if a > 2**128 or b > 2**128:
        a = mpz(a)
        b = mpz(b)
        return gmpy2.gcd(a, b)
    else:
        return math.gcd(a, b)

def modulo(number, mod):
    """Tính số dư của number chia cho mod"""
    if number > 2 ** 128 or mod > 2 ** 128:
        number = mpz(number)
        mod = mpz(mod)
    return number % mod

def moduloPower(number, power, mod):
    """Tính modulo của number lũy thừa power
        Hay tính a^b mod m"""
    if number > 2 ** 128 or mod > 2 ** 128 or power > 2 ** 128 :
        return gmpy2.powmod(number, power, mod)
    else:
        return pow(number, power, mod)

def inverseModulo(a, mod):
    """Tìm nghịch đảo của a theo modulo mod
        Hay tìm x để a * x = 1 (mod m)
        Hay x = a^(-1) mod m"""
    if a > 2 ** 128 or mod > 2 ** 128:
        return int(gmpy2.invert(a, mod))
    else:
        return pow(a, -1, mod)


def linearCongruence(a, b, m):
    """Giải phương trình đồng dư ax = b (mod m)"""
    d = gcd(a, m)

    a = a // d
    b = b // d
    m = m // d

    inv = inverseModulo(a, m)
    if inv is None:
        return None

    x = modulo(b * inv, m)
    res = np.empty(d, dtype=int)
    for i in np.arange(d-1):
        res = np.add(res, x + i * m)
    return res

def part_primitive_root(p, limit = 50, use_parallel=True):
    """Tìm căn nguyên thủy của số nguyên tố p"""
    firstNumber = find_primitive_root(p, use_parallel)
    res = [firstNumber]
    n = p - 1

    for i in range(2, p):
        if limit == 1:
            break
        if gcd(i, n) == 1:
            res.append(i)
            limit = limit - 1

    return res

