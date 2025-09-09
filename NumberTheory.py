import numpy as np
from Prime_All import primitiveRoot_prime

def gcd(a, b):
    """Tìm ước chung lớn nhất của 2 số a và b"""
    r = a % b
    while (r != 0):
        a = b
        b = r
        r = a % b
    return b

def modulo(number, mod):
    """Tính số dư của number chia cho mod"""
    return number % mod

def moduloPower(number, power, mod):
    """Tính modulo của number lũy thừa power
        Hay tính a^b mod m"""
    x = 1
    if (number < 0 and power < 0):
        number = -number
        power = -power
    elif number >= 0 and power < 0:
        number = inverseModulo(number, mod)
        power = -power

    tmp = number
    while power != 0:
        r = power % 2
        if r == 1:
            x = (x * tmp) % mod
        tmp = (tmp * tmp) % mod
        power = power // 2

    while x < 0:
        x = x + mod

    return x

def inverseModulo(a, mod):
    """Tìm nghịch đảo của a theo modulo mod
        Hay tìm x để a * x = 1 (mod m)
        Hay x = a^(-1) mod m"""
    a0 = mod
    b0 = a
    s0, s = 1, 0
    t0, t = 0, 1
    q = a0 // b0
    r = a0 - q * b0

    while r > 0:
        temp = t0 - q * t
        t0 = t
        t = temp
        temp = s0 - q * s
        s0 = s
        s = temp
        a0 = b0
        b0 = r
        q = a0 // b0
        r = a0 - q * b0

    if b0 != 1:
        return None
    else:
        while t < 0:
            t = t + mod
        return t

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

def primitiveRoot(p):
    """Tìm căn nguyên thủy của p"""
    n = p - 1
    # Tính nhân tử của p - 1
    primes = primitiveRoot_prime(n)
    factor = [p for p in primes if n % p == 0]

    first = 0
    factor_count = len(factor)
    for alpha in np.arange(2, n):
        if gcd(alpha, p) != 1:
            continue
        for i, pi in enumerate(factor):
            x = n / pi
            k = moduloPower(alpha, x, p)
            if k == 1:
                break
            if i == factor_count - 1:
                first = alpha

        if first != 0:
            break

    res = np.array([first])
    for i in np.arange(2, n):
        if gcd(i, n) == 1:
            res = np.append(res, moduloPower(first, i, p))

    res.sort()

    return res

# print(np.random.choice(primitiveRoot(127)))
#
print(moduloPower(32343, 4575657222473777152, 12323))

print(inverseModulo(7, 26))