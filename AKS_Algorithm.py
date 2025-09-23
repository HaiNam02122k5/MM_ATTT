import math
from sympy import Poly, symbols
from sympy.polys.domains import ZZ

""" Kiểm tra nguyên tố AKS - Thuật toán chính xác
    Có ý nghĩa trên lý thuyết - Thực tế quá phức tạp để triển khai"""
def gcd_in_prime(a, b):
    """Tính ước chung lớn nhất của a và b (dùng thuật toán Euclid)"""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def is_perfect_power(n):
    """Kiểm tra xem n có phải là lũy thừa hoàn hảo (m^k với k >= 2) hay không"""
    if n < 2:
        return False
    for b in range(2, int(math.log2(n)) + 1):
        root = round(n ** (1.0 / b))
        if root ** b == n:
            return True
    return False


def multiplicative_order(n, r):
    """Tìm trật tự nhân của n trong Z/rZ"""
    if gcd_in_prime(n, r) != 1:
        return 0
    k = 1
    temp = n % r
    while temp != 1:
        temp = (temp * n) % r
        k += 1
        if k > r:  # Tránh vòng lặp vô hạn
            return float('inf')
    return k


def find_r(n):
    """Tìm r sao cho trật tự nhân của n trong Z/rZ lớn hơn log^2(n)"""
    log_n = math.log2(n)
    max_r = int(log_n ** 2) + 1
    for r in range(2, max_r + 1):
        if gcd_in_prime(r, n) != 1:
            continue
        if multiplicative_order(n, r) > log_n ** 2:
            return r
    return max_r


def poly_check(n, r, a):
    """Kiểm tra (x - a)^n ≡ x^n - a (mod n, x^r - 1)"""
    x = symbols('x')
    # Đa thức (x - a)^n
    p1 = Poly((x - a) ** n, x, domain=ZZ)
    # Đa thức x^n - a
    p2 = Poly(x ** n - a, x, domain=ZZ)
    # Lấy modulo theo x^r - 1
    modulus = Poly(x ** r - 1, x, domain=ZZ)
    p1_mod = p1 % modulus
    p2_mod = p2 % modulus
    # Lấy modulo n cho các hệ số
    p1_coeffs = [c % n for c in p1_mod.all_coeffs()]
    p2_coeffs = [c % n for c in p2_mod.all_coeffs()]
    # Tái tạo đa thức với các hệ số đã modulo
    p1_mod = Poly(p1_coeffs, x, domain=ZZ)
    p2_mod = Poly(p2_coeffs, x, domain=ZZ)
    # So sánh hai đa thức
    return p1_mod == p2_mod


def is_prime_aks(n):
    """Kiểm tra tính nguyên tố của n bằng thuật toán AKS"""
    # Bước 1: Kiểm tra cơ bản
    if n <= 1:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    if is_perfect_power(n):
        return False

    # Bước 2: Tìm r
    r = find_r(n)

    # Bước 3: Kiểm tra ước chung
    for a in range(2, min(r + 1, n)):
        g = gcd_in_prime(a, n)
        if 1 < g < n:
            return False

    # Bước 4: Nếu n <= r, n là số nguyên tố
    if n <= r:
        return True

    # Bước 5: Kiểm tra đa thức
    for a in range(1, int(math.sqrt(r) * math.log2(n)) + 1):
        if not poly_check(n, r, a):
            return False

    # Bước 6: Nếu vượt qua tất cả, n là số nguyên tố
    return True
