import math
import os
import random

import numpy as np
from sympy import Poly, symbols
from sympy.polys.domains import ZZ

"""
Sàng nguyên tố. Trả về dãy các số nguyên tố nhỏ hơn n
"""
def simple_sieve(limit):
    """Tìm các số nguyên tố nhỏ hơn limit bằng Sàng Eratosthenes."""
    # Khởi tạo mảng boolean cho các số lẻ từ 3 đến n
    size = (limit + 1) // 2
    is_prime = np.ones(size, dtype=np.bool_)

    # Sàng Eratosthenes
    for i in range(int(np.sqrt(limit)) // 2):
        if is_prime[i]:  # Nếu 2*i + 3 là số nguyên tố
            p = 2 * i + 3
            # Đánh dấu các bội số của p từ p^2 trở đi
            start = (p * p - 3) // 2
            # Tránh truy cập ngoài mảng
            if start < size:
                is_prime[start::p] = False

    # Tạo danh sách các số nguyên tố với số đầu là 2
    primes = [2]
    primes.extend(2 * i + 3 for i in range(len(is_prime)) if is_prime[i])
    return primes


def segmented_sieve(n):
    """Tìm tất cả số nguyên tố nhỏ hơn n bằng Sàng phân đoạn."""
    if n < 2:
        return []

    # Tìm các số nguyên tố nhỏ hơn sqrt(n)
    sqrt_n = int(math.sqrt(n))
    base_primes = simple_sieve(sqrt_n)

    # Kích thước đoạn (có thể điều chỉnh tùy máy)
    segment_size = max(sqrt_n, 10 ** 6)
    primes = base_primes.copy()

    # Xử lý từng đoạn [low, high]
    for low in range(sqrt_n + 1, n + 1, segment_size):
        high = min(low + segment_size - 1, n)

        # Khởi tạo mảng cho đoạn [low, high]
        sieve = np.ones(high - low + 1, dtype=np.bool_)

        # Đánh dấu bội số của các số nguyên tố trong đoạn
        for p in base_primes:
            # Tìm bội số nhỏ nhất của p >= low
            start = (low + p - 1) // p * p
            if start < low:
                start += p
            sieve[start - low:high - low + 1:p] = False

        # Thêm các số nguyên tố trong đoạn vào danh sách
        for i in range(high - low + 1):
            if sieve[i]:
                primes.append(low + i)

    return primes


def segmented_sieve_to_file(n, output_file="primes.txt"):
    """Tìm tất cả số nguyên tố nhỏ hơn n và lưu vào file."""
    if n < 2:
        with open(output_file, 'w') as f:
            f.write("")
        return 0

    # Tìm các số nguyên tố nhỏ hơn sqrt(n)
    sqrt_n = int(math.sqrt(n))

    existing_max_prime = 0
    existing_primes = []

    # Kiểm tra file primes.txt
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        with open(output_file, 'r') as f:
            lines = f.readlines()
            if lines:
                try:
                    existing_max_prime = int(lines[0].strip())
                    existing_primes = [int(p.strip()) for p in lines[1:] if p.strip().isdigit()]
                except (ValueError, IndexError):
                    existing_max_prime = 0
                    existing_primes = []

    # Nếu n nhỏ hơn hoặc bằng số nguyên tố lớn nhất, giữ nguyên file
    if n <= existing_max_prime:
        primes = [p for p in existing_primes if p <= n]
        return len(primes), primes

    base_primes = simple_sieve(sqrt_n)
    primes = base_primes.copy()

    # Kích thước đoạn
    segment_size = max(sqrt_n, 10 ** 6)

    # Bắt đầu từ max(existing_max_prime + 1, sqrt_n + 1)
    start_low = existing_max_prime
    for low in range(start_low, n + 1, segment_size):
        high = min(low + segment_size - 1, n)
        sieve = np.ones(high - low + 1, dtype=np.bool_)

        for p in base_primes:
            start = (low + p - 1) // p * p
            if start < low:
                start += p
            sieve[start - low:high - low + 1:p] = False

        for i in range(high - low + 1):
            if sieve[i]:
                primes.append(low + i)

    # Kết hợp các số nguyên tố hiện có và mới
    if existing_primes:
        primes = sorted(list(set(existing_primes + primes)))  # Loại bỏ trùng lặp và sắp xếp
        primes = [p for p in primes if p <= n]  # Lọc các số nguyên tố <= n

    prime_count = len(primes)

    # Ghi vào file
    with open(output_file, 'w') as f:
        if primes:
            f.write(f"{primes[-1]}\n")
        for p in primes:
            f.write(f"{p}\n")

    return prime_count, primes

def primitiveRoot_prime(n):
    """ Sử dụng cho tìm căn nguyên thủy
        Tìm các số nguyên tố nhỏ hơn hoặc bằng căn n"""
    prime_count, primes = segmented_sieve_to_file(n)
    sqrt_n = int(math.sqrt(n))
    return [p for p in primes if p <= sqrt_n]

""" Kiểm tra số nguyên tố bằng Miller Rabin"""
def Miller_Rabin_algorithm(n):
    """ Kiểm tra số nguyên tố bằng thuật toán Miller Rabin
        Thuật toán "Xác suất" có thể sai nhưng rất nhỏ """
    if type(n) != int:
        return False

    if n == 2 or n == 3 or n == 5 or n == 7:
        return True

    if n % 2 == 0:
        return False

    k = 0
    m = n - 1

    while m % 2 == 0:
        m = m / 2
        k = k + 1

    def test(a, n, k, m):

        mod = pow(int(a), int(m), int(n))
        if (mod == 1 or mod == n - 1):
            return True

        for i in range(1, k):
            mod = (mod * mod) % n
            if mod == n - 1:
                return True

        return False


    repeatTime = 10
    for i in np.arange(repeatTime):
        a = np.random.rand() % (n - 3) + 2
        if not test(a, n, k, m):
            return False

    return True

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

""" Kiểm tra số nguyên tố - Chưa sử dụng AKS"""
def prime_check(n):
    """ Kiểm tra số nguyên tố """
    if Miller_Rabin_algorithm(n) == False:
        return False

    # if is_prime_aks(n) == False:
    #     return False

    return True

""" Tạo số nguyên tố n bit"""
def general_prime_bit(bit):
    """ Tạo một số nguyên tố n bit"""
    res = 0
    while True:
        number = random.getrandbits(bit)
        if prime_check(number) == True:
            res = number
            break

    return res

def general_prime_in_range(low, high):
    res = 0
    while True:
        number = random.randint(low, high)
        if prime_check(number) == True:
            res = number
            break

    return res


# # Kiểm tra ví dụ
# test_numbers = 13248734287293484729327
# print(is_prime_aks(test_numbers))

# # Kiểm tra ví dụ
# print(general_prime(17))
#
# print(is_prime_aks(66337))

# # Ví dụ sử dụng
# n = 2**20
# import time
# start_time = time.time()
# count, primes = segmented_sieve_to_file(n, "primes.txt")
# end_time = time.time()
#
# print(f"Số lượng số nguyên tố nhỏ hơn {n}: {count}")
# print(f"Thời gian thực thi: {end_time - start_time:.4f} giây")
# print("Đã lưu các số nguyên tố vào primes.txt")


# # Ví dụ sử dụng
# n = 133  # n > 2^32
# import time
#
# start_time = time.time()
# primes = segmented_sieve(n)
# end_time = time.time()
#
# print(f"Số lượng số nguyên tố nhỏ hơn {n}: {len(primes)}")
# print(f"Thời gian thực thi: {end_time - start_time:.4f} giây")
# print(f"Một vài số nguyên tố đầu tiên: {primes[:100]}")
# print(f"Số nguyên tố cuối cùng: {primes[-1]}")