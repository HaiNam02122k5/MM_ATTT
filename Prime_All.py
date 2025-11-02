import math
import os
import random
import multiprocessing

import numpy as np
import gmpy2

from gmpy2 import mpz, random_state, mpz_random

"""
Sàng nguyên tố. Trả về dãy các số nguyên tố nhỏ hơn n sd khi nhỏ hơn 32 bit
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
    segment_size = max(sqrt_n, 10 ** 7)
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
        primes = [p for p in primes if p <= n and p > 1]  # Lọc các số nguyên tố <= n

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
    return [p for p in primes if p <= sqrt_n and p > 1]

# Pre-compute small primes cho trial division (dùng cho tất cả cases)
SMALL_PRIME_LIMIT = 10 ** 6
SMALL_PRIMES = simple_sieve(SMALL_PRIME_LIMIT)


def is_divisible_by_small_primes(n):
    """Chia thử để loại nhanh composites."""
    if n < 2:
        return True
    if n < 2 ** 128:  # Small/medium: Dùng Python thuần
        sqrt_n = int(math.sqrt(n)) + 1
        for p in SMALL_PRIMES:
            if p > sqrt_n:
                break
            if n % p == 0:
                return True
    else:  # Large: Dùng gmpy2.isqrt và mpz
        n = mpz(n)
        sqrt_n = gmpy2.isqrt(n) + 1
        for p in SMALL_PRIMES:
            if p > sqrt_n:
                break
            if n % p == 0:
                return True
    return False

# Miller-Rabin pure Python với fixed bases cho small/medium bits
def miller_rabin_pure(n):
    """Miller-Rabin deterministic cho n < 2^128."""
    if n < 2:
        return False
    if n in (2, 3, 5, 7, 11, 13, 17, 23):
        return True
    if n % 2 == 0:
        return False

    # Write n-1 as 2^r * s
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2

    # Fixed bases đủ cho n < 2^128 (dựa trên known deterministic sets)
    bases = [2, 3, 5, 7, 11, 13, 17, 23, 29, 31, 37]  # Chính xác 100% cho n < 3.317e38 (~128 bit)

    for a in bases:
        if a >= n:
            break
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


"""Kiểm tra số nguyên tố - Chia case"""
def prime_check(n):
    if n < 2 ** 128:  # Small/medium: Pure Python
        if is_divisible_by_small_primes(n):
            return False
        return miller_rabin_pure(n)
    else:  # Large: gmpy2
        n = mpz(n)
        return gmpy2.is_prime(n)

# Hàm check_candidate tách ra module level
def check_candidate_large(low, high, seed):
    """Check candidate cho large bits (gmpy2)."""
    rstate = random_state(seed)
    while True:
        candidate = mpz_random(rstate, high - low + 1) + low
        if not is_divisible_by_small_primes(int(candidate)):
            if prime_check(int(candidate)):
                return int(candidate)


"""Tạo số nguyên tố n bit - Chia case"""
def generate_prime_bit(bit, num_processes=8):
    """Tạo số nguyên tố n bit - Chia case."""
    low = 2 ** (bit - 1)
    high = 2 ** bit - 1

    if bit <= 128:  # Small/medium: Pure Python
        while True:
            candidate = random.randint(low, high)
            if not is_divisible_by_small_primes(candidate):
                if miller_rabin_pure(candidate):
                    return candidate
    else:  # Large: gmpy2 + parallel
        multiprocessing.set_start_method('spawn', force=True)  # Force spawn cho Windows
        low = mpz(low)
        high = mpz(high)

        # Windows cần if __name__ == '__main__' bảo vệ Pool
        # Nhưng vì đây là hàm, ta đảm bảo check_candidate_large ở module level
        with multiprocessing.Pool(processes=num_processes) as pool:
            seeds = [random.randint(1, 1000000) for _ in range(num_processes)]
            results = [pool.apply_async(check_candidate_large, args=(low, high, seed))
                       for seed in seeds]
            for res in results:
                try:
                    prime = res.get(timeout=5)
                    pool.terminate()
                    return prime
                except multiprocessing.TimeoutError:
                    pass
        return int(gmpy2.next_prime(low))


def generate_prime_in_range(low, high, num_processes=8):
    """Tạo số nguyên tố trong range - Chia case."""
    bit_est = math.log2(high)
    if bit_est <= 128:
        while True:
            candidate = random.randint(low, high)
            if not is_divisible_by_small_primes(candidate):
                if miller_rabin_pure(candidate):
                    return candidate
    else:
        low_mpz = mpz(low)
        high_mpz = mpz(high)

        with multiprocessing.Pool(processes=num_processes) as pool:
            seeds = [random.randint(1, 1000000) for _ in range(num_processes)]
            results = [pool.apply_async(check_candidate_large, args=(low_mpz, high_mpz, seed))
                       for seed in seeds]
            for res in results:
                try:
                    prime = res.get(timeout=5)
                    pool.terminate()
                    return prime
                except multiprocessing.TimeoutError:
                    pass
        return int(gmpy2.next_prime(low_mpz))

