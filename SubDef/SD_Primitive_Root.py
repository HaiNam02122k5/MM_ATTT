from random import randint
import gmpy2
from gmpy2 import mpz
import multiprocessing as mp


def pollard_rho(n, seed=2, max_steps=10000):
    """Thuật toán Pollard's Rho để tìm một thừa số nguyên tố của n."""
    n = mpz(n)
    if n == 1:
        return None
    if gmpy2.is_prime(n):
        return n
    if gmpy2.is_even(n):
        return mpz(2)

    def f(x):
        return (x * x + 1) % n

    x = mpz(seed)
    y = x
    c = mpz(randint(1, n - 1))
    d = mpz(1)
    steps = 0

    while d == 1 and steps < max_steps:
        x = f(x)
        y = f(f(y))
        d = gmpy2.gcd(abs(x - y), n)
        if d == n:
            return None  # Thử lại với seed khác
        steps += 1

    if gmpy2.is_prime(d):
        return d
    return pollard_rho(d, seed + 1, max_steps)  # Đệ quy trên thừa số tìm được

def factorize(n):
    """Phân tích n thành các nhân tử nguyên tố sử dụng Pollard's Rho và trial division."""
    # start_time = time.time()
    n = mpz(n)
    factors = []

    # Xử lý thừa số 2 trước
    while gmpy2.is_even(n):
        factors.append(mpz(2))
        n //= 2

    # Thử chia cho các số nhỏ (tối ưu cho các thừa số nhỏ)
    for i in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        while n % i == 0:
            factors.append(mpz(i))
            n //= i

    # Nếu n vẫn lớn, dùng Pollard's Rho
    while n > 1:
        if gmpy2.is_prime(n):
            factors.append(n)
            break
        factor = pollard_rho(n)
        if factor is None:
            # Nếu Pollard's Rho thất bại, thử seed khác hoặc tăng max_steps
            factor = pollard_rho(n, seed=randint(1, 100), max_steps=100000)
        if factor:
            factors.append(factor)
            n //= factor
        else:
            # Nếu không tìm được, có thể n là nguyên tố hoặc cần thuật toán khác
            factors.append(n)
            break

    # Trả về danh sách các nhân tử nguyên tố duy nhất
    factors = list(set(factors))
    # print(f"Factoring took {time.time() - start_time:.2f} seconds. Factors: {factors}")
    return factors

def is_primitive_root(g, p, factors):
    """Kiểm tra xem g có phải là căn nguyên thủy của p không."""
    g = mpz(g)
    p = mpz(p)
    for q in factors:
        exponent = (p - 1) // q
        if gmpy2.powmod(g, exponent, p) == 1:
            return False
    return True

def worker(args):
    """Worker cho parallel checking."""
    g, p, factors = args
    return g if is_primitive_root(g, p, factors) else None

def find_primitive_root(p, use_parallel=True):
    """Tìm căn nguyên thủy của số nguyên tố p."""
    p = mpz(p)

    # Kiểm tra xem p có phải là số nguyên tố
    if not gmpy2.is_prime(p):
        return "p must be prime"

    # Phân tích p-1 thành các nhân tử nguyên tố
    p_minus_1 = p - 1
    factors = factorize(p_minus_1)

    # Thử các g nhỏ trước
    small_gs = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    for g in small_gs:
        if is_primitive_root(g, p, factors):
            return mpz(g)

    # Thử ngẫu nhiên với parallel nếu bật
    trials = 100000
    if use_parallel:
        pool = mp.Pool(mp.cpu_count())  # Sử dụng tất cả lõi của Ryzen
        candidates = [(randint(2, int(p - 1)), p, factors) for _ in range(trials)]
        results = pool.map(worker, candidates)
        pool.close()
        for res in results:
            if res is not None:
                return mpz(res)
    else:
        # Sequential
        for _ in range(trials):
            g = randint(2, int(p - 1))
            if is_primitive_root(g, p, factors):
                return mpz(g)

    return None, "No primitive root found after trials (increase trials or check factoring)"

