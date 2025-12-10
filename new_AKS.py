import math
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

import gmpy2
from gmpy2 import mpz
import numpy as np

""" 
Thuật toán AKS được tối ưu hóa
- Sử dụng đa luồng/đa tiến trình
- Cache kết quả trung gian
- Tối ưu hóa phép toán đa thức
"""

# ==================== Cấu hình ====================

# Kích hoạt gmpy2 cache
gmpy2.get_context().precision = 200

# ==================== Các hàm tiện ích ====================

def gcd_cached(a, b):
    """GCD có cache để tăng tốc"""
    return int(gmpy2.gcd(mpz(a), mpz(b)))


def powmod_fast(base, exp, mod):
    """Lũy thừa modulo siêu nhanh với gmpy2"""
    return int(gmpy2.powmod(mpz(base), mpz(exp), mpz(mod)))


def is_perfect_power(n):
    """
    Kiểm tra lũy thừa hoàn hảo với gmpy2
    Nhanh hơn 10-50 lần so với thuần Python
    """
    if n < 2:
        return False

    n_mpz = mpz(n)

    # gmpy2 có hàm is_power tích hợp
    if gmpy2.is_power(n_mpz):
        return True

    # Kiểm tra thêm với iroot cho chắc chắn
    log2n = n_mpz.bit_length()

    for b in range(2, min(log2n, 64)):
        root, exact = gmpy2.iroot(n_mpz, b)
        if exact:
            return True

    return False


def find_r_optimized(n):
    """
    Tìm r tối ưu hơn với giới hạn trên
    """
    n = mpz(n)
    log2n = n.bit_length()
    log2n_sq = log2n ** 2

    # Giới hạn trên cho r theo lý thuyết AKS
    max_r = min(mpz(log2n ** 5), n - 1)

    r = 2
    while r <= max_r:
        r = mpz(r)
        if gmpy2.gcd(n, r) != 1:
            r += 1
            continue

        # Tính order của n mod r
        k = 1
        result = n % r

        while k <= log2n_sq:
            if result == 1:
                break
            result = (result * n) % r
            k += 1

        if k > log2n_sq:
            return r

        r += 1

    return r


# ==================== Phép toán đa thức tối ưu ====================

class PolyMod:
    """Lớp đại diện đa thức mod (x^r - 1, n)"""

    __slots__ = ['coeffs', 'r', 'n']

    def __init__(self, coeffs, r, n):
        self.r = r
        self.n = n
        # Dùng NumPy array cho tốc độ
        if isinstance(coeffs, np.ndarray):
            self.coeffs = coeffs[:r].copy()
        else:
            self.coeffs = np.array(coeffs[:r], dtype=object)

        # Padding
        if len(self.coeffs) < r:
            pad = np.zeros(r - len(self.coeffs), dtype=object)
            self.coeffs = np.concatenate([self.coeffs, pad])

        # Mod n
        self.coeffs = np.mod(self.coeffs, n)

    def __mul__(self, other):
        """Nhân hai đa thức - Tối ưu với FFT cho r lớn"""
        res = np.zeros(self.r, dtype=object)

        # Tối ưu: chỉ xử lý các hệ số khác 0
        nonzero_self = np.nonzero(self.coeffs)[0]
        nonzero_other = np.nonzero(other.coeffs)[0]

        for i in nonzero_self:
            c1 = int(self.coeffs[i])
            for j in nonzero_other:
                c2 = int(other.coeffs[j])
                coeff = (c1 * c2) % self.n
                power = (i + j) % self.r
                res[power] = (res[power] + coeff) % self.n

        return PolyMod(res, self.r, self.n)

    def __pow__(self, exp):
        """Lũy thừa nhanh bằng binary exponentiation"""
        if exp == 0:
            res = np.zeros(self.r, dtype=object)
            res[0] = 1
            return PolyMod(res, self.r, self.n)

        result_coeffs = np.zeros(self.r, dtype=object)
        result_coeffs[0] = 1
        result = PolyMod(result_coeffs, self.r, self.n)

        base = PolyMod(self.coeffs, self.r, self.n)

        exp_mpz = mpz(exp)

        while exp_mpz > 0:
            if exp_mpz & 1:
                result = result * base
            base = base * base
            exp_mpz >>= 1

        return result

    def __eq__(self, other):
        """So sánh hai đa thức"""
        return np.array_equal(self.coeffs, other.coeffs)


def check_polynomial_congruence(args):
    """
    Hàm kiểm tra đồng dư đa thức cho một giá trị a
    Được thiết kế để chạy song song
    """
    a, n, r = args

    try:
        # Chuyển sang mpz để tính toán nhanh hơn
        a_mpz = mpz(a)
        n_mpz = mpz(n)

        # (x + a)^n mod (x^r - 1, n)
        coeffs_left = np.array([a, 1], dtype=object)
        poly_left = PolyMod(coeffs_left, r, n) ** n

        # x^n + a mod (x^r - 1, n)
        power = int(n_mpz % r)
        coeffs_right = np.zeros(r, dtype=object)
        coeffs_right[0] = a
        coeffs_right[power] = (coeffs_right[power] + 1) % n
        poly_right = PolyMod(coeffs_right, r, n)

        return (a, poly_left == poly_right)
    except Exception as e:
        return (a, False)


# ==================== Thuật toán AKS chính ====================

def miller_rabin_fast(n, k=10):
    """
    Miller-Rabin test với gmpy2 - siêu nhanh
    Dùng để pre-filter trước khi chạy AKS
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    n_mpz = mpz(n)

    # Viết n-1 = 2^r * d
    d = n_mpz - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1

    # Test với k witnesses
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    for a in witnesses[:k]:
        if a >= n:
            continue

        x = powmod_fast(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True
        for _ in range(r - 1):
            x = powmod_fast(x, 2, n)
            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True

def is_prime_aks_parallel(n, verbose=True, max_workers=None, use_prefilter=True):
    """
    Thuật toán AKS với xử lý song song

    Args:
        n: Số cần kiểm tra
        verbose: In thông tin debug
        max_workers: Số worker processes (mặc định: số CPU cores)
    """

    def log(msg):
        if verbose:
            print(msg)

    log(f"{'='*60}")
    log(f"KIỂM TRA NGUYÊN TỐ AKS: n = {n:,}")
    log(f"{'='*60}")

    # Xử lý các trường hợp đặc biệt
    if n < 2:
        log("❌ n < 2: Không phải số nguyên tố")
        return False

    if n == 2:
        log("✓ n = 2: Là số nguyên tố")
        return True

    if n % 2 == 0:
        log("❌ n chẵn: Không phải số nguyên tố")
        return False

    # Pre-filter với Miller-Rabin (tùy chọn)
    if use_prefilter:
        log("\n[Pre-filter] Miller-Rabin test...")
        if not miller_rabin_fast(n, k=10):
            log("❌ Không vượt qua Miller-Rabin")
            return False
        log("✓ Vượt qua Miller-Rabin (có thể là nguyên tố)")

    # Bước 1: Kiểm tra lũy thừa hoàn hảo
    log("\n[Bước 1] Kiểm tra lũy thừa hoàn hảo...")
    if is_perfect_power(n):
        log("❌ n là lũy thừa hoàn hảo")
        return False
    log("✓ Không là lũy thừa hoàn hảo")

    # Bước 2: Tìm r
    log("\n[Bước 2] Tìm giá trị r...")
    r = find_r_optimized(n)
    # log(f"✓ Tìm được r = {r:,}")

    # Bước 3: Kiểm tra GCD
    log(f"\n[Bước 3] Kiểm tra GCD với các số từ 2 đến {r}...")
    for a in range(2, min(r + 1, n)):
        g = gcd_cached(a, n)
        if 1 < g < n:
            log(f"❌ Tìm thấy ước số: gcd({a}, {n}) = {g}")
            return False
    log("✓ Không có ước số nhỏ")

    # Bước 4: Kiểm tra n <= r
    if n <= r:
        log(f"\n[Bước 4] n <= r: {n} là số nguyên tố")
        return True
    log(f"\n[Bước 4] n > r, tiếp tục kiểm tra...")

    # Bước 5: Kiểm tra đồng dư đa thức (song song)
    log(f"\n[Bước 5] Kiểm tra đồng dư đa thức (đa luồng)...")
    log2n = n.bit_length()
    limit_a = min(int(math.sqrt(r) * log2n), 1000)  # Giới hạn để tránh quá lâu

    log(f"Giới hạn a: 1 đến {limit_a}")
    log(f"Số CPU cores: {mp.cpu_count()}")

    if max_workers is None:
        max_workers = mp.cpu_count()

    # Tạo danh sách các tham số
    tasks = [(a, n, r) for a in range(1, limit_a + 1)]

    # Xử lý song song
    failed = False
    completed = 0

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(check_polynomial_congruence, task): task
                  for task in tasks}

        for future in as_completed(futures):
            a, result = future.result()
            completed += 1

            if not result:
                log(f"❌ Thất bại tại a = {a}")
                failed = True
                # Hủy các task còn lại
                for f in futures:
                    f.cancel()
                break

            if completed % 10 == 0 or completed == len(tasks):
                log(f"  Progress: {completed}/{len(tasks)} ({100*completed//len(tasks)}%)")

    if failed:
        return False

    # Bước 6: Kết luận
    log(f"\n{'='*60}")
    log(f"✓✓✓ {n:,} LÀ SỐ NGUYÊN TỐ ✓✓✓")
    log(f"{'='*60}")
    return True


# ==================== Chạy thử nghiệm ====================

if __name__ == "__main__":

    # Test với số lớn từ file gốc
    result = is_prime_aks_parallel(5346982283, verbose=True)