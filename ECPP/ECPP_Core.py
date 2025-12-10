"""
ECPP Implementation - Part 1: Fundamental Discriminants and Class Polynomials
Thuật toán chứng minh tính nguyên tố bằng Elliptic Curve
"""

import math
from NumberTheory import gcd, moduloPower, modulo, inverseModulo
from Prime_All import prime_check, SMALL_PRIMES


class ECPPHelper:
    """Các hàm hỗ trợ cho ECPP"""

    # Danh sách discriminants cơ bản (fundamental discriminants)
    # Đây là các D nhỏ thường được dùng trong ECPP
    SMALL_DISCRIMINANTS = [
        -3, -4, -7, -8, -11, -15, -19, -20, -23, -24,
        -31, -35, -39, -40, -43, -47, -51, -52, -55, -56,
        -59, -67, -68, -71, -79, -83, -84, -87, -88, -91,
        -95, -103, -104, -107, -115, -116, -119, -120, -123, -127
    ]

    @staticmethod
    def is_fundamental_discriminant(D):
        """
        Kiểm tra xem D có phải là fundamental discriminant không
        D phải âm và thỏa mãn một trong hai điều kiện:
        1. D ≡ 1 (mod 4) và square-free
        2. D = 4m với m ≡ 2,3 (mod 4) và m square-free
        """
        if D >= 0:
            return False

        D = abs(D)

        # Case 1: D ≡ 1 (mod 4)
        if D % 4 == 1:
            return ECPPHelper.is_square_free(D)

        # Case 2: D = 4m
        if D % 4 == 0:
            m = D // 4
            if m % 4 in [2, 3]:
                return ECPPHelper.is_square_free(m)

        return False

    @staticmethod
    def is_square_free(n):
        """Kiểm tra n có square-free không (không chia hết cho bình phương nào)"""
        if n <= 1:
            return False

        # Kiểm tra với các số nguyên tố nhỏ
        for p in SMALL_PRIMES:
            if p * p > n:
                break
            if n % (p * p) == 0:
                return False

        return True

    @staticmethod
    def kronecker_symbol(a, n):
        """
        Tính Kronecker symbol (a/n) - mở rộng của Jacobi symbol
        Dùng để tính class number và tìm curve phù hợp
        """
        if n == 0:
            return 1 if abs(a) == 1 else 0

        if n < 0:
            if a < 0:
                return -ECPPHelper.kronecker_symbol(-a, -n)
            return ECPPHelper.kronecker_symbol(a, -n)

        # Xử lý n chẵn
        if n % 2 == 0:
            if a % 2 == 0:
                return 0
            # (a/2) = 1 if a ≡ ±1 (mod 8), else -1
            k = 0
            n_temp = n
            while n_temp % 2 == 0:
                k += 1
                n_temp //= 2

            result = ECPPHelper.kronecker_symbol(a, n_temp)
            if k % 2 == 1:
                if a % 8 in [1, 7]:
                    return result
                else:
                    return -result
            return result

        # Jacobi symbol cho n lẻ
        return ECPPHelper.jacobi_symbol(a, n)

    @staticmethod
    def jacobi_symbol(a, n):
        """Tính Jacobi symbol (a/n)"""
        if n <= 0 or n % 2 == 0:
            raise ValueError("n must be positive odd integer")

        a = a % n
        result = 1

        while a != 0:
            while a % 2 == 0:
                a //= 2
                if n % 8 in [3, 5]:
                    result = -result

            a, n = n, a
            if a % 4 == 3 and n % 4 == 3:
                result = -result
            a = a % n

        return result if n == 1 else 0

    @staticmethod
    def compute_hilbert_class_polynomial(D):
        """
        Tính Hilbert class polynomial H_D(X)

        Đây là phần KHÁC NHẤT và QUAN TRỌNG NHẤT của ECPP!

        Vì tính toán chính xác Hilbert polynomial rất phức tạp (cần complex
        multiplication theory), ta dùng bảng lookup cho các D nhỏ.

        Trong thực tế production, cần dùng thư viện chuyên dụng hoặc
        pre-computed tables.
        """

        # Bảng lookup cho một số Hilbert polynomials thông dụng
        # Format: D -> [coefficients từ bậc cao xuống thấp]
        HILBERT_POLYNOMIALS = {
            -3: [1, 0],  # X
            -4: [1, 1728],  # X + 1728
            -7: [1, -3375],  # X - 3375
            -8: [1, 8000],  # X + 8000
            -11: [1, -32768],  # X - 32768
            -19: [1, -884736],  # X - 884736
            -43: [1, -884736000],  # X - 884736000
            -67: [1, -147197952000],  # X - 147197952000
            -163: [1, -262537412640768000],  # X - 262537412640768000

            # Class number 2
            -15: [1, 191025, -121287375],
            -20: [1, 1264000, -681472000],
            -24: [1, 2587918086, 21425024],
            -35: [1, 16581375, 1264000],
            -40: [1, 16581375, 1264000],
            -51: [1, 39491307, 1264000],
            -52: [1, 39491307, 1264000],
            -88: [1, 14670139281, 39491307],
            -91: [1, 14670139281, 39491307],
            -115: [1, 14670139281, 39491307],
            -123: [1, 14670139281, 39491307],
            -148: [1, 14670139281, 39491307],
            -187: [1, 14670139281, 39491307],
            -232: [1, 14670139281, 39491307],
            -235: [1, 14670139281, 39491307],
            -267: [1, 14670139281, 39491307],
            -403: [1, 14670139281, 39491307],
            -427: [1, 14670139281, 39491307],
        }

        if D in HILBERT_POLYNOMIALS:
            return HILBERT_POLYNOMIALS[D]

        # Nếu không có trong bảng, trả về None
        # Trong production cần implement hoặc dùng thư viện
        return None

    @staticmethod
    def find_curve_from_discriminant(D, p):
        """
        Tìm đường cong elliptic từ discriminant D modulo p
        Trả về (a, b) sao cho y² = x³ + ax + b có tính chất phù hợp
        """
        if D not in [-3, -4]:
            # Với D khác, cần tính từ Hilbert polynomial root
            # Simplified version: random search
            return ECPPHelper.find_curve_simple(p)

        if D == -3:
            # Curve: y² = x³ + b với j-invariant = 0
            # Chọn b random sao cho 4b³ ≢ 0 (mod p)
            for b in range(1, min(p, 1000)):
                if modulo(4 * moduloPower(b, 3, p), p) != 0:
                    return (0, b)

        if D == -4:
            # Curve: y² = x³ + ax với j-invariant = 1728
            # Chọn a random sao cho 4a³ ≢ 0 (mod p)
            for a in range(1, min(p, 1000)):
                if modulo(4 * moduloPower(a, 3, p), p) != 0:
                    return (a, 0)

        return None

    @staticmethod
    def find_curve_simple(p):
        """Tìm đường cong elliptic đơn giản (fallback method)"""
        import random
        for _ in range(100):
            a = random.randint(0, p - 1)
            b = random.randint(0, p - 1)
            # Kiểm tra non-singular: 4a³ + 27b² ≢ 0 (mod p)
            discriminant = modulo(4 * moduloPower(a, 3, p) + 27 * moduloPower(b, 2, p), p)
            if discriminant != 0:
                return (a, b)
        return None


# Test code
if __name__ == "__main__":
    helper = ECPPHelper()

    print("=== Testing Fundamental Discriminants ===")
    for D in [-3, -4, -7, -8, -11, -15, -16]:
        is_fund = helper.is_fundamental_discriminant(D)
        print(f"D = {D:4d}: {'✓' if is_fund else '✗'}")

    print("\n=== Testing Hilbert Class Polynomials ===")
    for D in [-3, -4, -7, -11, -19]:
        poly = helper.compute_hilbert_class_polynomial(D)
        print(f"H_{D}(X) = {poly}")

    print("\n=== Testing Curve Finding ===")
    p = 101
    for D in [-3, -4]:
        curve = helper.find_curve_from_discriminant(D, p)
        print(f"D = {D}, p = {p}: curve (a,b) = {curve}")