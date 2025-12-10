"""
ECPP Utilities - Cornacchia Algorithm và các hàm hỗ trợ
Đây là các thuật toán quan trọng cho CM (Complex Multiplication)
"""

from NumberTheory import gcd, moduloPower, modulo, inverseModulo
from ECPP_Types import ECPPCertificate  # ✓ Import từ file types riêng
import random


class ECPPUtils:
    """Các utility functions cho ECPP"""

    @staticmethod
    def tonelli_shanks(n, p):
        """
        Thuật toán Tonelli-Shanks: Tìm r sao cho r² ≡ n (mod p)

        Args:
            n: số cần tìm căn
            p: số nguyên tố (modulus)

        Returns:
            r sao cho r² ≡ n (mod p), hoặc None nếu không tồn tại
        """
        # Kiểm tra n có phải quadratic residue không
        if moduloPower(n, (p - 1) // 2, p) != 1:
            return None

        # Trường hợp đặc biệt: p ≡ 3 (mod 4)
        if p % 4 == 3:
            r = moduloPower(n, (p + 1) // 4, p)
            if modulo(r * r, p) == modulo(n, p):
                return r
            return None

        # Trường hợp tổng quát: p ≡ 1 (mod 4)
        # Viết p - 1 = Q * 2^S với Q lẻ
        Q = p - 1
        S = 0
        while Q % 2 == 0:
            Q //= 2
            S += 1

        # Tìm quadratic non-residue z
        z = 2
        while moduloPower(z, (p - 1) // 2, p) != p - 1:
            z += 1

        # Khởi tạo
        M = S
        c = moduloPower(z, Q, p)
        t = moduloPower(n, Q, p)
        R = moduloPower(n, (Q + 1) // 2, p)

        while True:
            if t == 0:
                return 0
            if t == 1:
                return R

            # Tìm i nhỏ nhất sao cho t^(2^i) = 1
            i = 1
            temp = modulo(t * t, p)
            while temp != 1 and i < M:
                temp = modulo(temp * temp, p)
                i += 1

            # Update
            b = moduloPower(c, 1 << (M - i - 1), p)
            M = i
            c = modulo(b * b, p)
            t = modulo(t * c, p)
            R = modulo(R * b, p)

    @staticmethod
    def cornacchia(D, p):
        """
        Thuật toán Cornacchia: Giải phương trình x² + |D|y² = 4p

        Đây là thuật toán QUAN TRỌNG NHẤT cho Complex Multiplication!
        Dùng để tính order của elliptic curve.

        Args:
            D: discriminant (âm)
            p: số nguyên tố

        Returns:
            (x, y) nếu tìm được, None nếu không
        """
        if D >= 0:
            return None

        # Special cases
        if p == 2:
            if D == -1:
                return (2, 2) if 8 % 4 == 0 else None
            return None

        # Kiểm tra (D/p) = 1 (Kronecker symbol)
        from ECPP_Core import ECPPHelper
        if ECPPHelper.kronecker_symbol(D, p) != 1:
            return None

        # Tìm x₀ sao cho x₀² ≡ D (mod p)
        x0 = ECPPUtils.tonelli_shanks(modulo(D, p), p)
        if x0 is None:
            return None

        # Đảm bảo x₀ lẻ (nếu D lẻ)
        if abs(D) % 2 == 1 and x0 % 2 == 0:
            x0 = p - x0

        # Euclidean algorithm với bound
        a = 2 * p
        b = x0
        limit = int((4 * p) ** 0.5)

        while b > limit:
            a, b = b, a % b

        # x = b
        x = b

        # Tính y từ phương trình
        # 4p = x² + |D|y²
        # y² = (4p - x²) / |D|
        remainder = 4 * p - x * x

        if remainder < 0 or remainder % abs(D) != 0:
            return None

        y_squared = remainder // abs(D)
        y = int(y_squared ** 0.5)

        # Verify
        if y * y != y_squared:
            return None

        if x * x + abs(D) * y * y != 4 * p:
            return None

        return (x, y)

    @staticmethod
    def compute_trace_from_cornacchia(D, p):
        """
        Tính trace của Frobenius từ kết quả Cornacchia

        Nếu x² + |D|y² = 4p thì trace t = x (hoặc ±x tùy vào D)
        Order của curve = p + 1 - t
        """
        result = ECPPUtils.cornacchia(D, p)
        if result is None:
            return None

        x, y = result

        # Với D = -3, -4 có quy tắc đặc biệt
        if D == -3:
            # #E = p + 1 - t với t = x (nếu p ≡ 1 mod 3)
            if p % 3 == 1:
                return x
            elif p % 3 == 2:
                return -x

        if D == -4:
            # #E = p + 1 - t với t = x (nếu p ≡ 1 mod 4)
            if p % 4 == 1:
                return x
            elif p % 4 == 3:
                return 0  # No points with order p+1

        # General case
        return x

    @staticmethod
    def compute_curve_order_cm(D, p):
        """
        Tính order của elliptic curve bằng Complex Multiplication

        Đây là cách CHÍNH XÁC để tính order!

        Returns:
            order của curve, hoặc None nếu thất bại
        """
        trace = ECPPUtils.compute_trace_from_cornacchia(D, p)
        if trace is None:
            return None

        # Order = p + 1 - trace
        order = p + 1 - trace

        # Verify bằng Hasse bound
        hasse_lower = p + 1 - 2 * int(p ** 0.5)
        hasse_upper = p + 1 + 2 * int(p ** 0.5)

        if not (hasse_lower <= order <= hasse_upper):
            return None

        return order

    @staticmethod
    def find_all_discriminants_for_prime(p, max_D=200):
        """
        Tìm tất cả discriminants phù hợp với p
        (tức là các D mà Cornacchia có solution)
        """
        from ECPP_Core import ECPPHelper

        valid_discriminants = []

        for D in ECPPHelper.SMALL_DISCRIMINANTS:
            if abs(D) > max_D:
                break

            # Kiểm tra Kronecker symbol
            if ECPPHelper.kronecker_symbol(D, p) == 1:
                # Thử Cornacchia
                result = ECPPUtils.cornacchia(D, p)
                if result is not None:
                    order = ECPPUtils.compute_curve_order_cm(D, p)
                    if order is not None:
                        valid_discriminants.append((D, order))

        return valid_discriminants

    @staticmethod
    def factor_simple(n, max_tries=1000):
        """
        Phân tích n thành các thừa số (simplified)
        Trả về (small_factors, large_remainder)
        """
        from Prime_All import SMALL_PRIMES

        factors = {}
        remaining = n

        # Trial division với small primes
        for p in SMALL_PRIMES[:200]:
            if p * p > remaining:
                break
            while remaining % p == 0:
                factors[p] = factors.get(p, 0) + 1
                remaining //= p

        return factors, remaining

    @staticmethod
    def certificate_to_string(cert):
        """Format certificate thành string dễ đọc"""
        from ECPP_Types import format_certificate  # ✓ Import function từ types
        return format_certificate(cert, verbose=True)

    @staticmethod
    def save_certificates(certificates, filename):
        """Lưu certificates ra file"""
        import json

        data = []
        for cert in certificates:
            data.append(cert.to_dict())  # ✓ Dùng method từ dataclass

        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)

    @staticmethod
    def load_certificates(filename):
        """Load certificates từ file"""
        import json

        with open(filename, 'r') as f:
            data = json.load(f)

        certificates = []
        for item in data:
            cert = ECPPCertificate.from_dict(item)  # ✓ Dùng classmethod
            certificates.append(cert)

        return certificates


# Test functions
if __name__ == "__main__":
    print("=== Testing Tonelli-Shanks ===")
    test_cases = [
        (10, 13),  # 10 mod 13
        (5, 41),   # 5 mod 41
        (2, 7),    # 2 mod 7
    ]

    for n, p in test_cases:
        r = ECPPUtils.tonelli_shanks(n, p)
        if r is not None:
            print(f"√{n} ≡ {r} (mod {p})")
            print(f"  Verify: {r}² ≡ {modulo(r*r, p)} (mod {p})")
        else:
            print(f"√{n} does not exist (mod {p})")

    print("\n=== Testing Cornacchia ===")
    test_primes = [7, 13, 31, 97]

    for p in test_primes:
        print(f"\nFor p = {p}:")
        discriminants = ECPPUtils.find_all_discriminants_for_prime(p, max_D=50)

        if discriminants:
            for D, order in discriminants:
                print(f"  D = {D:4d} -> order = {order}")
        else:
            print("  No valid discriminants found")

    print("\n=== Testing Curve Order Computation ===")
    p = 101
    D = -3

    print(f"Computing order for p = {p}, D = {D}")
    result = ECPPUtils.cornacchia(D, p)
    if result:
        x, y = result
        print(f"  Cornacchia: x={x}, y={y}")
        print(f"  Verify: {x}² + {abs(D)}·{y}² = {x*x + abs(D)*y*y} = {4*p}")

        order = ECPPUtils.compute_curve_order_cm(D, p)
        print(f"  Curve order: {order}")
        print(f"  Hasse bound: [{p + 1 - 2*int(p**0.5)}, {p + 1 + 2*int(p**0.5)}]")