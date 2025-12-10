"""
ECPP - Elliptic Curve Primality Proving
Main implementation của thuật toán Goldwasser-Kilian-Atkin

CHÚ Ý: Đây là educational implementation.
Production version cần:
- Hilbert class polynomial computation (dùng PARI/GP hoặc SageMath)
- Cornacchia algorithm cho CM
- Schoof-Elkies-Atkin cho order counting
"""

import time
from typing import Optional, List, Tuple

from NumberTheory import gcd, moduloPower, modulo
from Prime_All import prime_check, SMALL_PRIMES, generate_prime_bit
from ECC.Elliptic_Curve import EllipticCurve
from ECPP_Types import ECPPCertificate  # ✓ Import từ file riêng


class ECPP:
    """
    Elliptic Curve Primality Proving

    Thuật toán:
    1. Nếu n nhỏ, dùng trial division hoặc Miller-Rabin
    2. Tìm elliptic curve E với order m = h*q (q prime lớn)
    3. Tìm point P sao cho:
       - m*P = O (điểm vô cực)
       - (m/q)*P ≠ O
    4. Đệ quy chứng minh q là prime
    5. Xây dựng certificate chain
    """

    # Threshold để dừng đệ quy
    SMALL_PRIME_THRESHOLD = 10 ** 6

    def __init__(self, verbose=True):
        self.verbose = verbose
        self.certificates = []

    def prove_prime(self, n) -> Tuple[bool, List[ECPPCertificate]]:
        """
        Chứng minh n là số nguyên tố

        Returns:
            (is_prime, certificates)
        """
        self.certificates = []

        # Step 1: Kiểm tra trivial cases
        if n < 2:
            return False, []

        if n == 2:
            return True, []

        if n % 2 == 0:
            return False, []

        # Step 2: Kiểm tra với small primes
        for p in SMALL_PRIMES[:1000]:
            if n == p:
                return True, []
            if n % p == 0:
                return False, []

        # Step 3: Nếu n nhỏ, dùng deterministic test
        if n < self.SMALL_PRIME_THRESHOLD:
            is_prime = prime_check(n)
            return is_prime, []

        # Step 4: Chạy ECPP đệ quy
        self._log(f"\n{'='*60}")
        self._log(f"Starting ECPP for n = {n}")
        self._log(f"{'='*60}")

        try:
            result = self._ecpp_recursive(n, depth=0)
            return result, self.certificates
        except Exception as e:
            self._log(f"ECPP failed: {e}")
            return False, []

    def _ecpp_recursive(self, n, depth=0) -> bool:
        """
        Recursive ECPP

        Goldwasser-Kilian criterion:
        Nếu tồn tại elliptic curve E với order m = h*q, q prime, q > 4√n,
        và point P thỏa mãn:
        1. m*P = O
        2. (m/q)*P ≠ O
        3. gcd(m/q * P_y, n) = 1
        Thì n là prime
        """
        indent = "  " * depth
        self._log(f"{indent}Proving n = {n} (depth {depth})")

        # Base case: n đủ nhỏ
        if n < self.SMALL_PRIME_THRESHOLD:
            is_prime = prime_check(n)
            self._log(f"{indent}  Base case: {'PRIME' if is_prime else 'COMPOSITE'}")
            return is_prime

        # Step 1: Tìm suitable curve
        self._log(f"{indent}  Finding suitable curve...")
        curve_data = self._find_suitable_curve(n, depth)

        if curve_data is None:
            self._log(f"{indent}  ✗ Could not find suitable curve")
            # Fallback to probabilistic test
            return prime_check(n)

        a, b, m, q = curve_data
        self._log(f"{indent}  ✓ Found curve: y² = x³ + {a}x + {b}")
        self._log(f"{indent}    Order m = {m}, q = {q}")

        # Step 2: Tìm witness point
        self._log(f"{indent}  Finding witness point...")
        point = self._find_witness_point(a, b, m, q, n)

        if point is None:
            self._log(f"{indent}  ✗ Could not find witness point")
            return False

        self._log(f"{indent}  ✓ Found witness P = {point}")

        # Step 3: Verify Goldwasser-Kilian conditions
        self._log(f"{indent}  Verifying GK conditions...")
        if not self._verify_gk_conditions(a, b, m, q, point, n):
            self._log(f"{indent}  ✗ GK conditions failed")
            return False

        self._log(f"{indent}  ✓ GK conditions passed")

        # Step 4: Tạo certificate
        cert = ECPPCertificate(
            n=n, a=a, b=b, m=m, q=q,
            point_x=point[0], point_y=point[1]
        )
        self.certificates.append(cert)

        # Step 5: Đệ quy chứng minh q
        self._log(f"{indent}  Recursing to prove q = {q}")
        return self._ecpp_recursive(q, depth + 1)

    def _find_suitable_curve(self, n, depth=0) -> Optional[Tuple[int, int, int, int]]:
        """
        Tìm curve E với order m = h*q thỏa mãn:
        - q là prime
        - q > 4√n
        """
        from ECPP_Part2 import CurveOrderCounter

        min_q = int(4 * (n ** 0.5))

        # Thử tìm với các discriminants
        result = CurveOrderCounter.find_curve_with_good_order(n, min_q)

        if result is not None:
            a, b, m, q = result
            # Verify q > 4√n
            if q > min_q:
                return (a, b, m, q)

        return None

    def _find_witness_point(self, a, b, m, q, p) -> Optional[Tuple[int, int]]:
        """
        Tìm point P thỏa mãn:
        1. m*P = O
        2. (m/q)*P ≠ O
        """
        try:
            curve = EllipticCurve(a, b, p)
        except ValueError:
            return None

        h = m // q

        # Thử nhiều random points
        for attempt in range(100):
            try:
                # Generate random point
                P = curve.generate_point()
                if P is None:
                    continue

                # Check m*P = O
                mP = curve.point_multiply(m, P)
                if mP is not None:  # Phải là điểm vô cực
                    continue

                # Check (m/q)*P ≠ O
                hP = curve.point_multiply(h, P)
                if hP is None:  # Không được là điểm vô cực
                    continue

                # Check gcd(y_coord, n) = 1
                if gcd(hP[1], p) != 1:
                    continue

                return P

            except Exception as e:
                continue

        return None

    def _verify_gk_conditions(self, a, b, m, q, point, n) -> bool:
        """Verify Goldwasser-Kilian conditions"""
        try:
            curve = EllipticCurve(a, b, n)

            # 1. Check m*P = O
            mP = curve.point_multiply(m, point)
            if mP is not None:
                return False

            # 2. Check (m/q)*P ≠ O
            h = m // q
            hP = curve.point_multiply(h, point)
            if hP is None:
                return False

            # 3. Check gcd(y, n) = 1
            if gcd(hP[1], n) != 1:
                return False

            # 4. Check q > 4√n
            if q <= 4 * int(n ** 0.5):
                return False

            return True

        except Exception as e:
            return False

    def _log(self, message):
        """Log message if verbose"""
        if self.verbose:
            print(message)

    def verify_certificate_chain(self, n, certificates) -> bool:
        """
        Verify chuỗi chứng chỉ ECPP

        Returns True nếu certificates hợp lệ
        """
        if not certificates:
            return prime_check(n)

        # Verify từng certificate
        for i, cert in enumerate(certificates):
            if not self._verify_single_certificate(cert):
                print(f"Certificate {i} failed verification")
                return False

        # Verify chain consistency
        if certificates[0].n != n:
            return False

        for i in range(len(certificates) - 1):
            if certificates[i].q != certificates[i+1].n:
                return False

        # Verify base case
        last_q = certificates[-1].q
        if last_q >= self.SMALL_PRIME_THRESHOLD:
            return False

        return prime_check(last_q)

    def _verify_single_certificate(self, cert: ECPPCertificate) -> bool:
        """Verify một certificate"""
        try:
            curve = EllipticCurve(cert.a, cert.b, cert.n)
            point = (cert.point_x, cert.point_y)

            # Verify GK conditions
            return self._verify_gk_conditions(
                cert.a, cert.b, cert.m, cert.q, point, cert.n
            )
        except:
            return False


# Main test
if __name__ == "__main__":
    ecpp = ECPP(verbose=True)

    test_n = generate_prime_bit(100)

    print(ecpp.prove_prime(test_n))
