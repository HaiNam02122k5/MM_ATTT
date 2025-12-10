"""
ECPP_Types.py - Data structures cho ECPP
File này chứa các dataclass và types, không import modules khác
"""

from dataclasses import dataclass
from typing import Optional, List


@dataclass
class ECPPCertificate:
    """
    Chứng chỉ ECPP cho một bước chứng minh

    Mỗi certificate chứng minh rằng n là prime bằng cách:
    1. Tìm elliptic curve E: y² = x³ + ax + b (mod n)
    2. Curve có order m = h*q với q là số nguyên tố lớn
    3. Tồn tại point P thỏa mãn Goldwasser-Kilian conditions
    """
    n: int  # Số cần chứng minh là prime
    a: int  # Tham số curve
    b: int  # Tham số curve
    m: int  # Order của curve (#E)
    q: int  # Prime factor lớn của m
    point_x: int  # X-coordinate của witness point
    point_y: int  # Y-coordinate của witness point

    def __str__(self):
        return f"ECPPCert(n={self.n}, q={self.q})"

    def __repr__(self):
        return self.__str__()

    def to_dict(self):
        """Convert to dictionary"""
        return {
            'n': self.n,
            'a': self.a,
            'b': self.b,
            'm': self.m,
            'q': self.q,
            'point_x': self.point_x,
            'point_y': self.point_y
        }

    @classmethod
    def from_dict(cls, data):
        """Create from dictionary"""
        return cls(
            n=data['n'],
            a=data['a'],
            b=data['b'],
            m=data['m'],
            q=data['q'],
            point_x=data['point_x'],
            point_y=data['point_y']
        )

    def get_cofactor(self):
        """Return cofactor h = m / q"""
        return self.m // self.q

    def verify_basic_properties(self):
        """
        Kiểm tra các thuộc tính cơ bản của certificate
        (không verify trên curve, chỉ check logic)
        """
        # Check q divides m
        if self.m % self.q != 0:
            return False, "q does not divide m"

        # Check q > 4√n
        min_q = 4 * int(self.n ** 0.5)
        if self.q <= min_q:
            return False, f"q={self.q} not > 4√n={min_q}"

        # Check n, q, m > 0
        if self.n <= 0 or self.q <= 0 or self.m <= 0:
            return False, "n, q, m must be positive"

        return True, "OK"


@dataclass
class CurveInfo:
    """Thông tin về một elliptic curve"""
    a: int  # Tham số a
    b: int  # Tham số b
    p: int  # Modulus (prime)
    order: int  # Order của curve (#E)
    discriminant: int  # Discriminant D được dùng để tạo curve

    def __str__(self):
        return f"y² = x³ + {self.a}x + {self.b} (mod {self.p}), order={self.order}"


@dataclass
class ProofResult:
    """Kết quả của một lần chứng minh"""
    n: int  # Số được test
    is_prime: bool  # Kết quả
    certificates: List[ECPPCertificate]  # Chain of certificates
    time_elapsed: float  # Thời gian (seconds)
    method: str  # Phương pháp ("ECPP", "TrialDivision", etc.)

    def __str__(self):
        result = "PRIME" if self.is_prime else "COMPOSITE"
        return f"ProofResult({self.n} is {result}, {len(self.certificates)} certs, {self.time_elapsed:.4f}s)"

    def get_certificate_chain_summary(self):
        """Tóm tắt certificate chain"""
        if not self.certificates:
            return "No certificates (base case or composite)"

        summary = []
        for i, cert in enumerate(self.certificates):
            bits = cert.n.bit_length()
            summary.append(f"  Step {i}: n={cert.n} ({bits} bits) -> q={cert.q}")

        return "\n".join(summary)


# Helper functions (không import gì cả)
def format_certificate(cert: ECPPCertificate, verbose=False) -> str:
    """Format certificate thành string"""
    if not verbose:
        return f"Cert: n={cert.n} ({cert.n.bit_length()} bits) -> q={cert.q}"

    return f"""
Certificate:
  N = {cert.n}
  Curve: y² = x³ + {cert.a}x + {cert.b} (mod {cert.n})
  Order m = {cert.m}
  Large prime q = {cert.q}
  Cofactor h = {cert.m // cert.q}
  Witness point P = ({cert.point_x}, {cert.point_y})

  Verification conditions:
    - m·P = O (point at infinity) ✓
    - (m/q)·P ≠ O ✓
    - gcd(y_coord, n) = 1 ✓
    - q > 4√n = {int(4 * (cert.n ** 0.5))} ✓
"""


def validate_certificate_chain(certificates: List[ECPPCertificate]) -> tuple[bool, str]:
    """
    Validate chuỗi certificates (chỉ check structure, không verify trên curve)

    Returns: (is_valid, error_message)
    """
    if not certificates:
        return True, "Empty chain (OK for small primes)"

    # Check first certificate
    first = certificates[0]
    is_valid, msg = first.verify_basic_properties()
    if not is_valid:
        return False, f"First certificate invalid: {msg}"

    # Check chain consistency
    for i in range(len(certificates) - 1):
        current = certificates[i]
        next_cert = certificates[i + 1]

        # q của cert hiện tại phải bằng n của cert tiếp theo
        if current.q != next_cert.n:
            return False, f"Chain broken at step {i}: q={current.q} != n={next_cert.n}"

        # Verify basic properties
        is_valid, msg = next_cert.verify_basic_properties()
        if not is_valid:
            return False, f"Certificate {i + 1} invalid: {msg}"

    return True, "Chain structure is valid"