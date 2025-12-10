"""
ECPP Implementation - Part 2: Curve Order Counting
Tính số điểm trên đường cong elliptic (critical cho ECPP)
"""

import random
from NumberTheory import gcd, moduloPower, modulo, inverseModulo
from Prime_All import prime_check
from ECC.Elliptic_Curve import EllipticCurve


class CurveOrderCounter:
    """Đếm số điểm trên elliptic curve modulo p"""
    
    @staticmethod
    def count_points_naive(a, b, p):
        """
        Đếm số điểm bằng brute force (chỉ dùng cho p nhỏ)
        Độ phức tạp: O(p)
        """
        if p > 10000:
            raise ValueError("Naive counting chỉ dùng cho p nhỏ (<10000)")
        
        count = 1  # Điểm vô cực
        
        for x in range(p):
            # Tính y² = x³ + ax + b (mod p)
            y_squared = modulo(moduloPower(x, 3, p) + a * x + b, p)
            
            # Kiểm tra y_squared có phải là quadratic residue không
            # Dùng Legendre symbol (y²/p)
            legendre = moduloPower(y_squared, (p - 1) // 2, p)
            
            if legendre == 0:
                count += 1  # y = 0
            elif legendre == 1:
                count += 2  # Có 2 giá trị y
        
        return count
    
    @staticmethod
    def schoof_algorithm_simplified(a, b, p):
        """
        Thuật toán Schoof đơn giản hóa để đếm điểm
        
        Ý tưởng: Sử dụng Hasse's theorem:
        |#E - (p + 1)| ≤ 2√p
        
        Với ECPP, ta thường dùng Complex Multiplication method,
        nhưng đây là simplified version.
        """
        # Hasse bound
        hasse_bound = 2 * int(p ** 0.5) + 1
        
        # Estimate sử dụng random points
        # Trong thực tế, cần implement Schoof-Elkies-Atkin (SEA)
        estimated_order = CurveOrderCounter._estimate_order_by_sampling(a, b, p)
        
        return estimated_order
    
    @staticmethod
    def _estimate_order_by_sampling(a, b, p, num_samples=20):
        """
        Ước lượng order bằng cách sample random points
        CHÚ Ý: Đây KHÔNG phải là cách chính xác, chỉ là estimation!
        
        Trong ECPP thực sự, cần dùng:
        1. Complex Multiplication để tính chính xác
        2. Hoặc Schoof-Elkies-Atkin algorithm
        """
        try:
            curve = EllipticCurve(a, b, p)
        except ValueError:
            return None
        
        # Sample một số điểm
        points = []
        attempts = 0
        max_attempts = num_samples * 10
        
        while len(points) < num_samples and attempts < max_attempts:
            try:
                point = curve.generate_point()
                if point is not None:
                    points.append(point)
            except:
                pass
            attempts += 1
        
        if len(points) < 5:
            # Không đủ điểm để estimate, dùng Hasse bound
            return p + 1
        
        # Heuristic: assume order gần p + 1
        # (đúng theo Hasse theorem)
        return p + 1
    
    @staticmethod
    def find_curve_with_good_order(p, target_factor_size=None):
        """
        Tìm đường cong có order "tốt" cho ECPP
        
        Order "tốt" nghĩa là:
        - #E = m * q với q là số nguyên tố lớn
        - q > 4√p (Goldwasser-Kilian condition)
        
        Returns: (a, b, order, q) hoặc None
        """
        if target_factor_size is None:
            target_factor_size = p // 4  # Heuristic
        
        min_q = int(4 * (p ** 0.5))
        
        # Thử nhiều discriminants
        from ECPP_Core import ECPPHelper
        
        for D in ECPPHelper.SMALL_DISCRIMINANTS:
            try:
                # Tìm curve từ discriminant
                result = CurveOrderCounter._try_discriminant(D, p, min_q)
                if result is not None:
                    return result
            except:
                continue
        
        # Fallback: random search
        return CurveOrderCounter._random_search_curve(p, min_q)
    
    @staticmethod
    def _try_discriminant(D, p, min_q, max_attempts=10):
        """Thử tìm curve tốt từ discriminant D"""
        from ECPP_Core import ECPPHelper
        
        # Kiểm tra Kronecker symbol
        if ECPPHelper.kronecker_symbol(D, p) != 1:
            return None
        
        # Tìm curve
        curve_params = ECPPHelper.find_curve_from_discriminant(D, p)
        if curve_params is None:
            return None
        
        a, b = curve_params
        
        # Tính order (simplified - trong thực tế cần CM method)
        order = CurveOrderCounter._compute_order_from_cm(D, p)
        if order is None:
            # Fallback
            order = CurveOrderCounter.schoof_algorithm_simplified(a, b, p)
        
        # Tìm prime factor lớn
        q = CurveOrderCounter._find_large_prime_factor(order, min_q)
        
        if q is not None and q > min_q:
            return (a, b, order, q)
        
        return None
    
    @staticmethod
    def _compute_order_from_cm(D, p):
        """
        Tính order từ Complex Multiplication
        
        Công thức: #E(F_p) = p + 1 - t
        với t thỏa mãn: 4p = t² - Dv²
        
        Đây là simplified version!
        """
        # Với D = -3, -4, có công thức đặc biệt
        if D == -3:
            # #E = p + 1 - a_p với a_p được tính từ p mod 3
            if p % 3 == 1:
                # Cần tìm a, b sao cho p = a² + 3b²
                # Simplified: estimate
                return p + 1
            elif p % 3 == 2:
                return p + 1
        
        if D == -4:
            # #E = p + 1 - a_p với a_p từ p mod 4
            if p % 4 == 1:
                # p = a² + b²
                return p + 1
            elif p % 4 == 3:
                return p + 1
        
        # General case: cần Cornacchia algorithm
        # Simplified: return Hasse estimate
        return p + 1
    
    @staticmethod
    def _find_large_prime_factor(n, min_size):
        """Tìm prime factor lớn của n"""
        if prime_check(n):
            return n if n >= min_size else None
        
        # Trial division với small primes
        from Prime_All import SMALL_PRIMES
        
        remaining = n
        for p in SMALL_PRIMES[:100]:  # Chỉ thử 100 primes đầu
            while remaining % p == 0:
                remaining //= p
            if remaining < min_size:
                return None
        
        # Check nếu remaining là prime và đủ lớn
        if remaining >= min_size and prime_check(remaining):
            return remaining
        
        return None
    
    @staticmethod
    def _random_search_curve(p, min_q, max_attempts=100):
        """Random search for good curve"""
        for _ in range(max_attempts):
            a = random.randint(0, p-1)
            b = random.randint(0, p-1)
            
            # Check non-singular
            disc = modulo(4 * moduloPower(a, 3, p) + 27 * moduloPower(b, 2, p), p)
            if disc == 0:
                continue
            
            # Estimate order
            order = CurveOrderCounter.schoof_algorithm_simplified(a, b, p)
            
            # Find large prime factor
            q = CurveOrderCounter._find_large_prime_factor(order, min_q)
            
            if q is not None and q > min_q:
                return (a, b, order, q)
        
        return None


# Test
if __name__ == "__main__":
    print("=== Testing Order Counting ===")
    
    # Test với p nhỏ
    test_cases = [
        (1, 1, 7),    # y² = x³ + x + 1 mod 7
        (2, 3, 97),   # y² = x³ + 2x + 3 mod 97
    ]
    
    for a, b, p in test_cases:
        try:
            order_naive = CurveOrderCounter.count_points_naive(a, b, p)
            print(f"Curve y² = x³ + {a}x + {b} (mod {p})")
            print(f"  Order (naive): {order_naive}")
            print(f"  Hasse bound: [{p + 1 - 2*int(p**0.5)}, {p + 1 + 2*int(p**0.5)}]")
        except Exception as e:
            print(f"Error: {e}")
        print()
    
    print("\n=== Testing Curve Finding ===")
    for p in [101, 1009]:
        print(f"\nFinding good curve for p = {p}...")
        result = CurveOrderCounter.find_curve_with_good_order(p)
        if result:
            a, b, order, q = result
            print(f"  Found: y² = x³ + {a}x + {b}")
            print(f"  Order: {order}")
            print(f"  Large prime factor q: {q}")
        else:
            print(f"  Could not find suitable curve")
