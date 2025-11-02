import random

from NumberTheory import inverseModulo, moduloPower, modulo


class EllipticCurve:
    """Đường cong Elliptic dạng Weierstrass: y^2 = x^3 + ax + b (mod p)"""

    def __init__(self, a, b, p):
        """
        Khởi tạo đường cong Elliptic
        Args:
            a, b: Tham số đường cong
            p: Số nguyên tố (trường hữu hạn)
        """
        self.a = a
        self.b = b
        self.p = p

        # Kiểm tra điều kiện không suy biến: 4a^3 + 27b^2 != 0 (mod p)
        if modulo(4 * moduloPower(a, 3, p) + 27 * moduloPower(b, 2, p), p) == 0:
            raise ValueError("Đường cong suy biến! 4a^3 + 27b^2 ≡ 0 (mod p)")

    def is_on_curve(self, point):
        """Kiểm tra điểm có nằm trên đường cong không"""
        if point is None:  # Điểm vô cực
            return True

        x, y = point
        return modulo(moduloPower(y, 2, self.p) - (moduloPower(x, 3, self.p) + self.a * x + self.b), self.p) == 0

    def point_add(self, P, Q):
        """
        Phép cộng hai điểm trên đường cong Elliptic
        Returns: P + Q
        """
        # Trường hợp đặc biệt: điểm vô cực (O)
        if P is None:
            return Q
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q

        # Trường hợp P = -Q (cùng x, y đối nhau)
        if x1 == x2 and modulo(y1 + y2, self.p) == 0:
            return None  # Kết quả là điểm vô cực

        # Tính hệ số góc lambda
        if P == Q:  # Phép nhân đôi điểm
            # lambda = (3x1^2 + a) / (2y1)
            numerator = modulo(3 * moduloPower(x1, 2, self.p) + self.a, self.p)
            denominator = (2 * y1) % self.p
        else:  # Phép cộng hai điểm khác nhau
            # lambda = (y2 - y1) / (x2 - x1)
            numerator = modulo(y2 - y1, self.p)
            denominator = modulo(x2 - x1, self.p)

        # Tính lambda = numerator * denominator^(-1) mod p
        denominator_inv = inverseModulo(denominator, self.p)
        lambda_val = modulo(numerator * denominator_inv, self.p)

        # Tính tọa độ điểm mới
        x3 = modulo(moduloPower(lambda_val, 2, self.p) - x1 - x2, self.p)
        y3 = modulo(lambda_val * (x1 - x3) - y1,  self.p)

        return (x3, y3)

    def point_multiply(self, k, P):
        """
        Phép nhân vô hướng: k * P (cộng P với chính nó k lần)
        Sử dụng thuật toán "double and add" để tối ưu
        """
        if k == 0:
            return None  # Điểm vô cực

        if k < 0:
            # k * P = -k * (-P)
            k = -k
            P = (P[0], (-P[1]) % self.p)

        result = None  # Điểm vô cực
        addend = P

        # Thuật toán double-and-add
        while k:
            if k & 1:  # Bit thấp nhất là 1
                result = self.point_add(result, addend)
            addend = self.point_add(addend, addend)  # Nhân đôi
            k >>= 1  # Dịch phải 1 bit

        return result

    def generate_point(self):
        """Tìm một điểm ngẫu nhiên trên đường cong"""
        for _ in range(1000):  # Thử tối đa 1000 lần
            x = random.randint(0, self.p - 1)
            # Tính y^2 = x^3 + ax + b (mod p)
            y_squared = modulo(moduloPower(x, 3, self.p) + self.a * x + self.b, self.p)

            # Kiểm tra xem y_squared có phải là số chính phương không
            y = self._mod_sqrt(y_squared, self.p)
            if y is not None:
                return (x, y)

        raise ValueError("Không tìm thấy điểm trên đường cong sau 1000 lần thử!")

    def _mod_sqrt(self, a, p):
        """Tính căn bậc hai modulo p (thuật toán Tonelli-Shanks đơn giản)"""
        # Chỉ làm việc với p ≡ 3 (mod 4) để đơn giản
        if p % 4 == 3:
            y = moduloPower(a, (p + 1) // 4, p)
            if moduloPower(y, 2, p) == a:
                return y

        # Thử brute force cho số nhỏ
        if p < 1000000:
            for y in range(p):
                if moduloPower(y, 2, p) == a:
                    return y

        return None

    def __str__(self):
        return f"EllipticCurve: y^2 = x^3 + {self.a}x + {self.b} (mod {self.p})"


class Curve25519:
    """
    Curve25519: Đường cong Montgomery
    Dạng: By² = x³ + Ax² + x (mod p)

    Tham số chuẩn:
    - p = 2^255 - 19
    - A = 486662
    - B = 1
    - Base point: x = 9
    - Order (bậc): n = 2^252 + 27742317777372353535851937790883648493
    """

    def __init__(self):
        # Tham số đường cong Curve25519
        self.p = pow(2, 255) - 19
        self.A = 486662
        self.B = 1

        # Điểm sinh G (base point)
        # Curve25519 thường dùng x-coordinate only
        self.G_x = 9

        # Tính y từ x = 9
        # By² = x³ + Ax² + x (mod p)
        # y² = (x³ + Ax² + x) / B (mod p)
        rhs = (moduloPower(self.G_x, 3, self.p) + self.A * moduloPower(self.G_x, 2, self.p) + self.G_x) % self.p
        rhs = (rhs * inverseModulo(self.B, self.p)) % self.p

        # Tính căn bậc hai (với p ≡ 3 mod 4, dùng y = rhs^((p+1)/4))
        # Nhưng p = 2^255 - 19 ≡ 5 (mod 8), cần phương pháp khác
        self.G_y = self._mod_sqrt_p_5_mod_8(rhs)

        self.G = (self.G_x, self.G_y)

        # Bậc của đường cong (order/cofactor)
        # Curve25519: #E = 8 * n (cofactor h = 8)
        self.n = pow(2, 252) + 27742317777372353535851937790883648493
        self.cofactor = 8

        # Bậc thực tế của điểm sinh G
        self.order = self.n

        print("=" * 70)
        print("CURVE25519 - Thông số chuẩn")
        print("=" * 70)
        print(f"Đường cong Montgomery: By² = x³ + Ax² + x (mod p)")
        print(f"\nTham số:")
        print(f"  p = 2^255 - 19")
        print(f"    = {self.p}")
        print(f"  A = {self.A}")
        print(f"  B = {self.B}")
        print(f"\nĐiểm sinh G:")
        print(f"  G.x = {self.G_x}")
        print(f"  G.y = {self.G_y}")
        print(f"\nBậc (Order):")
        print(f"  n = 2^252 + 27742317777372353535851937790883648493")
        print(f"    = {self.n}")
        print(f"  Cofactor = {self.cofactor}")
        print(f"  #E(Fp) = {self.cofactor * self.n}")
        print("=" * 70)

    def _mod_sqrt_p_5_mod_8(self, a):
        """
        Tính căn bậc hai modulo p khi p ≡ 5 (mod 8)
        Thuật toán Tonelli-Shanks cho trường hợp đặc biệt
        """
        p = self.p

        # Công thức: y = a^((p+3)/8) hoặc y = a^((p+3)/8) * 2^((p-1)/4)
        exp1 = (p + 3) // 8
        y = moduloPower(a, exp1, p)

        # Kiểm tra
        if moduloPower(y, 2, p) == a:
            return y

        # Thử phương án 2
        exp2 = (p - 1) // 4
        y = (y * moduloPower(2, exp2, p)) % p

        if moduloPower(y, 2, p) == a:
            return y

        raise ValueError("Không tìm thấy căn bậc hai!")

    def is_on_curve(self, point):
        """Kiểm tra điểm có nằm trên đường cong không"""
        if point is None:
            return True

        x, y = point
        # By² = x³ + Ax² + x (mod p)
        lhs = (self.B * moduloPower(y, 2, self.p)) % self.p
        rhs = (moduloPower(x, 3, self.p) + self.A * moduloPower(x, 2, self.p) + x) % self.p

        return lhs == rhs

    def point_add(self, P, Q):
        """
        Phép cộng điểm trên đường cong Montgomery
        Công thức khác với Weierstrass!
        """
        if P is None:
            return Q
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q

        # Trường hợp P = -Q
        if x1 == x2 and (y1 + y2) % self.p == 0:
            return None

        # Tính hệ số góc
        if P == Q:
            # Phép nhân đôi
            # λ = (3x1² + 2Ax1 + 1) / (2By1)
            numerator = (3 * moduloPower(x1, 2, self.p) + 2 * self.A * x1 + 1) % self.p
            denominator = (2 * self.B * y1) % self.p
        else:
            # Phép cộng
            # λ = (y2 - y1) / (x2 - x1)
            numerator = (y2 - y1) % self.p
            denominator = (x2 - x1) % self.p

        lambda_val = (numerator * inverseModulo(denominator, self.p)) % self.p

        # Tính x3, y3
        # x3 = Bλ² - A - x1 - x2
        x3 = (self.B * moduloPower(lambda_val, 2, self.p) - self.A - x1 - x2) % self.p

        # y3 = λ(x1 - x3) - y1
        y3 = (lambda_val * (x1 - x3) - y1) % self.p

        return (x3, y3)

    def point_multiply(self, k, P):
        """
        Phép nhân vô hướng: k * P
        Sử dụng Montgomery ladder (an toàn với timing attack)
        """
        if k == 0:
            return None

        if k < 0:
            k = -k
            P = (P[0], (-P[1]) % self.p)

        # Montgomery ladder
        R0 = None  # O
        R1 = P

        # Lấy bit của k từ trái sang phải
        bits = bin(k)[2:]  # Bỏ '0b'

        for bit in bits:
            if bit == '0':
                R1 = self.point_add(R0, R1)
                R0 = self.point_add(R0, R0)
            else:
                R0 = self.point_add(R0, R1)
                R1 = self.point_add(R1, R1)

        return R0
