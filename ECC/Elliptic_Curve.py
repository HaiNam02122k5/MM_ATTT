import hashlib
import os
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
    """Lớp đại diện cho đường cong Curve25519/Ed25519"""

    def __init__(self, curve_type="ed25519"):
        """
        Khởi tạo đường cong
        curve_type: "curve25519" (Montgomery) hoặc "ed25519" (Edwards)
        """
        self.curve_type = curve_type.lower()

        # Tham số chung
        self.p = 2 ** 255 - 19  # Prime field modulus
        self.P = self.p  # Giữ P để tương thích với code cũ
        self.n = 2 ** 252 + 27742317777372353535851937790883648493  # Order (cấp của đường cong)
        self.L = self.n  # Giữ L để tương thích với code cũ

        if self.curve_type == "ed25519":
            # Ed25519 parameters (Edwards curve - dùng cho signature)
            # Phương trình: -x^2 + y^2 = 1 + d*x^2*y^2
            self.a = -1  # Hệ số a trong twisted Edwards curve
            self.d = -121665 * int(inverseModulo(121666, self.p)) % self.p  # Hệ số d
            self.D = self.d  # Giữ D để tương thích

            # Base point (điểm sinh)
            self.Bx = 15112221349535400772501151409588531511454012693041857206046113283949847762202
            self.By = 46316835694926478169428394003475163141307993866256225615783033603165251855960
            self.B = (self.Bx, self.By)
            self.G = self.B  # Generator point

            # Không có b cho Edwards curve (dùng a và d)
            self.b = None
        else:
            # Curve25519 parameters (Montgomery curve - dùng cho ECDH)
            # Phương trình: By^2 = x^3 + Ax^2 + x
            self.A = 486662  # Hệ số A
            self.a = self.A
            self.B_coeff = 1  # Hệ số B trong phương trình Montgomery
            self.b = self.B_coeff

            # Base point u-coordinate (điểm sinh)
            self.Bu = 9
            # v-coordinate tương ứng (tính được từ phương trình)
            self.Bv = self._compute_v_from_u(self.Bu)
            self.G = self.Bu  # Generator u-coordinate
            self.B = (self.Bu, self.Bv)  # Base point đầy đủ

    def _compute_v_from_u(self, u):
        """Tính v-coordinate từ u-coordinate cho Montgomery curve"""
        # By^2 = x^3 + Ax^2 + x
        # v^2 = u^3 + A*u^2 + u (mod p)
        rhs = (pow(u, 3, self.p) + self.A * pow(u, 2, self.p) + u) % self.p
        # Tính căn bậc hai modulo p (Tonelli-Shanks có thể dùng)
        # Với p ≡ 3 (mod 4), có thể dùng v = rhs^((p+1)/4) mod p
        # Nhưng p = 2^255 - 19 ≡ 5 (mod 8), cần thuật toán phức tạp hơn
        # Để đơn giản, trả về None hoặc giá trị đã biết
        # v cho u=9 được biết trước
        if u == 9:
            return 14781619447589544791020593568409986887264606134616475288964881837755586237401
        return None

    def get_curve_params(self):
        """Trả về dictionary chứa các tham số của đường cong"""
        params = {
            'type': self.curve_type,
            'p': self.p if hasattr(self, 'p') else self.P,
            'n': self.L,
        }

        if self.curve_type == "ed25519":
            params.update({
                'equation': '-x^2 + y^2 = 1 + d*x^2*y^2',
                'a': -1,
                'd': self.D,
                'base_point': f"G = ({self.Bx}, {self.By})"
            })
        else:
            params.update({
                'equation': 'By^2 = x^3 + Ax^2 + x',
                'A': self.A,
                'B': 1,
                'base_point': f"G = (u={self.Bu})"
            })

        return params

    def print_curve_info(self):
        """In thông tin đầy đủ về đường cong"""
        p = self.p if hasattr(self, 'p') else self.P

        print(f"\n{'=' * 70}")
        print(f"Curve Type: {self.curve_type.upper()}")
        print(f"{'=' * 70}")

        if self.curve_type == "ed25519":
            print(f"Equation: -x² + y² = 1 + d·x²·y²")
            print(f"\nParameters:")
            print(f"  p (prime modulus) = 2^255 - 19")
            print(f"                    = {p}")
            print(f"  a = -1")
            print(f"  d = {self.D}")
            print(f"  n (order)        = {self.L}")
            print(f"\nBase Point (Generator G):")
            print(f"  Gx = {self.Bx}")
            print(f"  Gy = {self.By}")
        else:
            print(f"Equation: B·y² = x³ + A·x² + x")
            print(f"\nParameters:")
            print(f"  p (prime modulus) = 2^255 - 19")
            print(f"                    = {p}")
            print(f"  A = {self.A}")
            print(f"  B = 1")
            print(f"  n (order)        = {self.L}")
            print(f"\nBase Point (Generator G):")
            print(f"  Gu = {self.Bu}")

        print(f"{'=' * 70}\n")

    def point_add(self, P1, P2):
        """Cộng hai điểm trên Ed25519"""
        if self.curve_type != "ed25519":
            raise NotImplementedError("Point addition chỉ implement cho Ed25519")

        if P1 is None:
            return P2
        if P2 is None:
            return P1

        x1, y1 = P1
        x2, y2 = P2

        x1x2 = x1 * x2
        y1y2 = y1 * y2
        dx1x2y1y2 = self.D * x1x2 * y1y2

        x3 = ((x1 * y2 + y1 * x2) * int(inverseModulo(1 + dx1x2y1y2, self.P))) % self.P
        y3 = ((y1y2 + x1x2) * int(inverseModulo(1 - dx1x2y1y2, self.P))) % self.P

        return (x3, y3)

    def scalar_mult(self, k, point=None):
        """Nhân vô hướng điểm"""
        if self.curve_type == "ed25519":
            return self._ed25519_scalar_mult(k, point)
        else:
            return self._curve25519_scalar_mult(k, point)

    def _ed25519_scalar_mult(self, k, point=None):
        """Nhân vô hướng cho Ed25519"""
        if point is None:
            point = self.B

        if k == 0:
            return None
        if k == 1:
            return point

        result = None
        addend = point

        while k:
            if k & 1:
                result = self.point_add(result, addend)
            addend = self.point_add(addend, addend)
            k >>= 1

        return result

    def _curve25519_scalar_mult(self, k, u=None):
        """
        Nhân vô hướng cho Curve25519 (Montgomery ladder)
        Trả về u-coordinate của k*P
        """
        if u is None:
            u = self.Bu

        # Montgomery ladder để tránh timing attacks
        u2, w2 = (1, 0)
        u3, w3 = (u, 1)

        k_bits = bin(k)[2:]

        for bit in k_bits:
            if bit == '1':
                u2, w2, u3, w3 = self._differential_add(u2, w2, u3, w3, u)
            else:
                u3, w3, u2, w2 = self._differential_add(u3, w3, u2, w2, u)

        return (u2 * int(inverseModulo(w2, self.P))) % self.P

    def _differential_add(self, u2, w2, u3, w3, u):
        """Phép cộng differential cho Montgomery ladder"""
        A = (u2 + w2) % self.P
        AA = (A * A) % self.P
        B = (u2 - w2) % self.P
        BB = (B * B) % self.P
        E = (AA - BB) % self.P
        C = (u3 + w3) % self.P
        D = (u3 - w3) % self.P
        DA = (D * A) % self.P
        CB = (C * B) % self.P

        u_sum = pow(DA + CB, 2, self.P)
        u_diff = (u * pow(DA - CB, 2, self.P)) % self.P

        u2_new = (AA * BB) % self.P
        w2_new = (E * (AA + ((self.A - 2) // 4) * E)) % self.P

        return u2_new, w2_new, u_sum, u_diff

    def generate_private_key(self):
        """Tạo private key ngẫu nhiên"""
        private_key = int.from_bytes(os.urandom(32), 'little')

        if self.curve_type == "ed25519":
            private_key = private_key % self.L
            if private_key == 0:
                private_key = 1
        else:
            # Clamping cho Curve25519
            private_key &= ~7  # Clear 3 bits thấp nhất
            private_key &= ~(128 << 8 * 31)  # Clear bit cao nhất
            private_key |= 64 << 8 * 31  # Set bit thứ 2 cao nhất

        return private_key

    def get_public_key(self, private_key):
        """Tính public key từ private key"""
        return self.scalar_mult(private_key)

    def derive_signing_key(self, shared_secret, context=b"signature"):
        """
        Derive signing key từ shared secret
        shared_secret: shared secret từ ECDH
        context: context string để phân biệt các mục đích sử dụng
        """
        if self.curve_type != "ed25519":
            # Chuyển sang Ed25519 để signing
            ed_curve = Curve25519("ed25519")
            return ed_curve.derive_signing_key(shared_secret, context)

        # Hash shared secret với context để tạo signing key
        if isinstance(shared_secret, int):
            shared_secret = shared_secret.to_bytes(32, 'little')

        material = context + shared_secret
        derived = hashlib.sha512(material).digest()

        # Lấy 32 bytes đầu làm private key
        signing_key = int.from_bytes(derived[:32], 'little') % self.L
        if signing_key == 0:
            signing_key = 1

        return signing_key

    def __str__(self):
        return f"Curve25519({self.curve_type})"

    def __repr__(self):
        return self.__str__()