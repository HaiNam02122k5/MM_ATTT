import hashlib
import multiprocessing
import random

from ECC.Elliptic_Curve import Curve25519 as EllipticCurve
from NumberTheory import inverseModulo


def get_input_by_key(file_name):
    """Lấy dữ liệu từ file. Trả về: 1 dict theo chỉ mục"""
    try:
        with open(file_name, 'r') as file:
            data = {}
            for line in file:
                if ':' in line and not line.strip().startswith('#'):
                    key, value = line.strip().split(':', 1)
                    key = key.strip()
                    value = value.strip()

                    # Xử lý tuple (điểm trên đường cong)
                    if value.startswith('(') and value.endswith(')'):
                        value = eval(value)
                    else:
                        try:
                            value = int(value)
                        except ValueError:
                            pass

                    data[key] = value
            return data
    except FileNotFoundError:
        print("File không tồn tại!")
        return {}
    except ValueError as e:
        print(f"Định dạng file không đúng! {e}")
        return {}


def ECDSA_setup(bit=256, output_file="ECDSA_information.txt"):
    """
    Khởi tạo hệ thống ECDSA (Elliptic Curve Digital Signature Algorithm)

    Returns:
        curve: Đường cong Elliptic
        G: Điểm sinh (base point)
        n: Bậc của điểm G (thứ tự của nhóm con sinh bởi G)
        private_key: Khóa bí mật d (số nguyên ngẫu nhiên)
        public_key: Khóa công khai Q = d * G
    """
    import os

    # Tạo thư mục nếu chưa tồn tại
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Đã tạo thư mục: {output_dir}")

    curve = EllipticCurve()
    p = curve.p
    a = curve.A
    b = curve.B
    G = curve.G
    n = curve.n

    # Tạo khóa bí mật: số nguyên ngẫu nhiên trong [1, n-1]
    print("Đang tạo cặp khóa...")
    private_key = random.randint(1, n - 1)

    # Tạo khóa công khai: public_key = private_key * G
    public_key = curve.point_multiply(private_key, G)
    print(f"✓ Đã tạo khóa công khai")

    # Lưu thông tin
    print(f"Đang lưu thông tin vào: {output_file}")
    try:
        with open(output_file, "w", encoding="utf-8") as file:
            print(f"#Elliptic: y^2 = x^3 + {a}x + {b} (mod {p})", file=file)
            print(f"p: {p}", file=file)
            print(f"a: {a}", file=file)
            print(f"b: {b}", file=file)
            print(f"G: {G}", file=file)
            print(f"n: {n}", file=file)
            print(f"private_key: {private_key}", file=file)
            print(f"public_key: {public_key}", file=file)
        print(f"✓ Đã lưu thành công vào: {output_file}")
    except Exception as e:
        print(f"✗ LỖI khi lưu file: {e}")
        import traceback
        traceback.print_exc()

    return curve, G, n, private_key, public_key


def hash_message(message):
    """
    Hash thông điệp bằng SHA-256

    Args:
        message: Thông điệp dạng string hoặc bytes

    Returns:
        int: Giá trị hash dạng số nguyên
    """
    if isinstance(message, str):
        message = message.encode('utf-8')

    hash_obj = hashlib.sha256(message)
    hash_hex = hash_obj.hexdigest()
    hash_int = int(hash_hex, 16)

    return hash_int


def ECDSA_sign(message, curve=None, G=None, n=None, private_key=None,
               input_file="ECDSA_information.txt",
               output_file="./ECDSA_Signature/ECDSA_signature.txt"):
    """
    Ký số ECDSA

    Args:
        message: Thông điệp cần ký (string)

    Returns:
        (r, s): Chữ ký số

    Thuật toán:
        1. Tính e = hash(message)
        2. Chọn k ngẫu nhiên trong [1, n-1]
        3. Tính (x1, y1) = k * G
        4. Tính r = x1 mod n (nếu r = 0 thì chọn lại k)
        5. Tính s = k^(-1) * (e + d*r) mod n (nếu s = 0 thì chọn lại k)
        6. Chữ ký là (r, s)
    """
    import os

    # Tạo thư mục nếu chưa tồn tại
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Đọc thông tin nếu cần
    if curve is None or G is None or n is None or private_key is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve()
        G = data["G"]
        n = data["n"]
        private_key = data["private_key"]

    # Bước 1: Hash thông điệp
    e = hash_message(message)
    # Đảm bảo e nằm trong [0, n-1]
    e = e % n

    # Bước 2-5: Tính chữ ký
    while True:
        # Chọn k ngẫu nhiên
        k = random.randint(1, n - 1)

        # Tính (x1, y1) = k * G
        point = curve.point_multiply(k, G)
        if point is None:  # Điểm vô cực
            continue

        x1, y1 = point

        # Tính r = x1 mod n
        r = x1 % n
        if r == 0:
            continue

        print(k, n)
        # Tính s = k^(-1) * (e + d*r) mod n
        k_inv = inverseModulo(k, n)
        s = (k_inv * (e + private_key * r)) % n

        if s != 0:
            break  # Chữ ký hợp lệ

    # Lưu chữ ký và thông điệp
    try:
        with open(output_file, "w", encoding="utf-8") as file:
            print(f"message: {message}", file=file)
            print(f"hash: {e}", file=file)
            print(f"r: {r}", file=file)
            print(f"s: {s}", file=file)
        print(f"✓ Đã lưu chữ ký vào: {output_file}")
    except Exception as e:
        print(f"✗ Lỗi khi lưu file: {e}")

    return (r, s)


def ECDSA_verify(message, r = None, s = None, curve=None, G=None, n=None, public_key=None,
                 input_file="ECDSA_information.txt",
                 signature_file="./ECDSA_Signature/ECDSA_signature.txt",
                 output_file="./ECDSA_Signature/ECDSA_verify.txt"):
    """
    Xác minh chữ ký ECDSA

    Args:
        message: Thông điệp gốc
        r, s: Chữ ký cần xác minh

    Returns:
        bool: True nếu chữ ký hợp lệ, False nếu không

    Thuật toán:
        1. Kiểm tra r, s trong [1, n-1]
        2. Tính e = hash(message)
        3. Tính w = s^(-1) mod n
        4. Tính u1 = e*w mod n
        5. Tính u2 = r*w mod n
        6. Tính (x1, y1) = u1*G + u2*Q
        7. Nếu (x1, y1) = O thì từ chối
        8. Tính v = x1 mod n
        9. Chấp nhận nếu v = r
    """
    import os

    # Tạo thư mục nếu chưa tồn tại
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Đọc thông tin nếu cần
    if curve is None or G is None or n is None or public_key is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve()
        G = data["G"]
        n = data["n"]
        public_key = data["public_key"]

    if message is None or r is None or s is None:
        data = get_input_by_key(signature_file)
        message = data.get("message")
        r = data.get("r")
        s = data.get("s")

    # Bước 1: Kiểm tra r, s trong [1, n-1]
    if not (1 <= r < n and 1 <= s < n):
        print("✗ Chữ ký không hợp lệ: r hoặc s nằm ngoài phạm vi")
        result = False
    else:
        # Bước 2: Hash thông điệp
        e = hash_message(message)
        e = e % n

        print(s, n)
        # Bước 3: Tính w = s^(-1) mod n
        w = inverseModulo(s, n)

        # Bước 4-5: Tính u1, u2
        u1 = (e * w) % n
        u2 = (r * w) % n

        # Bước 6: Tính (x1, y1) = u1*G + u2*Q
        point1 = curve.point_multiply(u1, G)
        point2 = curve.point_multiply(u2, public_key)
        point = curve.point_add(point1, point2)

        # Bước 7: Kiểm tra điểm vô cực
        if point is None:
            print("✗ Chữ ký không hợp lệ: điểm kết quả là O")
            result = False
        else:
            x1, y1 = point

            # Bước 8-9: Tính v và so sánh với r
            v = x1 % n
            result = (v == r)

            if result:
                print("✓ Chữ ký hợp lệ!")
            else:
                print(f"✗ Chữ ký không hợp lệ: v={v} ≠ r={r}")

    # Lưu kết quả
    try:
        with open(output_file, "w", encoding="utf-8") as file:
            print(f"message: {message}", file=file)
            print(f"r: {r}", file=file)
            print(f"s: {s}", file=file)
            print(f"valid: {result}", file=file)
            if result:
                print("The signature is VALID!", file=file)
            else:
                print("The signature is INVALID!", file=file)
        print(f"✓ Đã lưu kết quả xác minh vào: {output_file}")
    except Exception as e:
        print(f"✗ Lỗi khi lưu file: {e}")

    return result


# ===== DEMO =====
if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    print("=" * 70)
    print("DEMO: ECDSA (Elliptic Curve Digital Signature Algorithm)")
    print("=" * 70)

    # # 1. Setup hệ thống
    # print("\n=== BƯỚC 1: Khởi tạo hệ thống ECDSA ===")
    # bit = 256
    # curve, G, n, private_key, public_key = ECDSA_setup(bit)
    # print(f"Đường cong: {curve}")
    # print(f"Điểm sinh G: {G}")
    # print(f"Khóa bí mật d: {private_key}")
    # print(f"Khóa công khai Q: {public_key}")

    # 2. Ký thông điệp
    print("\n=== BƯỚC 2: Ký thông điệp ===")
    message = "NguyenHaiNam23021643"
    print(f"Thông điệp: {message}")

    r, s = ECDSA_sign(message)
    print(f"Chữ ký (r, s): ({r}, {s})")

    # 3. Xác minh chữ ký (đúng)
    print("\n=== BƯỚC 3: Xác minh chữ ký (thông điệp đúng) ===")
    is_valid = ECDSA_verify(message)
    print(f"Kết quả: {'HỢP LỆ' if is_valid else 'KHÔNG HỢP LỆ'}")

    # 4. Xác minh chữ ký (sai)
    print("\n=== BƯỚC 4: Xác minh chữ ký (thông điệp sai) ===")
    fake_message = "NguyenHaiNam99999999"
    print(f"Thông điệp giả mạo: {fake_message}")
    is_valid_fake = ECDSA_verify(fake_message)
    print(f"Kết quả: {'HỢP LỆ' if is_valid_fake else 'KHÔNG HỢP LỆ'}")

    # 5. Tổng kết
    print("\n=== TỔNG KẾT ===")
    print(f"✓ Chữ ký đúng: {'Đã xác minh' if is_valid else 'Thất bại'}")
    print(f"✓ Chữ ký giả mạo: {'Đã phát hiện' if not is_valid_fake else 'Không phát hiện được'}")

    print("=" * 70)