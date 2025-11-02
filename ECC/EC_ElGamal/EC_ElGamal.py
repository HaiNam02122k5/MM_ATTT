import multiprocessing
import random

from Prime_All import generate_prime_bit
from ECC.Elliptic_Curve import EllipticCurve
from SubDef.String_Point import string_to_point, point_to_string

def get_input_by_key(file_name):
    """Lấy dữ liệu từ file. Trả về: 1 dict theo chỉ mục"""
    try:
        with open(file_name, 'r') as file:
            data = {}
            for line in file:
                if ':' in line:
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


def ElGamal_ECC_setup(bit=256, output_file="EC_ElGamal_information.txt"):
    """
    Khởi tạo hệ thống ElGamal trên đường cong Elliptic

    Returns:
        curve: Đường cong Elliptic
        G: Điểm sinh (base point)
        n: Bậc của điểm G (số điểm trong nhóm cyclic sinh bởi G)
        private_key: Khóa bí mật (số nguyên ngẫu nhiên)
        public_key: Khóa công khai (điểm trên đường cong)
    """
    # Tạo số nguyên tố p
    p = generate_prime_bit(bit)

    # Chọn tham số đường cong a, b ngẫu nhiên
    a = random.randint(1, p - 1)
    b = random.randint(1, p - 1)

    # Tạo đường cong
    while True:
        try:
            curve = EllipticCurve(a, b, p)
            break
        except ValueError:
            # Nếu đường cong suy biến, chọn lại a, b
            a = random.randint(1, p - 1)
            b = random.randint(1, p - 1)

    # Tìm điểm sinh G trên đường cong
    G = curve.generate_point()

    # Ước lượng bậc của G (trong thực tế cần tính chính xác hơn)
    # Theo định lý Hasse: |n - (p+1)| <= 2*sqrt(p)
    n = p + 1  # Ước lượng đơn giản

    # Tạo khóa bí mật: số nguyên ngẫu nhiên trong [1, n-1]
    private_key = random.randint(1, n - 1)

    # Tạo khóa công khai: public_key = private_key * G
    public_key = curve.point_multiply(private_key, G)

    # Lưu thông tin
    try:
        with open(output_file, "w") as file:
            print(f"#Elliptic: y^2 = x^3 + {a}x + {b} (mod {p})", file=file)
            print(f"p: {p}", file=file)
            print(f"a: {a}", file=file)
            print(f"b: {b}", file=file)
            print(f"G: {G}", file=file)
            print(f"n: {n}", file=file)
            print(f"private_key: {private_key}", file=file)
            print(f"public_key: {public_key}", file=file)
    except Exception as e:
        print(f"Lỗi khi lưu file: {e}")

    return curve, G, n, private_key, public_key


def ElGamal_ECC_encrypt(message_point, curve=None, G=None, public_key=None, n=None,
                        input_file="EC_ElGamal_information.txt",
                        output_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_encrypted.txt"):
    """
    Mã hóa ElGamal trên đường cong Elliptic

    Args:
        message_point: Điểm biểu diễn thông điệp trên đường cong

    Returns:
        (C1, C2): Cặp điểm mã hóa
            C1 = k * G
            C2 = M + k * public_key
    """
    # Đọc thông tin nếu cần
    if curve is None or G is None or public_key is None or n is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve(a, b, p)
        G = data["G"]
        public_key = data["public_key"]
        n = data["n"]

    # Kiểm tra message_point có nằm trên đường cong không
    if not curve.is_on_curve(message_point):
        raise ValueError("Message point không nằm trên đường cong!")

    # Chọn số ngẫu nhiên k (ephemeral key)
    k = random.randint(1, n - 1)

    # Tính C1 = k * G
    C1 = curve.point_multiply(k, G)

    # Tính C2 = M + k * public_key
    k_public = curve.point_multiply(k, public_key)
    C2 = curve.point_add(message_point, k_public)

    # Lưu kết quả
    try:
        with open(output_file, "w") as file:
            print(f"C1: {C1}", file=file)
            print(f"C2: {C2}", file=file)
            print(f"# k (ephemeral key - chỉ để debug): {k}", file=file)
    except Exception as e:
        print(f"Lỗi khi lưu file: {e}")

    return (C1, C2)


def ElGamal_ECC_decrypt(C1, C2, curve=None, private_key=None,
                        input_file="EC_ElGamal_information.txt",
                        cypher_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_encrypted.txt",
                        output_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_decrypted.txt"):
    """
    Giải mã ElGamal trên đường cong Elliptic

    Args:
        C1, C2: Cặp điểm mã hóa

    Returns:
        M: Điểm biểu diễn thông điệp gốc
            M = C2 - private_key * C1
    """
    # Đọc thông tin nếu cần
    if curve is None or private_key is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve(a, b, p)
        private_key = data["private_key"]

    if C1 is None or C2 is None:
        data = get_input_by_key(cypher_file)
        C1 = data["C1"]
        C2 = data["C2"]

    # Tính S = private_key * C1
    S = curve.point_multiply(private_key, C1)

    # Tính M = C2 - S = C2 + (-S)
    # Điểm -S có cùng x, y đối dấu
    neg_S = (S[0], (-S[1]) % curve.p)
    M = curve.point_add(C2, neg_S)

    # Lưu kết quả
    try:
        with open(output_file, "w") as file:
            print(f"plaintext_point: {M}", file=file)
    except Exception as e:
        print(f"Lỗi khi lưu file: {e}")

    return M


def ElGamal_ECC_encrypt_text(text, curve=None, G=None, public_key=None, n=None,
                             input_file="EC_ElGamal_information.txt",
                             output_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_encrypted.txt"):
    """Mã hóa chuỗi văn bản"""
    # Đọc thông tin nếu cần
    if curve is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve(a, b, p)
        G = data["G"]
        public_key = data["public_key"]
        n = data["n"]

    # Chuyển text thành điểm
    message_point, offset = string_to_point(text, curve)

    # Mã hóa
    C1, C2 = ElGamal_ECC_encrypt(message_point, curve, G, public_key, n, input_file, output_file)

    # Lưu thêm offset
    try:
        with open(output_file, "a") as file:
            print(f"offset: {offset}", file=file)
    except Exception as e:
        print(f"Lỗi khi lưu file: {e}")

    return (C1, C2, offset)


def ElGamal_ECC_decrypt_text(C1=None, C2=None, offset=None, curve=None, private_key=None,
                             input_file="EC_ElGamal_information.txt",
                             cypher_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_encrypted.txt",
                             output_file="./EC_ElGamal_Cryptosystem/EC_ElGamal_decrypted.txt"):
    """Giải mã chuỗi văn bản"""
    # Đọc thông tin nếu cần
    if curve is None:
        data = get_input_by_key(input_file)
        p = data["p"]
        a = data["a"]
        b = data["b"]
        curve = EllipticCurve(a, b, p)
        private_key = data["private_key"]

    if C1 is None or C2 is None or offset is None:
        data = get_input_by_key(cypher_file)
        C1 = data["C1"]
        C2 = data["C2"]
        offset = data["offset"]

    # Giải mã
    message_point = ElGamal_ECC_decrypt(C1, C2, curve, private_key, input_file, cypher_file, output_file)

    # Chuyển điểm thành text
    text = point_to_string(message_point, offset, curve)

    # Lưu kết quả
    try:
        with open(output_file, "a") as file:
            print(f"plaintext: {text}", file=file)
    except Exception as e:
        print(f"Lỗi khi lưu file: {e}")

    return text


# ===== DEMO =====
if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    print("=" * 70)
    print("DEMO: ElGamal Encryption trên Đường cong Elliptic")
    print("=" * 70)

    # 1. Setup hệ thống
    # print("\n=== BƯỚC 1: Khởi tạo hệ thống ===")
    # bit = 256
    # curve, G, n, private_key, public_key = ElGamal_ECC_setup(bit)
    # print(f"Đường cong: {curve}")
    # print(f"Điểm sinh G: {G}")
    # print(f"Khóa bí mật: {private_key}")
    # print(f"Khóa công khai: {public_key}")

    # 2. Mã hóa văn bản
    print("\n=== BƯỚC 2: Mã hóa văn bản ===")
    plaintext = "NguyenHaiNam23021643"
    print(f"Văn bản gốc: {plaintext}")

    C1, C2, offset = ElGamal_ECC_encrypt_text(plaintext)
    print(f"C1: {C1}")
    print(f"C2: {C2}")
    print(offset)

    # 3. Giải mã
    print("\n=== BƯỚC 3: Giải mã ===")
    decrypted_text = ElGamal_ECC_decrypt_text()
    print(f"Văn bản giải mã: {decrypted_text}")

    # 4. Kiểm tra
    print("\n=== BƯỚC 4: Kiểm tra ===")
    if plaintext == decrypted_text:
        print("✓ Mã hóa và giải mã thành công!")
    else:
        print("✗ Có lỗi xảy ra!")

    print("=" * 70)