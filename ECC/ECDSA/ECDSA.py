import os
import hashlib
import multiprocessing
from gmpy2 import mpz
from ECC.Elliptic_Curve import Curve25519
from NumberTheory import inverseModulo


def get_input_by_key(file_name):
    """ Lấy dữ liệu từ file
        Trả về: 1 dict theo chỉ mục"""
    try:
        with open(file_name, 'r') as file:
            data = {}
            for line in file:
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    try:
                        data[key.strip()] = int(value.strip())
                    except ValueError:
                        data[key.strip()] = value.strip()
    except FileNotFoundError:
        print("File không tồn tại!")
        return {}
    except ValueError:
        print("Định dạng file không đúng!")
        return {}
    return data


def ECDSA_generate_keys(curve=None, output_file="ECDSA_information.txt"):
    """Khởi tạo các tham số cho ECDSA với Ed25519"""
    if curve is None:
        curve = Curve25519("ed25519")

    # Tạo private key
    private_key = curve.generate_private_key()

    # Tạo public key
    public_key = curve.get_public_key(private_key)

    try:
        with open(output_file, "w") as file:
            print("private_key:", private_key, file=file)
            if curve.curve_type == "ed25519":
                print("public_key_x:", public_key[0], file=file)
                print("public_key_y:", public_key[1], file=file)
            else:
                print("public_key_u:", public_key, file=file)
            print("curve:", curve.curve_type, file=file)
    except Exception as e:
        print(f"Lỗi khi ghi file: {e}")

    return private_key, public_key, curve


def ECDSA_sign(message, private_key=None, curve=None,
               input_file="ECDSA_information.txt",
               message_file="./ECDSA_Signature/message.txt",
               output_file="./ECDSA_Signature/ECDSA_signed.txt"):
    """Ký message bằng ECDSA"""

    # Tạo thư mục nếu chưa có
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(message_file), exist_ok=True)

    no_mess = True

    # Khởi tạo curve nếu chưa có
    if curve is None:
        curve = Curve25519("ed25519")

    # Lấy private key nếu chưa có
    if private_key is None:
        data = get_input_by_key(input_file)
        private_key = data.get("private_key")
        if private_key is None:
            print("Không tìm thấy private key!")
            return None

    # Lấy message
    if message is None:
        data = get_input_by_key(message_file)
        message = data.get("message", "")
        no_mess = False

    # Hash message
    if isinstance(message, str):
        message_bytes = message.encode()
    else:
        message_bytes = str(message).encode()

    h = int.from_bytes(hashlib.sha512(message_bytes).digest(), 'little') % curve.L

    # Tạo nonce ngẫu nhiên
    k = int.from_bytes(os.urandom(32), 'little') % curve.L
    if k == 0:
        k = 1

    # Tính r = (k * G).x mod L
    R = curve.scalar_mult(k)
    r = R[0] % curve.L

    # Tính s = k^(-1) * (h + r * private_key) mod L
    k_inv = int(inverseModulo(k, curve.L))
    s = (k_inv * (h + r * private_key)) % curve.L

    try:
        with open(output_file, "w") as file:
            print("signature_r:", r, file=file)
            print("signature_s:", s, file=file)

        if no_mess:
            with open(message_file, "w") as file:
                print("message:", message, file=file)
    except Exception as e:
        print(f"Lỗi khi ghi file: {e}")

    return (r, s)


def ECDSA_verify(message=None, signature=None, public_key=None, curve=None,
                 input_file="ECDSA_information.txt",
                 message_file="./ECDSA_Signature/message.txt",
                 signature_file="./ECDSA_Signature/signature.txt",
                 output_file="./ECDSA_Signature/ECDSA_verify.txt"):
    """Xác thực chữ ký ECDSA"""

    # Tạo thư mục nếu chưa có
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Khởi tạo curve nếu chưa có
    if curve is None:
        curve = Curve25519("ed25519")

    # Lấy public key
    if public_key is None:
        data = get_input_by_key(input_file)
        pub_x = data.get("public_key_x")
        pub_y = data.get("public_key_y")
        if pub_x is None or pub_y is None:
            print("Không tìm thấy public key!")
            return False
        public_key = (pub_x, pub_y)

    # Lấy message
    if message is None:
        data = get_input_by_key(message_file)
        message = data.get("message", "")

    # Lấy signature
    if signature is None:
        data = get_input_by_key(signature_file)
        r = data.get("signature_r")
        s = data.get("signature_s")
        if r is None or s is None:
            print("Không tìm thấy signature!")
            return False
        signature = (r, s)

    r, s = signature

    # Kiểm tra r, s trong khoảng hợp lệ
    if not (0 < r < curve.L and 0 < s < curve.L):
        print("Signature không hợp lệ!")
        return False

    # Hash message
    if isinstance(message, str):
        message_bytes = message.encode()
    else:
        message_bytes = str(message).encode()

    h = int.from_bytes(hashlib.sha512(message_bytes).digest(), 'little') % curve.L

    # Tính w = s^(-1) mod L
    w = int(inverseModulo(s, curve.L))

    # Tính u1 = h * w mod L
    u1 = (h * w) % curve.L

    # Tính u2 = r * w mod L
    u2 = (r * w) % curve.L

    # Tính điểm P = u1 * G + u2 * public_key
    P1 = curve.scalar_mult(u1)
    P2 = curve.scalar_mult(u2, public_key)
    P = curve.point_add(P1, P2)

    if P is None:
        result = False
    else:
        # Kiểm tra P.x mod L == r
        result = (P[0] % curve.L) == r

    try:
        with open(output_file, "w") as file:
            print("Result:", result, file=file)
            if result:
                print("The signature is True!", file=file)
            else:
                print("The signature is False!", file=file)
    except Exception as e:
        print(f"Lỗi khi ghi file: {e}")

    return result


if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    # print("=== ECDSA with Curve25519 Class ===")
    #
    # Tạo đường cong Ed25519
    curve = Curve25519("ed25519")
    print(f"\nCurve: {curve}")
    curve.print_curve_info()

    # Tạo cặp khóa
    print("\n1. Generating keys...")
    private_key, public_key, curve = ECDSA_generate_keys(curve)
    print(f"Private key: {private_key}")
    print(f"Public key: ({public_key[0]}, {public_key[1]})")

    # Ký message
    print("\n2. Signing message...")
    message = "NguyenHaiNam23021643"
    signature = ECDSA_sign(message)
    print(f"Message: {message}")
    print(f"Signature: r={signature[0]}, s={signature[1]}")

    # Xác thực chữ ký
    print("\n3. Verifying signature...")
    is_valid = ECDSA_verify(message, signature)
    print(f"Verification result: {is_valid}")

    # # Test với message sai
    # print("\n4. Testing with wrong message...")
    # wrong_message = "WrongMessage"
    # is_valid_wrong = ECDSA_verify(wrong_message, signature, public_key, curve)
    # print(f"Verification result (wrong message): {is_valid_wrong}")

    # print("=== ECDH Key Exchange ===")
    # curve_ecdh = Curve25519("curve25519")
    # curve_ecdh.print_curve_info()
    #
    # # Alice và Bob tạo khóa ECDH
    # print("Hai Nam")
    # alice_private = 100805
    # alice_public = curve_ecdh.get_public_key(alice_private)
    # print(alice_public)
    #
    # print("Hoang Anh Dep Zai")
    # bob_private = 19102005
    # bob_public = curve_ecdh.get_public_key(bob_private)
    # print(bob_public)
    #
    # # Tính shared secret
    # alice_shared = curve_ecdh.scalar_mult(alice_private, bob_public)
    # bob_shared = curve_ecdh.scalar_mult(bob_private, alice_public)
    # print(f"Shared secret established: {alice_shared == bob_shared}")
    #
    # # Derive signing keys
    # curve_sign = Curve25519("ed25519")
    # alice_signing_key = curve_sign.derive_signing_key(alice_shared, b"alice_signs")
    # alice_signing_public = curve_sign.get_public_key(alice_signing_key)
    #
    # bob_signing_key = curve_sign.derive_signing_key(bob_shared, b"bob_signs")
    # bob_signing_public = curve_sign.get_public_key(bob_signing_key)
    #
    # print("\n=== Alice Signs ===")
    # alice_message = "EMYEUTHAYLEPHEDO"
    # alice_sig = ECDSA_sign(alice_message, alice_signing_key, curve_sign,
    #                        message_file="./ECDSA_Signature/alice_message.txt",
    #                        output_file="./ECDSA_Signature/alice_signature.txt")
    # print(f"Message: {alice_message}")
    # print(f"Signature: (r={alice_sig[0]}, s={alice_sig[1]})")
    #
    # # Bob verify
    # bob_derived_alice_public = curve_sign.get_public_key(
    #     curve_sign.derive_signing_key(bob_shared, b"alice_signs"))
    # alice_valid = ECDSA_verify(alice_message, alice_sig, bob_derived_alice_public, curve_sign)
    # print(f"Bob verifies: {alice_valid}")
    #
    # print("\n=== Bob Signs ===")
    # bob_message = "TOISINHVIENK68"
    # bob_sig = ECDSA_sign(bob_message, bob_signing_key, curve_sign,
    #                      message_file="./ECDSA_Signature/bob_message.txt",
    #                      output_file="./ECDSA_Signature/bob_signature.txt")
    # print(f"Message: {bob_message}")
    # print(f"Signature: (r={bob_sig[0]}, s={bob_sig[1]})")
    #
    # # Alice verify
    # alice_derived_bob_public = curve_sign.get_public_key(
    #     curve_sign.derive_signing_key(alice_shared, b"bob_signs"))
    # bob_valid = ECDSA_verify(bob_message, bob_sig, alice_derived_bob_public, curve_sign)
    # print(f"Alice verifies: {bob_valid}")


