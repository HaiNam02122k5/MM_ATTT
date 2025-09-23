import multiprocessing
from gmpy2 import mpz

from Prime_All import general_prime_bit, general_prime_in_range
from NumberTheory import moduloPower, inverseModulo, gcd


def get_input_by_key(file_name):
    """ Lấy dữ liệu từ file
        Trả về: 1 dict theo chỉ mục"""
    try:
        with open(file_name, 'r') as file:
            data = {}
            for line in file:
                key, value = line.strip().split(':')
                data[key.strip()] = int(value.strip())
    except FileNotFoundError:
        print("File không tồn tại!")
    except ValueError:
        print("Định dạng file không đúng!")

    return data

def RSA_algorithm(bit=40, output_file = "RSA_information.txt"):
    """ Khởi tạo các tham số cho RSA """
    p = general_prime_bit(bit)
    q = general_prime_bit(bit)

    # Tạo n
    n = p * q

    # Tạo phi_n là số lượng số dương nhỏ hơn và nguyên tố với n
    # Do n là số nguyên tố nên phi_n bằng tích p-1 và q-1
    phi_n = (p - 1) * (q - 1)

    # Chọn e nguyên tố cùng nhau với phi_n và nhỏ hơn phi_n
    e = general_prime_in_range(1, phi_n)
    while gcd(e, phi_n) != 1:
        e = general_prime_bit(2 * bit - bit / 2)

    d = inverseModulo(e, phi_n)

    try:
        with open(output_file, "w") as file:
            print("p:", p, file=file)
            print("q:", q, file=file)
            print("n:", n, file=file)
            print("e:", e, file=file)
            print("d:", d, file=file)
    except FileNotFoundError:
        print("File không tồn tại!")

    return n, p, q, phi_n, e, d

def RSA_encrypt(x, n = None, e = None, input_file = "RSA_information.txt", output_file = "RSA_encrypted.txt"):
    data = get_input_by_key(input_file)
    if n is None:
        n = mpz(data["n"])
    if e is None:
        e = mpz(data["e"])

    result = moduloPower(x, e, n)

    try:
        with open(output_file, "w") as file:
            print("cypher_text:", result, file=file)
    except FileNotFoundError:
        print("File dell tồn tại!")

    return result

def RSA_decrypt(y = None, n = None, d = None, input_file = "RSA_information.txt", cypher_text = "cypher_text.txt", output_file = "RSA_decrypted.txt"):
    data = get_input_by_key(input_file)
    if n is None:
        n = mpz(data["n"])
    if d is None:
        d = mpz(data["d"])
    if y is None:
        data = get_input_by_key(cypher_text)
        y = mpz(data["cypher_text"])
    result = moduloPower(y, d, n)
    try:
        with open(output_file, "w") as file:
            print("plaintext:", result, file=file)

    except FileNotFoundError:
        print("File dell tồn tại!")
    return result

if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    bit = 2048

    print(RSA_decrypt())















