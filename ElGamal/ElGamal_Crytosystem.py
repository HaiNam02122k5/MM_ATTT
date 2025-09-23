import multiprocessing
import random

from gmpy2 import mpz

from NumberTheory import moduloPower, part_primitive_root
from Prime_All import general_prime_bit


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

def ElGamal_information(a, bit=40, output_file = "ElGamal_information.txt"):
    p = mpz(general_prime_bit(bit))

    list_alpha = part_primitive_root(p)
    alpha = mpz(random.choice(list_alpha))

    beta = mpz(moduloPower(alpha, a, p))

    try:
        with open(output_file, "w") as file:
            print("p:", p, file=file)
            print("alpha:", alpha, file=file)
            print("beta:", beta, file=file)
            print("a:", a, file=file)
    except FileNotFoundError:
        print("File không tồn tại!")

    return a, alpha, beta, p


def ElGamal_encrypt(x, k, alpha = None, beta = None, p = None, input_file = "ElGamal_information.txt", output_file = "ElGamal_encrypted.txt"):
    data = get_input_by_key(input_file)
    if alpha is None:
        alpha = mpz(data["alpha"])
    if beta is None:
        beta = mpz(data["beta"])
    if p is None:
        p = mpz(data["p"])

    y1 = mpz(moduloPower(alpha, k, p))
    y2 = mpz(moduloPower(x * moduloPower(beta, k, p), 1, p))

    try:
        with open(output_file, "w") as file:
            print("y1:", y1, file=file)
            print("y2:", y2, file=file)
    except FileNotFoundError:
        print("File dell tồn tại!")

    return y1, y2

def ElGamal_decrypt(y1 = None, y2 = None, a = None, p = None, input_file = "ElGamal_information.txt", cypher_text = "cypher_text.txt", output_file = "ElGamal_decrypted.txt"):
    data = get_input_by_key(input_file)
    if a is None:
        a = mpz(data["a"])
    if p is None:
        p = mpz(data["p"])
    if y1 is None or y2 is None:
        data = get_input_by_key(cypher_text)
        y1 = mpz(data["y1"])
        y2 = mpz(data["y2"])

    dk = moduloPower(y2 * moduloPower(y1, p - a -1, p), 1, p)

    try:
        with open(output_file, "w") as file:
            print("plaintext:", dk, file=file)

    except FileNotFoundError:
        print("File dell tồn tại!")
    return dk


if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    ElGamal_decrypt()

