import multiprocessing
import random

from gmpy2 import mpz

from NumberTheory import moduloPower, part_primitive_root, modulo, inverseModulo
from Prime_All import generate_prime_bit
from SubDef.String_Int import text_to_int, int_to_text


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

def ElGamal_information(a, bit=40, p = None, output_file = "ElGamal_information.txt"):
    if p is None:
        p = mpz(generate_prime_bit(bit))

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


def ElGamal_encrypt(x, k, alpha = None, beta = None, p = None,
                    input_file = "ElGamal_information.txt",
                    output_file = "./ElGamal_Crytosystem/ElGamal_encrypted.txt"):
    x = text_to_int(x)
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

def ElGamal_decrypt(y1 = None, y2 = None, a = None, p = None,
                    input_file = "ElGamal_information.txt",
                    cypher_text = "./ElGamal_Crytosystem/cypher_text.txt",
                    output_file = "./ElGamal_Crytosystem/ElGamal_decrypted.txt"):
    data = get_input_by_key(input_file)
    if a is None:
        a = mpz(data["a"])
    if p is None:
        p = mpz(data["p"])
    if y1 is None or y2 is None:
        data = get_input_by_key(cypher_text)
        y1 = mpz(data["y1"])
        y2 = mpz(data["y2"])

    dk = int_to_text(moduloPower(y2 * moduloPower(y1, p - a -1, p), 1, p))

    try:
        with open(output_file, "w") as file:
            print("plaintext:", dk, file=file)

    except FileNotFoundError:
        print("File dell tồn tại!")
    return dk

def ElGamal_sign(x = None, k = None, a = None, alpha = None, p = None,
                 input_file = "ElGamal_information.txt",
                 message = "./ElGamal_Signature_Scheme/message.txt",
                 output_file = "./ElGamal_Signature_Scheme/ElGamal_signed.txt"):
    no_mess = True  # truyền x từ đầu vào
    data = get_input_by_key(input_file)
    if alpha is None:
        alpha = mpz(data["alpha"])
    if a is None:
        a = mpz(data["a"])

    if p is None:
        p = mpz(data["p"])

    if x is None:
        data = get_input_by_key(message)
        x = data["message"]
        no_mess = False

    if k is None:
        k = random.randint(0, 2**32 - 1)

    gama = moduloPower(alpha, k, p)
    delta = modulo(modulo(x - a * gama, p - 1) * inverseModulo(k, p - 1), p-1)

    try:
        with open(output_file, "w") as file:
            print("gama:", gama, file=file)
            print("delta:", delta, file=file)
        if no_mess:
            with open(message, "w") as file:
                print("message:", x, file=file)
    except FileNotFoundError:
        print("File dell tồn tại!")

    return gama, delta

def ElGamal_ver(x = None, alpha = None, beta = None, delta = None, gama = None, p = None,
            input_file = "ElGamal_information.txt",
            message = "./ElGamal_Signature_Scheme/message.txt",
            signature = "./ElGamal_Signature_Scheme/signature.txt",
            output_file = "./ElGamal_Signature_Scheme/ElGamal_verify.txt"):
    data = get_input_by_key(input_file)
    if alpha is None:
        alpha = data["alpha"]
    if beta is None:
        beta = data["beta"]
    if p is None:
        p = data["p"]

    data = get_input_by_key(signature)
    if delta is None:
        delta = data["delta"]
    if gama is None:
        gama = data["gama"]

    if x is None:
        data = get_input_by_key(message)
        x = data["message"]

    VT = modulo(moduloPower(beta, gama, p) * moduloPower(gama, delta, p), p)
    VP = moduloPower(alpha, x, p)
    if VT == VP:
        try:
            with open(output_file, "w") as file:
                print("True", file=file)
                print("The signature is True!", file=file)
        except FileNotFoundError:
            print("File dell tồn tại!")
        return True
    if VT == VP:
        try:
            with open(output_file, "w") as file:
                print("False", file=file)
                print("The signature is False!", file=file)
        except FileNotFoundError:
            print("File dell tồn tại!")
        return False

if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    ElGamal_encrypt("NguyenHaiNam23021643", 1008)
    ElGamal_decrypt()

