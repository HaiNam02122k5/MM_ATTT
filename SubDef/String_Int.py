from gmpy2 import mpz
from NumberTheory import modulo
ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
BASE = len(ALPHABET)

def text_to_int(text):
    """Chuyển chuỗi (A-Z) sang số nguyên trong hệ 10, coi là cơ số 26"""
    if not text:
        return 0
    result = mpz(0)
    for char in text:
        result = result * 256 + ord(char)
    return result

def int_to_text(number):
    """Chuyển số nguyên về chuỗi A-Z theo cơ số 26"""
    if number == 0:
        return ""
    number = mpz(number)
    chars = []
    while number > 0:
        chars.append(chr(modulo(number, 256)))
        number //= 256
    return ''.join(reversed(chars))