import multiprocessing
import random
import matplotlib.pyplot as plt

from gmpy2 import mpz

from Prime_All import general_prime_bit, prime_check
from NumberTheory import moduloPower, modulo

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def x(self):
        return self.x
    def y(self):
        return self.y

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

def cal_a_b(a, b, p):
    res = modulo(4 * moduloPower(a, 3, p) + 27 * moduloPower(b, 2, p), p)
    if res == 0:
        return "lam lai"
    return a, b

def build_elliptic_curve(a, b, p):
    """Tạo 1 đường cong Elliptic y^2 = x^3 + ax + b với số điểm nguyên tố"""
    a = mpz(a)
    b = mpz(b)
    p = mpz(p)
    if modulo(4 * moduloPower(a, 3, p) + 27 * moduloPower(b, 2, p), p) == 0:
        return "4x^3 + 27b^2 = 0 (mod p)"

    q_p = {}
    q_p[0] = [0]
    for i in range (1, (p + 1) // 2):
        q_p[moduloPower(i, 2, p)] = [i, p - i]

    list_point = {}
    for i in range(0, p):
        r = modulo(moduloPower(i, 3, p) + a * i + b, p)
        if r in q_p:
            list_point[i] = q_p[r]

    num_point = len(list_point) * 2 - 1
    if not prime_check(len(list_point) * 2 - 1):
        print(str(len(list_point)) + " - Số điểm không nguyên tố!")
        return {}

    return num_point, list_point

def draw_Elliptic_point(list_point):
    for key in list_point:
        for y in list_point[key]:
            plt.scatter(key, y, s = 5, color="black")

    plt.grid(True)

    plt.show()


def ECC_information():
    return None

def ECC_encrypt():
    return None

def ECC_decrypt():
    return None

if __name__ == '__main__':
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.freeze_support()

    num_point,list_point = build_elliptic_curve(2, 12, 127)
    print(num_point)
    draw_Elliptic_point(list_point)
