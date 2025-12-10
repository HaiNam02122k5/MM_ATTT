import math
from sympy import Poly, symbols
from sympy.polys.domains import ZZ

""" Kiểm tra nguyên tố AKS - Thuật toán chính xác
    Có ý nghĩa trên lý thuyết - Thực tế quá phức tạp để triển khai"""
def gcd_in_prime(a, b):
    """Tính ước chung lớn nhất của a và b (dùng thuật toán Euclid)"""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def is_perfect_power(n):
    """
    Bước 1: Kiểm tra n có phải là dạng a^b (b > 1) không.
    """
    if n < 2:
        return False

    # Giới hạn trên của b là log2(n)
    log2n = math.log2(n)

    # Thử mọi b từ 2 đến log2(n)
    for b in range(2, int(log2n) + 1):
        # Tính a = n^(1/b)
        # Dùng (val + 0.5) để làm tròn số nguyên gần nhất
        a = int(n ** (1 / b) + 0.5)

        # Kiểm tra lại xem a^b có chính xác bằng n không
        if a ** b == n:
            return True

    return False


def multiplicative_order(n, r):
    """Tìm trật tự nhân của n trong Z/rZ"""
    if gcd_in_prime(n, r) != 1:
        return 0
    k = 1
    temp = n % r
    while temp != 1:
        temp = (temp * n) % r
        k += 1
        if k > r:  # Tránh vòng lặp vô hạn
            return float('inf')
    return k


def find_r(n):
    """
    Bước 2: Tìm r nhỏ nhất sao cho O_r(n) > (log2 n)^2.
    """
    log2n = math.log2(n)
    log2n_sq = log2n ** 2

    r = 2
    while True:
        # Bỏ qua r nếu r là ước của n
        if math.gcd(n, r) != 1:
            r += 1
            continue

        # Tính O_r(n) (cấp của n mod r)
        k = 1
        result = n % r

        # Chỉ cần kiểm tra đến giới hạn log2n_sq
        while k <= log2n_sq:
            if result == 1:
                break  # Order là k, quá nhỏ. Thử r tiếp theo.

            result = (result * n) % r
            k += 1

        # Nếu k > log2n_sq, nghĩa là order > log2n_sq
        if k > log2n_sq:
            return r  # Tìm thấy r phù hợp

        r += 1  # Thử r tiếp theo


def poly_mul(p1, p2, r, n):
    """
    Nhân hai đa thức (mod x^r - 1, n).
    p1, p2 là list các hệ số.
    """
    deg1 = len(p1)
    deg2 = len(p2)
    # Kết quả luôn có bậc < r
    res = [0] * r

    for i in range(deg1):
        for j in range(deg2):
            # Lấy hệ số mod n
            coeff = (p1[i] * p2[j]) % n
            # Lấy số mũ mod r (vì x^r = 1)
            power = (i + j) % r

            res[power] = (res[power] + coeff) % n

    # Cắt bớt các số 0 ở cuối list cho gọn
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def poly_pow(poly, e, r, n):
    """
    Tính (poly)^e (mod x^r - 1, n)
    Sử dụng thuật toán Bình phương và Nhân.
    """
    # Đa thức đơn vị (hằng số 1)
    res = [1]
    temp = poly[:]  # Tạo bản sao

    while e > 0:
        if e % 2 == 1:
            # res = res * temp
            res = poly_mul(res, temp, r, n)

        # temp = temp * temp
        temp = poly_mul(temp, temp, r, n)
        e //= 2

    return res


def is_prime_aks(n):
    """
    Thuật toán AKS đầy đủ (phiên bản giáo dục).
    """
    print(f"--- BẮT ĐẦU KIỂM TRA: n = {n} ---")

    # [Bước 1] Kiểm tra lũy thừa hoàn hảo
    if is_perfect_power(n):
        print(f"Bước 1: {n} là lũy thừa hoàn hảo.")
        return False
    print("Bước 1: OK (không là lũy thừa hoàn hảo)")

    # [Bước 2] Tìm r
    r = find_r(n)
    print(f"Bước 2: Tìm được r = {r}")

    # [Bước 3] Kiểm tra gcd(a, n)
    for a in range(1, r + 1):
        g = math.gcd(a, n)
        if 1 < g < n:
            print(f"Bước 3: {n} có ước chung {g} <= r.")
            return False
    print(f"Bước 3: OK (không có ước nhỏ <= {r})")

    # [Bước 4] Kiểm tra n nhỏ
    if n <= r:
        print(f"Bước 4: {n} <= r, kết luận {n} là SNT.")
        return True
    print(f"Bước 4: OK ({n} > {r})")

    # [Bước 5] Vòng lặp chính - Kiểm tra đồng dư đa thức
    print(f"Bước 5: Bắt đầu kiểm tra đồng dư đa thức...")
    log2n = math.log2(n)
    # Giới hạn của a. Dùng r thay cho phi(r) cho đơn giản
    limit_a = int(math.sqrt(r) * log2n)

    for a in range(1, limit_a + 1):
        # 1. Tính vế trái: (x + a)^n
        # Đa thức (x + a) là [a, 1] (tức a + 1*x)
        poly_left = poly_pow([a, 1], n, r, n)

        # 2. Tính vế phải: (x^n + a)
        # x^n (mod x^r - 1) là x^(n % r)
        power = n % r
        poly_right = [a]  # Khởi tạo là đa thức hằng số 'a'
        if len(poly_right) <= power:
            # Kéo dài list ra để gán poly_right[power]
            poly_right.extend([0] * (power - len(poly_right) + 1))
        # Cộng 1 vào hệ số của x^(n%r)
        poly_right[power] = (poly_right[power] + 1) % n

        # Cắt bớt list cho chuẩn
        while len(poly_right) > 1 and poly_right[-1] == 0:
            poly_right.pop()

        # 3. So sánh
        if poly_left != poly_right:
            print(f"Bước 5: THẤT BẠI với a = {a}")
            print(f"  Vế trái: {poly_left}")
            print(f"  Vế phải: {poly_right}")
            return False

        if a % 5 == 0 or a == limit_a:  # In log cho đỡ rối
            print(f"Bước 5: ... OK với a = {a} / {limit_a}")

    # [Bước 6] Kết luận
    print(f"Bước 6: Vượt qua tất cả, {n} là số nguyên tố.")
    return True

is_prime_aks(53469821)