from Prime_All import general_prime_bit, general_prime_in_range
from Prime_All import general_prime_bit, general_prime_in_range


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

def RSA_algorithm(bit=40, output_file = "RSA_information"):
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
        e = general_prime_in_range(1, phi_n)

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

def RSA_encrypt()








