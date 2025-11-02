
def string_to_point(text, curve):
    """
    Chuyển đổi chuỗi thành điểm trên đường cong
    Phương pháp: Encode text thành số, sau đó tìm điểm có x = số đó
    """
    # Chuyển text thành số
    number = 0
    for char in text:
        number = number * 256 + ord(char)

    # Tìm điểm có x gần với number
    for offset in range(100):  # Thử offset 0-99
        x = (number + offset) % curve.p
        y_squared = (pow(x, 3, curve.p) + curve.a * x + curve.b) % curve.p

        y = curve._mod_sqrt(y_squared, curve.p)
        if y is not None:
            return (x, y), offset

    raise ValueError("Không thể encode text thành điểm!")


def point_to_string(point, offset, curve):
    """
    Chuyển đổi điểm trên đường cong thành chuỗi
    """
    x, y = point
    # Trừ offset để lấy lại số gốc
    number = (x - offset) % curve.p

    # Chuyển số thành text
    chars = []
    while number > 0:
        chars.append(chr(number % 256))
        number //= 256

    return ''.join(reversed(chars))
