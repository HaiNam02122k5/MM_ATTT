import multiprocessing
import time

import Prime_All

"""
Dùng đa luồng ở Prime_All thì nhờ thêm dòng sau
multiprocessing.set_start_method('spawn', force=True)
"""

if __name__ == "__main__":
    # Đảm bảo code multiprocessing chạy trong if __name__ == '__main__'
    multiprocessing.set_start_method('spawn', force=True)  # Force spawn cho Windows

    print(Prime_All.general_prime_bit(2048))
