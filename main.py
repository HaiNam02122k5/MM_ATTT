import multiprocessing
import time

from mpmath.libmp import bitcount

import Prime_All

"""
Dùng đa luồng ở Prime_All thì nhờ thêm dòng sau
multiprocessing.set_start_method('spawn', force=True)
"""

if __name__ == "__main__":
    # Đảm bảo code multiprocessing chạy trong if __name__ == '__main__'
    multiprocessing.set_start_method('spawn', force=True)  # Force spawn cho Windows

    print(bitcount(78041869916343559787769434555109302616513091073697565918482195704926881125447))
