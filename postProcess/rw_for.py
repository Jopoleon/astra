import numpy as np

def wr_for(arr_in, fmt='%13.6E', n_lin=6):
    arr_flat = arr_in.T.ravel()
    lines = [
        ''.join(fmt % val for val in arr_flat[i:i + n_lin])
        for i in range(0, len(arr_flat), n_lin)
    ]
    return '\n'.join(lines) + '\n'
