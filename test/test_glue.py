#!/usr/bin/env python2
#
#   (c) Daniel Bershatsky <daniel.bershatsky@skolkovotech.ru>, 2015
#

from __future__ import print_function
from numpy import array, arange, diag, abs, max
from scipy.sparse import diags
from mod import mod


def test_make_banded_matrix_real():
    pass

def test_banded_matvec():
    pass

def test_make_diff2_radial_real():
    n, m, h = 8, 8, 0.1
    r = arange(0, h * m - 1e-10, h)
    D1 = diags((1, -8, 0, 8, -1), (-2, -1, 0, 1, 2), shape=(m, m)).toarray()
    D2 = diags((-1, 16, -30, 16, -1), (-2, -1, 0, 1, 2), shape=(m, m)).toarray()
    r[0] = 1
    D1[0, :] = 0
    D1[1, 1] += 1
    D2[0, :3] = [-60, 64, -4]
    D2[1, 1] += -1
    D = D2 / (24 * h ** 2) + diag(1.0 / r).dot(D1)  / (12 * h)
    L = mod.make_diff2_radial_real(n, 5, h)
    tolerance = 5.0e-5
    for k in xrange(0, 3):
        err = max(abs(diag(D, 2 - k) - L[k, 2 - k:]))
        print('Diagonal # ', 2 - k, ':    ', err)
        if err >= tolerance:
            print(diag(D, 2 - k))
            print(L[k, 2 - k:])
    for k in xrange(1, 3):
        err = max(abs(diag(D, -k) - L[2 + k, :-k]))
        print('Diagonal #', -k, ':    ', err)
        if err >= tolerance:
            print(diag(D, -k))
            print(L[2 + k, :-k])

def test_solve_nls():
    u0 = array([1, 1, 1, 1, 1, 1, 1, 1])
    op = mod.make_diff2_radial_complex(8, 5, 0.01)
    u = mod.solve_nls(10, 0.00001, 0.0, 0.01, u0, 8, 5, op)
    print(u)

if __name__ == '__main__':
    test_make_banded_matrix_real()
    test_banded_matvec()
    test_make_diff2_radial_real()
