#!/usr/bin/env python2
#   __init__.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from sys import path

# Extend path variable in order to import fortran routines
NLS_MODULE_PATH = '../../src/'

if NLS_MODULE_PATH not in path:
    path.append(NLS_MODULE_PATH)

from nls import Problem, GaussianPumping


def test():
    model = Problem().model(
        model = 'default',
        dx = 1.0e-1,
        dt = 1.0e-3,
        t0 = 0.0,
        u0 = 0.1,
        order = 5,
        num_nodes = 1000,
        num_iters = 10000,
        pumping = GaussianPumping(power=3.0, variation=6.84931506849),
        original_params = {
            'R': 0.05,
            'gamma': 0.566,
            'g': 1.0e-3,
            'tilde_g': 0.011,
            'gamma_R': 10,
        },
        dimless_params = {
        })

    solution = model.solve()
    solution.visualize()
    solution.show()
    solution.store()
    solution.report()

    animation = model.animate('condensation-point.mp4')


if __name__ == '__main__':
    test()
