#!/usr/bin/env python2
#   solve2d.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from sys import path

# Extend path variable in order to import fortran routines
NLS_MODULE_PATH = '../'

if NLS_MODULE_PATH not in path:
    path.append(NLS_MODULE_PATH)

from nls.animation import IterationIncreaseAnimation
from nls.model import Problem
from nls.pumping import GaussianPumping2D, GaussianRingPumping2D


def main():
    model = Problem().model(
        model = '2d',
        dx = 1.0e-1,
        dt = 1.0e-3,
        t0 = 0.0,
        u0 = 0.1,
        order = 5,
        num_nodes = 400,
        num_iters = 5000,
        pumping = GaussianRingPumping2D(power=20.0, radius=10.0, variation=3.14),
        original_params = {
            'R': 0.0242057488654,
            'gamma': 0.0242057488654,
            'g': 0.00162178517398,
            'tilde_g': 0.0169440242057,
            'gamma_R': 0.242057488654
        },
        dimless_params = {
        })

    # Obtain steady state solution
    solution = model.solve()
    solution.report()
    solution.visualize()
    solution.show()


if __name__ == '__main__':
    main()
