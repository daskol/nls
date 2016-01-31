#!/usr/bin/env python2
#   animate.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from sys import path

# Extend path variable in order to import fortran routines
NLS_MODULE_PATH = '../src/'

if NLS_MODULE_PATH not in path:
    path.append(NLS_MODULE_PATH)

from nls.model import Problem
from nls.pumping import GaussianPumping
from nls.animation import IterationIncreaseAnimation, PumpingRadiusIncreaseAnimation, PumpingPowerIncreaseAnimation


def main():
    model = Problem().model(
        model = '1d',
        dx = 1.0e-1,
        dt = 1.0e-3,
        t0 = 0.0,
        u0 = 0.1,
        order = 5,
        num_nodes = 400,
        num_iters = 10000,
        pumping = GaussianPumping(power=3.0, variation=6.84931506849),
               #+ GaussianPumping(power=5.0, x0=+5.0, variation=6.84931506849),
        original_params = {
            'R': 0.0242057488654,
            'gamma': 0.0242057488654,
            'g': 0.00162178517398,
            'tilde_g': 0.0169440242057,
            'gamma_R': 0.242057488654,
        },
        dimless_params = {
        })

    # Animate solution profile as iterations increase
    animation = IterationIncreaseAnimation(model, 10, 50)
    animation.render('point-iteration-increase.mp4')
    animation.report()

    # Animate solution profile as radius increase
    animation = PumpingRadiusIncreaseAnimation(model, 2, 0.1)
    animation.render('ring-radius-increase.mp4')
    animation.report()

    # Animate solution profile as power increase
    animation = PumpingPowerIncreaseAnimation(model, 6, 0.05)
    animation.render('point-power-increase.mp4')
    animation.report()


if __name__ == '__main__':
    main()