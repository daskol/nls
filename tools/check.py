#!/usr/bin/env python2
#   check.py
#   (c) Daniel Bershatsky, 2016
#   See LISENCE for datailes.

from __future__ import absolute_import, print_function
from argparse import ArgumentParser
from copy import deepcopy
from datetime import datetime
from pprint import pprint
from sys import path

from numpy import linspace, pi, zeros, mean, std
from numpy.fft import fft
from matplotlib.pyplot import figure, show

path.append('..')  # root of the repository

from nls.model import Problem, Solution


TOLERANCE = 1.0e-1


def timestamp():
    return datetime.now().replace(microsecond=0).isoformat()

def check_convergence(solution, num_steps=20, num_iters=5, tolerance=TOLERANCE):
    model = deepcopy(solution.getModel())
    profile = solution
    eigenvalues = zeros(num_steps, dtype=complex)
    integrals = zeros(num_steps, dtype=float)

    for i in xrange(num_steps):
        model.setInitialSolution(profile.getSolution())
        profile = model.solve(num_iters)
        eigenvalues[i] = model.getChemicalPotential(profile.getSolution())
        integrals[i] = profile.getDampingIntegral()

    fig = figure()
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(eigenvalues.real, label='real: {0:5.2f}'.format(eigenvalues.real[-1]))
    ax.plot(eigenvalues.imag, label='imag: {0:5.2f}'.format(eigenvalues.imag[-1]))
    ax.plot(abs(eigenvalues), label='abs: {0:5.2f}'.format(abs(eigenvalues[-1])))
    ax.legend(loc='best')
    ax.set_xlabel(str(num_iters) + '$\cdot$iters ')
    ax.set_ylabel('chemical potential, $\mu$')
    ax.set_title('Chemical potential on iterations(precalculations: {0} iters).'.format(solution.getModel().getNumberOfIterations()))
    ax = fig.add_subplot(1, 2, 2)
    ax.plot(integrals, label='value')
    ax.legend(loc='best')
    ax.set_xlabel(str(num_iters) + '$\cdot$iters ')
    ax.set_ylabel('integral value, I')
    show()

    value = mean(eigenvalues.real)
    error = std(eigenvalues.real)

    print(timestamp(), 'Mean value:', value)
    print(timestamp(), 'Standard error:', error)
    print()

    return std(eigenvalues) < tolerance

def check_integral(solution, tolerance=TOLERANCE):
    integral = solution.getDampingIntegral()

    print(timestamp(), 'Integral value:', integral)
    print()

    return abs(integral) <= tolerance

def main():
    parser = ArgumentParser(prog='NLSe Solution Checker',
                            description='Check convergence and consistance of saved solutions with nls package.')
    parser.add_argument('solution_path', metavar='solution', type=unicode, help='.mat file that contains solution')
    parser.add_argument('--1d', '-1', action='store_true', help='Force loading solution as one dimensional')
    parser.add_argument('--2d', '-2', action='store_true', help='Force loading solution as two dimensional')
    parser.add_argument('--verbose', '-V', action='store_true')
    parser.add_argument('--check-convergence', '-c', action='store_true', help='Check whether eigenvalue stabilzed')
    parser.add_argument('--check-integral', '-i', action='store_true', help='Check complex-valued integrand is zero')

    args = parser.parse_args()

    if not args.check_integral and not args.check_convergence:
        args.check_integral = args.check_convergence = True

    model = Problem().model(filename=args.solution_path)
    solution = Solution(model).restore(args.solution_path)

    print(timestamp(), 'Model description:')
    print(model)

    checks = {}

    if args.check_convergence:
        checks['convergence'] = check_convergence(solution)
    if args.check_integral:
        checks['integral'] = check_integral(solution)

    for key in sorted(checks.keys()):
        print(timestamp(), '{0} check: {1}'.format(key.capitalize(), checks[key]))


if __name__ == '__main__':
    main()