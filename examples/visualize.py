#!/usr/bin/env python2
#   visualize.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from argparse import ArgumentParser
from sys import path

try:
    import seaborn as sns
except ImportError:
    pass

# Extend path variable in order to import fortran routines
NLS_MODULE_PATH = '../'

if NLS_MODULE_PATH not in path:
    path.append(NLS_MODULE_PATH)

from nls.model import Solution


def main():
    parser = ArgumentParser(description='Visualize saved solutions.')
    parser.add_argument('solution_path', metavar='solution', type=unicode, help='.mat file that contains solution')
    parser.add_argument('--1d', '-1', action='store_true')
    parser.add_argument('--2d', '-2', action='store_true')
    parser.add_argument('--default', action='store_true')
    parser.add_argument('--contour', action='store_true')
    parser.add_argument('--stream', action='store_true')

    args = parser.parse_args()

    solution = Solution.load(args.solution_path)
    solution.visualize()
    solution.show()


if __name__ == '__main__':
    main()