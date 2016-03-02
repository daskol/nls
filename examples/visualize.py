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
    parser = ArgumentParser(prog='NLSe Soliution Visualizer', description='Visualize saved solutions with nls package.')
    parser.add_argument('solution_path', metavar='solution', type=unicode, help='.mat file that contains solution')
    parser.add_argument('--1d', '-1', action='store_true', help='Force loading solution as one dimensional')
    parser.add_argument('--2d', '-2', action='store_true', help='Force loading solution as two dimensional')
    parser.add_argument('--default', action='store_true', help='Visualize solution in default way')
    parser.add_argument('--contour', action='store_true', help='Show contour plot')
    parser.add_argument('--stream', action='store_true', help='Show stream plot')
    parser.add_argument('--filename', '-f', metavar='path', default=None, type=unicode, help='Filename to save figure')
    parser.add_argument('--show', '-s', action='store_true', help='Show visualized solution.')
    parser.add_argument('--no-show', '-ns', action='store_true', help='Do not show visualized solution.')

    args = parser.parse_args()

    solution = Solution.load(args.solution_path)
    solution.visualize(
        filename = args.filename,
        contour = args.contour,
        stream = args.stream,
        )
    
    if not args.no_show:
        solution.show()


if __name__ == '__main__':
    main()