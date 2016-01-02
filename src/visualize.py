#!/usr/bin/env python2
#
#   (c) Daniel Bershatsky <daniel.bershatsky@skolkovotech.ru>, 2015
#

from __future__ import print_function
from numpy import arange
from matplotlib.pyplot import plot, show, legend
from sys import argv


def visualize(filename='solution'):
    with open(filename) as f:
        u = [float(x) for x in f.read().split()]
    
    n, h = len(u), 0.1
    plot(arange(0.0, n * h, h), u, label=filename)

def main():
    visualize(argv[1] if len(argv) == 2 else 'solution')
    legend(loc='best')
    show()


if __name__ == '__main__':
    main()