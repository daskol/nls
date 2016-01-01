#!/usr/bin/env python2
#
#   (c) Daniel Bershatsky <daniel.bershatsky@skolkovotech.ru>, 2015
#

from __future__ import print_function
from matplotlib.pyplot import plot, show

with open('solution') as f:
    u = [float(x) for x in f.read().split()]

plot(u)
show()
