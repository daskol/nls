#   nls/pumping.py
#   The module defines abstract type and operation with pumping objects. Also it provides several predefined pumping
#   classes.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from numpy import exp


class OpSumPumping(object):

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def __call__(self, x, t=None):
        return self.lhs(x, t) + self.rhs(x, t)

    def __repr__(self):
        return repr(self.lhs) + u' + ' + repr(self.rhs)


class OpSubPumping(object):

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def __call__(self, x, t=None):
        return self.lhs(x, t) + self.rhs(x, t)

    def __repr__(self):
        return repr(self.lhs) + u' - ' + repr(self.rhs)


class AbstractPumping(object):
    """Base class of pumping tree that define commont interface of pumping objects behavior.
    """

    def __add__(self, other):
        return OpSumPumping(self, other)

    def __sub__(self, other):
        return OpSubPumping(self, other)

    def __str__(self):
        return str(repr(self))

    def __repr__(self):
        return u'AbstractPumping'

    def __unicode__ (self):
        return repr(self);


class GaussianPumping(AbstractPumping):
    """Steady state gaussian pumping with given origin, maximum power, and decay.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0):
        super(AbstractPumping, self).__init__()

        self.power = power
        self.x0 = x0
        self.y0 = y0
        self.variation = variation

    def __call__(self, x, y=None, t=None):
        y = 0.0 if y is None else y
        return self.power * exp( - ((x - self.x0) ** 2 + (y - self.y0) ** 2) / (2.0 * self.variation))

    def __repr__(self):
        return u'{0} exp(-{1} (x - {2})^2)'.format(self.power, 1.0 / (2.0 * self.variation), self.x0)
