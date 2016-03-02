#   nls/pumping.py
#   The module defines abstract type and operation with pumping objects. Also it provides several predefined pumping
#   classes.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from numpy import arctan2, exp, sqrt, cos, sin


class AbstractPumping(object):
    """Base class of pumping tree that define commont interface of pumping objects behavior.
    """

    def __init__(self, power=1.0):
        self.power = 1.0

    def setPower(self, power):
        """Set local parameter power of pumping. Word `local` means multiplier in case of arithmetic operation
        (`OpSumPumping`, `OpSubPumping`, etc) and real power of pumping in other cases.
        """
        self.power = power

    def __call__(self, *args, **kwargs):
        raise Exception("Nothing to call: abstract class could not represent pumping.")

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


class OpSumPumping(AbstractPumping):
    """Functor object that incapsulates two pumping functional objects and represents their sum.
    """

    def __init__(self, lhs, rhs):
        super(OpSumPumping, self).__init__()

        self.lhs = lhs
        self.rhs = rhs

    def __call__(self, *args, **kwargs):
        return self.power * (self.lhs(*args, **kwargs) + self.rhs(*args, **kwargs))

    def __repr__(self):
        return repr(self.lhs) + u' + ' + repr(self.rhs)


class OpSubPumping(AbstractPumping):
    """Functor object that incapsulates two pumping functional objects and represents their substract.
    """

    def __init__(self, lhs, rhs):
        super(OpSubPumping, self).__init__()

        self.lhs = lhs
        self.rhs = rhs

    def __call__(self, *args, **kwargs):
        return self.power * (self.lhs(*args, **kwargs) - self.rhs(*args, **kwargs))

    def __repr__(self):
        return repr(self.lhs) + u' - ' + repr(self.rhs)


class OpMulPumping(AbstractPumping):
    """Functor object that incapsulates two pumping functional objects and represents their cartesian multiplication.
    For example, if one passes two one dimension pumping, it will yeild one two dimensional pumping.
    """

    def __init__(self, lhs, rhs):
        super(OpMulPumping, self).__init__()

        self.lhs = lhs
        self.rhs = rhs

    def __call__(self, *args, **kwargs):
        raise Exception("Not supported yet!")

    def __repr__(self):
        return repr(self.lhs) + u' x ' + repr(self.rhs)


class GridPumping(AbstractPumping):
    """This class represents general grid approximation of pumping profile. Also it is used to represet custom pumpings
    as well as loaded from file pumping.
    """

    def __init__(self, pumping, desciption=None):
        super(GridPumping, self).__init__()

        self.pumping = pumping
        self.desciption = desciption if desciption else '<class GridPumping>'

    def __call__(self, *args, **kwargs):
        return self.pumping

    def __repr__(self):
        return self.desciption


class GaussianPumping(AbstractPumping):
    """Steady state gaussian pumping with given origin, maximum power, and decay.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0):
        super(GaussianPumping, self).__init__()

        self.power = power
        self.x0 = x0
        self.y0 = y0
        self.variation = variation

    def __call__(self, x, y, t=None):
        return self.power * exp( - ((x - self.x0) ** 2 + (y - self.y0) ** 2) / (2.0 * self.variation ** 2))

    def __repr__(self):
        pattern = u'{0} exp(-{1} ((x - {2})^2 - (y - {3})^2))'
        return pattern.format(self.power, 1.0 / (2.0 * self.variation ** 2), self.x0, self.y0)


class GaussianPumping1D(GaussianPumping):
    """Facade object for `GaussianPumping` class in one dimension case. This overrides only parent call-method.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0):
        super(GaussianPumping1D, self).__init__(power, x0, y0, variation)

    def __call__(self, x, t=None):
        return GaussianPumping.__call__(self, x, 0.0, t)


class GaussianPumping2D(GaussianPumping):
    """Facade(actually proxy) for `GaussianPumping` class in two dimension case.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0):
        super(GaussianPumping2D, self).__init__(power, x0, y0, variation)


class GaussianRingPumping1D(OpSumPumping):
    """Implemet pumping on ring with gaussian profile  The center of a ring is in origin.
    """

    def __init__(self, power=1.0, radius=0.0, variation=5.0):
        super(GaussianRingPumping1D, self).__init__(
            GaussianPumping1D(power, +radius, 0.0, variation),
            GaussianPumping1D(power, -radius, 0.0, variation))


class GaussianRingPumping2D(AbstractPumping):
    """Provide ring pumping in 2d with any origin of a ring.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0, radius=1.0):
        super(GaussianRingPumping2D, self).__init__()

        self.x0 = x0
        self.y0 = y0
        self.pumping = GaussianRingPumping1D(power, radius, variation)

    def __call__(self, x, y, t=None):
        radii = sqrt((x - self.x0) ** 2 + (y - self.y0) ** 2)
        return self.pumping(radii, t)

    def __repr__(self):
        return '<class GaussianRingPumping2D(AbstractPumping)>'


class GaussianElipticPumping2D(AbstractPumping):
    """Generalization of GaussianRingPumping2D but ring(circle) is replaced with elipse.
    """

    def __init__(self, power=1.0, x0=0.0, y0=0.0, variation=5.0, a=1.0, b=1.0):
        super(GaussianElipticPumping2D, self).__init__(power)

        self.x0 = x0
        self.y0 = y0
        self.var = variation
        self.a = a
        self.b = b

    def __call__(self, x, y, t=None):
        t = arctan2(x - self.x0, y - self.y0)
        re = sqrt((self.a * cos(t)) ** 2 + (self.b * sin(t)) ** 2)
        rp = sqrt(x ** 2 + y ** 2)
        return self.power * (exp(-(re - rp) ** 2 / (2 * self.var)) + \
                          exp(-(re + rp) ** 2 / (2 * self.var)))

    def __repr__(self):
        return u'<class GaussianElipticPumping2D(AbstractPumping)>'