#
#   __init__.py
#
#   (c) Daniel Bershatsky, 2016
#

from __future__ import print_function

__all__ = ['native']

from pprint import pprint
from time import time
from numpy import array, exp, arange
from scipy.io import loadmat, savemat
from matplotlib.pyplot import plot, show
from native import nls


class Problem(object):
    """
    """
    
    def __init__(self):
        pass

    def model(self, *args, **kwargs):
        """
        Piority of Arguments: Arguments passed in `kwargs` has the most piority, 'param' key in `kwargs` has less
        piority than `kwargs` and dictionary arguments in `args` have the least piority. Other arguments are ignored.
        Argument List:
            model - set model type, default value 'default';
            dx - default value '1.0e-1';
            dt - default value '1.0e-3';
            t0 - default value '0.0';
            u0 - default value '1.0e-1';
            order - default value '5';
            pumping - default value ``;
            !original_params - default value `{}`;
            !dimless_params - default value `{}`;
        """
        if 'params' in kwargs:
            params = kwargs.pop('params')
            #kwargs = {**params, **kwargs}  # python 3

        kwargs['model'] == 'default' if 'model' not in kwargs else kwargs['model']
        kwargs['dx'] = 1.0e-1 if 'dx' not in kwargs else kwargs['dx']
        kwargs['dt'] = 1.0e-3 if 'dt' not in kwargs else kwargs['dt']
        kwargs['t0'] = 0.0e+0 if 't0' not in kwargs else kwargs['t0']
        kwargs['u0'] = 1.0e-1 if 'u0' not in kwargs else kwargs['u0']
        kwargs['order'] = 5 if 'order' not in kwargs else kwargs['order']
        kwargs['pumping'] = GaussianPumping() if 'pumping' not in kwargs else kwargs['pumping']
        kwargs['num_nodes'] = 1000 if 'num_nodes' not in kwargs else kwargs['num_nodes']
        kwargs['num_iters'] = 100000 if 'num_iters' not in kwargs else kwargs['num_iters']

        if 'model' in kwargs and kwargs['model'] == 'default':
            return Model(**kwargs)
        else:
            raise Exception('Unknown model passed!')


class Model(object):
    """Default model that is NLS equation with reservoire in axe symmentic case.
    """

    def __init__(self, *args, **kwargs):
        pprint(kwargs)
        if 'original_params' in kwargs:
            pprint(kwargs['original_params'])
        self.solution = Solution(kwargs['dt'], kwargs['dx'], kwargs['num_nodes'], kwargs['order'], kwargs['num_iters'])
        self.solver = Solver(self.solution)

    def solve(self):
        return self.solver()

    def store(self):
        pass

    def restore(self):
        pass


class Solution(object):

    def __init__(self, dt, dx, num_nodes, order, num_iters):
        self.dt = dt
        self.dx = dx
        self.order = order
        self.num_nodes = num_nodes
        self.num_iters = num_iters
        self.solution = None

    def getTimeStep(self):
        return self.dt

    def getSpatialStep(self):
        return self.dx

    def getApproximationOrder(self):
        return self.order

    def getNumberOfNodes(self):
        return self.num_nodes

    def getNumberOfIterations(self):
        return self.num_iters

    def getElapsedTime(self):
        return self.elapsed_time

    def getSolution(self):
        return self.solution

    def setSolution(self, solution):
        self.solution = solution

    def setElapsedTime(self, seconds):
        self.elapsed_time = seconds

    def visualize(self):
        x = arange(0.0, self.dx * self.num_nodes, self.dx)
        y = (self.solution.conj() * self.solution).real
        plot(x, y)
        show()

    def store(self):
        pass

    def restore(self):
        pass


class Solver(object):
    """Wrapper(or Facade template) of native fortran subroutine `solve_nls`.
    """
    
    def __init__(self, solution):
        self.solution = solution
        self.elapsed_time = 0.0

    def __call__(self):
        self.elapsed_time = -time()
        self.solution.setSolution(nls.solve_nls(
            self.solution.getTimeStep(),
            self.solution.getSpatialStep(),
            self.solution.getNumberOfNodes(),
            self.solution.getApproximationOrder(),
            self.solution.getNumberOfIterations()))
        self.elapsed_time += time()
        return self.solution


class GaussianPumping(object):
    """Steady state gaussian pumping with given origin, maximum power, and decay.
    """

    def __init__(self, power=1.0, x0=0.0, variation=1.0):
        self.power = power
        self.x0 = x0
        self.variation = variation

    def __call__(self, x, t=None):
        return self.power * exp( - (x - self.x0) ** 2 / (2.0 * variation))

    def __str__(self):
        return str(repr(self))

    def __repr__(self):
        return u'{0} exp(-{1} (x - {2})^2)'.format(self.power, 1.0 / (2.0 * self.variation), self.x0)

    def __unicode__ (self):
        return repr(self);


def test():
    model = Problem().model(
        model = 'default',
        dx = 1.0e-1,
        dt = 1.0e-3,
        t0 = 0.0,
        u0 = 1.0,
        order = 5,
        num_nodes = 1000,
        num_iters = 1000,
        pumping = GaussianPumping(),
        original_params = {
            'R': 0.05,
            'gamma': 0.566,
            'g': 0.001,
            'tilde_g': 0.011,
            'gamma_R': 10,
        },
        dimless_params = {
            'a': 0.0,
        })
    solution = model.solve()
    solution.visualize()
    solution.store()


if __name__ == '__main__':
    test()
