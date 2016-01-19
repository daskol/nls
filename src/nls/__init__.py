#
#   __init__.py
#
#   (c) Daniel Bershatsky, 2016
#

from __future__ import print_function

__all__ = ['native']

from pprint import pprint
from time import time
from datetime import datetime
from numpy import array, exp, arange, ones, zeros
from scipy.io import loadmat, savemat
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, subplot, legend, xlim, ylim, contourf
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
        self.solution = Solution(kwargs['dt'], kwargs['dx'], kwargs['num_nodes'], kwargs['order'], kwargs['num_iters'],
                             kwargs['pumping'], kwargs['original_params'])
        self.solver = Solver(self.solution)

    def solve(self):
        return self.solver()

    def store(self, filename=None, label='', desc='', date=datetime.now()):
        """Store object to mat-file. TODO: determine format specification
        """
        filename = filename if filename else str(date).replace(' ', '_') + '.mat'

        matfile = {}
        matfile['desc'] = desc
        matfile['label'] = label
        matfile['originals'] = {}

        savemat(filename, matfile)

    def restore(self, filename):
        """Restore object from mat-file. TODO: determine format specification
        """
        matfile = loadmat(filename)

        self.desc = str(matfile['desc'][0]) if matfile['desc'].size else ''
        self.label = str(matfile['label'][0]) if matfile['label'].size else ''
        self.originals = {}

        return self


class Solution(object):

    t0 = 1.0e+0 # seconds

    def __init__(self, dt, dx, num_nodes, order, num_iters, pumping, originals):
        self.dt = dt
        self.dx = dx
        self.order = order
        self.num_nodes = num_nodes
        self.num_iters = num_iters
        self.pumping = pumping
        self.solution = None
        self.originals = originals
        self.coeffs = zeros(23)
        self.elapsed_time = 0.0

        # NLS equation coeficients
        self.coeffs[0] = 1.0  # \partial_t
        self.coeffs[1] = 1.0  # \nabla^2
        self.coeffs[2] = originals['R'] / (4.0 * originals['tilde_g'])  # 
        self.coeffs[3] = originals['gamma'] * Solution.t0 / 2  # linear damping
        self.coeffs[4] = 1.0  # nonlinearity
        self.coeffs[5] = 1.0  # interaction to reservoir

        # Reservoir equation coefficients
        self.coeffs[10] = 0.0  # \parital_t
        self.coeffs[11] = 2.0 * originals['tilde_g'] * Solution.t0 / originals['gamma_R']  # pumping coefficient
        self.coeffs[12] = 1.0  # damping
        self.coeffs[13] = originals['R'] / (originals['gamma_R'] * originals['g'])  # interaction term
        self.coeffs[14] = 0.0  # diffusive term

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

    def getPumping(self):
        return self.pumping(arange(0, self.dx * self.num_nodes, self.dx))

    def getCoefficients(self):
        return ones(23)

    def getSolution(self):
        return self.solution

    def getElapsedTime(self):
        return self.elapsed_time

    def setSolution(self, solution):
        self.solution = solution

    def setElapsedTime(self, seconds):
        self.elapsed_time = seconds

    def visualize(self):
        x = arange(0.0, self.dx * self.num_nodes, self.dx)
        p = self.pumping(x)  # pumping profile
        u = (self.solution.conj() * self.solution).real  # density profile
        n = self.coeffs[11] *  p / (self.coeffs[12] + self.coeffs[13] * u)

        def rect_plot(subplot_number, value, label, name, labelx, labely, xmax=20):
            subplot(2, 3, subplot_number)
            plot(x, p, label=label)
            xlim((0, xmax))
            legend(loc='best')
            title(name)
            xlabel(labelx)
            ylabel(labely)

        rect_plot(1, p, 'pumping', 'Pumping profile.', 'r', 'p')
        rect_plot(2, u, 'density', 'Density distribution of BEC.', 'r', 'u')
        rect_plot(3, n, 'reservoir', 'Density distribution of reservoir.', 'r', 'n')

        def polar_plot(subplot_number, value, xmax=20):
            subplot(2, 3, subplot_number, polar=True)
            theta = arange(0, 2 * 3.14 + 0.1, 0.1)
            contourf(theta, x, array([p for _ in theta]).T)
            ylim((0, xmax))

        polar_plot(4, p)
        polar_plot(5, u)
        polar_plot(6, n)

        show()

    def store(self, filename=None, label='', desc='', date=datetime.now()):
        """Store object to mat-file. TODO: determine format specification
        """
        filename = filename if filename else str(date).replace(' ', '_') + '.mat'

        matfile = {}
        matfile['desc'] = desc
        matfile['dimlesses'] = self.coeffs
        matfile['elapsed_time'] = self.elapsed_time
        matfile['label'] = label
        matfile['originals'] = self.originals
        matfile['pumping'] = self.getPumping()
        matfile['solution'] = self.solution

        savemat(filename, matfile)

    def restore(self, filename):
        """Restore object from mat-file. TODO: determine format specification
        """
        matfile = loadmat(filename)

        self.desc = str(matfile['desc'][0]) if matfile['desc'].size else ''
        self.coeffs = matfile['dimlesses'].T
        self.elapsed_time = matfile['elapsed_time']
        self.label = str(matfile['label'][0]) if matfile['label'].size else ''
        self.originals = {}
        self.solution = matfile['solution'].T

        return self

    def report(self):
        message = 'Elapsed in {0} seconds with {1} iteration on {2} grid nodes.'
        print(message.format(self.elapsed_time, self.num_iters, self.num_nodes))


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
            self.solution.getApproximationOrder(),
            self.solution.getNumberOfIterations(),
            self.solution.getPumping(),
            self.solution.getCoefficients()))
        self.elapsed_time += time()
        self.solution.setElapsedTime(self.elapsed_time)
        return self.solution


class GaussianPumping(object):
    """Steady state gaussian pumping with given origin, maximum power, and decay.
    """

    def __init__(self, power=1.0, x0=0.0, variation=5.0):
        self.power = power
        self.x0 = x0
        self.variation = variation

    def __call__(self, x, t=None):
        return self.power * exp( - (x - self.x0) ** 2 / (2.0 * self.variation))

    def __str__(self):
        return str(repr(self))

    def __repr__(self):
        return u'{0} exp(-{1} (x - {2})^2)'.format(self.power, 1.0 / (2.0 * self.variation), self.x0)

    def __unicode__ (self):
        return repr(self);