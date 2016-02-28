#   nls/model.py
#   This module define core abstractions that maps to problem and model definition.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function

from pprint import pprint
from time import time
from types import FunctionType
from datetime import datetime
from numpy import array, exp, sqrt, arange, ones, zeros, meshgrid, mgrid, linspace, angle, gradient
from scipy.integrate import simps
from scipy.io import loadmat, savemat
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, cm
from matplotlib.pyplot import figure, plot, show, title, xlabel, ylabel, subplot, legend, xlim, ylim, contourf, hold, colorbar

from .animation import *
from .native import *
from .pumping import *
from .solver import *


class Problem(object):
    """Entry point in any computation. It implements design pattern `Factory` that used to construct object of type
    `Model`.
    """

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

        kwargs['model'] == 'default' if 'model' not in kwargs else kwargs['model']

        if 'model' in kwargs and kwargs['model'] in ('1d', 'default'):
            return self.fabricateModel1D(*args, **kwargs)
        elif 'model' in kwargs and kwargs['model'] == '2d':
            return self.fabricateModel2D(*args, **kwargs)
        else:
            raise Exception('Unknown model passed!')

    def fabricateModel1D(*args, **kwargs):
        kwargs['dx'] = 1.0e-1 if 'dx' not in kwargs else kwargs['dx']
        kwargs['dt'] = 1.0e-3 if 'dt' not in kwargs else kwargs['dt']
        kwargs['t0'] = 0.0e+0 if 't0' not in kwargs else kwargs['t0']
        kwargs['u0'] = 1.0e-1 if 'u0' not in kwargs else kwargs['u0']
        kwargs['order'] = 5 if 'order' not in kwargs else kwargs['order']
        kwargs['pumping'] = GaussianPumping() if 'pumping' not in kwargs else kwargs['pumping']
        kwargs['num_nodes'] = 1000 if 'num_nodes' not in kwargs else kwargs['num_nodes']
        kwargs['num_iters'] = 100000 if 'num_iters' not in kwargs else kwargs['num_iters']

        if type(kwargs['u0']) in (int, float, complex):
            kwargs['u0'] = kwargs['u0'] * ones(kwargs['num_nodes']) 
        elif isinstance(kwargs['u0'], FunctionType):
            grid = linspace(0.0, kwargs['dx'] * kwargs['num_nodes'], kwargs['num_nodes'])
            kwargs['u0'] = kwargs['u0'](grid)

        return Model1D(**kwargs)

    def fabricateModel2D(self, *args, **kwargs):
        kwargs['dx'] = 1.0e-1 if 'dx' not in kwargs else kwargs['dx']
        kwargs['dt'] = 1.0e-3 if 'dt' not in kwargs else kwargs['dt']
        kwargs['t0'] = 0.0e+0 if 't0' not in kwargs else kwargs['t0']
        kwargs['u0'] = 1.0e-1 if 'u0' not in kwargs else kwargs['u0']
        kwargs['order'] = 3 if 'order' not in kwargs else kwargs['order']
        kwargs['pumping'] = GaussianPumping() if 'pumping' not in kwargs else kwargs['pumping']
        kwargs['num_nodes'] = 40 if 'num_nodes' not in kwargs else kwargs['num_nodes']
        kwargs['num_iters'] = 1000 if 'num_iters' not in kwargs else kwargs['num_iters']

        if type(kwargs['u0']) in (int, float, complex):
            kwargs['u0'] = kwargs['u0'] * ones((kwargs['num_nodes'], kwargs['num_nodes'])) 

        return Model2D(**kwargs)


class AbstractModel(object):
    """Base type for objects which constructed with `Problem` class. Child object of this class implements computation
    and other related routines. This class defines common routines of initialization, solving, and model storage.
    """

    def __init__(self, *args, **kwargs):
        pprint({
            'dt': kwargs['dt'],
            'dx': kwargs['dx'],
            'order': kwargs['order'],
            'num_nodes': kwargs['num_nodes'],
            'num_iters': kwargs['num_iters'],
            'pumping': kwargs['pumping'],
            'originals': kwargs['original_params'],
            })
        self.solution = Solution(kwargs['dt'], kwargs['dx'], kwargs['num_nodes'], kwargs['order'], kwargs['num_iters'],
                             kwargs['pumping'], kwargs['original_params'], kwargs['u0'])
        self.solver = None

    def solve(self, num_iters=None):
        """Call solver that is aggregated certain child objects.
        """
        return self.solver(num_iters)

    def getChemicalPotential(self):
        """Call solver in order to calculate chemical potential.
        """
        self.mu = self.solver.chemicalPotential()
        return self.mu

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


class Model1D(AbstractModel):
    """Default model that is NLS equation with reservoir in axe symmentic case.
    """

    def __init__(self, *args, **kwargs):
        super(Model1D, self).__init__(*args, **kwargs)

        self.solver = Solver1D(self.solution)


class Model2D(AbstractModel):
    """Model that is NLS equation with reservoir on two dimensional grid.
    """

    def __init__(self, *args, **kwargs):
        super(Model2D, self).__init__(*args, **kwargs)

        self.solver = Solver2D(self.solution)


class Solution(object):
    """Object that represents solution of a given model. Also it contains all model parameters and has ability to store
    and to load solution. TODO: improve design.
    """

    def __init__(self, dt, dx, num_nodes, order, num_iters, pumping, originals, init_solution):
        self.dt = dt
        self.dx = dx
        self.order = order
        self.num_nodes = num_nodes
        self.num_iters = num_iters
        self.pumping = pumping
        self.init_sol = init_solution
        self.solution = None
        self.originals = originals
        self.coeffs = zeros(23)
        self.elapsed_time = 0.0

        hbar = 6.61e-34
        m_e = 9.1e-31
        m_0 = 1.0e-5 * m_e

        phi0 = sqrt(originals['gamma'] / (2.0 * originals['g']))
        t0 = phi0
        x0 = sqrt(hbar * t0 / (2 * m_0))
        n0 = 2.0 / (originals['R'] * t0)

        # NLS equation coeficients
        self.coeffs[0] = 1.0  # \partial_t
        self.coeffs[1] = 1.0  # \nabla^2
        self.coeffs[2] = 1.0  #
        self.coeffs[3] = 1.0  # linear damping
        self.coeffs[4] = 1.0  # originals['g'] * phi0 ** 3  # nonlinearity
        self.coeffs[5] = 4.0 * originals['tilde_g'] / originals['R']  #* phi0 * n0  # interaction to reservoir

        # Reservoir equation coefficients
        self.coeffs[10] = 0.0  # \parital_t
        self.coeffs[11] = 1.0 / (n0 * originals['gamma_R'])  # pumping coefficient
        self.coeffs[12] = 1.0  # damping
        self.coeffs[13] = originals['R'] * phi0 ** 2 / originals['gamma_R']  # interaction term
        self.coeffs[14] = 0.0  # diffusive term

        print(self.coeffs)

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
        if len(self.init_sol.shape) == 1:
            right = self.num_nodes * self.dx
            left = 0.0
            x = linspace(left, right, self.num_nodes)
            grid = meshgrid(x)
            return self.pumping(*grid)
        else:
            right = self.num_nodes * self.dx / 2
            left = -right
            x = linspace(left, right, self.num_nodes)
            grid = meshgrid(x, x)
            return self.pumping(*grid)

    def getCoefficients(self):
        return self.coeffs

    def getInitialSolution(self):
        return self.init_sol

    def getSolution(self):
        return self.solution

    def getElapsedTime(self):
        return self.elapsed_time

    def getParticleNumber(self, method='simps'):
        return simps((self.solution.conj() * self.solution).real, dx=self.dx)

    def setNumberOfIterations(self, num_iters):
        self.num_iters = num_iters

    def setPumping(self, pumping):
        self.pumping = pumping

    def setInitialSolution(self, solution):
        self.init_sol = solution

    def setSolution(self, solution):
        self.solution = solution

    def setElapsedTime(self, seconds):
        self.elapsed_time = seconds

    def visualize(self, *args, **kwargs):
        if len(self.init_sol.shape) == 1:
            self.visualize1d(*args, **kwargs)
        else:
            self.visualize2d(*args, **kwargs)

    def visualize1d(self, *args, **kwargs):
        x = arange(0.0, self.dx * self.num_nodes, self.dx)
        p = self.pumping(x)  # pumping profile
        u = (self.solution.conj() * self.solution).real  # density profile
        n = self.coeffs[11] *  p / (self.coeffs[12] + self.coeffs[13] * u)

        def rect_plot(subplot_number, value, label, name, labelx, labely, xmax=20):
            subplot(2, 3, subplot_number)
            hold(False)
            plot(x, value, label=label)
            xlim((0, xmax))
            legend(loc='best')
            title(name)
            xlabel(labelx)
            ylabel(labely)

        rect_plot(1, p, 'pumping', 'Pumping profile.', 'r', 'p')
        rect_plot(2, u, 'density', 'Density distribution of BEC.', 'r', 'u')
        rect_plot(3, n, 'reservoir', 'Density distribution of reservoir.', 'r', 'n')

        def polar_plot(subplot_number, value, xmax=20):
            hold(False)
            subplot(2, 3, subplot_number, polar=True)
            theta = arange(0, 2 * 3.14 + 0.1, 0.1)
            contourf(theta, x, array([value for _ in theta]).T)
            ylim((0, xmax))

        polar_plot(4, p)
        polar_plot(5, u)
        polar_plot(6, n)

    def visualize2d(self, *args, **kwargs):
        right = self.num_nodes * self.dx / 2
        left = -right
        x = linspace(left, right, self.num_nodes)
        gx, gy = meshgrid(x, x)
        p = self.getPumping()
        u = (self.solution.conj() * self.solution).real  # density profile
        n = self.coeffs[11] *  p / (self.coeffs[12] + self.coeffs[13] * u)

        fig = kwargs['figure'] if 'figure' in kwargs else figure()

        def surface_plot(subplot_number, value, label, name, labels):
            ax = fig.add_subplot(130 + subplot_number, projection='3d')
            ax.plot_surface(gx, gy, value, label=label)
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            ax.set_zlabel(labels[2])
            ax.set_title(name)

        def contour_plot(subplot_number, value, label, name, labels):
            levels = linspace(0.0, value.max() + 1.0e-3, 11)
            extent = (gx[0, 0], gx[-1, -1], gy[0, 0], gy[-1, -1])

            ax = fig.add_subplot(130 + subplot_number, aspect='equal')
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            ax.set_title(name)
            
            cp = ax.contourf(gx, gy, value, levels, cmap=cm.get_cmap('Accent'), extent=extent)
            colorbar(cp, orientation='horizontal')

        def stream_plot(subplot_number, value, label, name, labels):
            """Plot stream of complex function.
            :param: value tuple Pair of absolute value and its angle.
            """
            jx, jy = value[0] * gradient(value[1])

            ax = fig.add_subplot(120 + subplot_number, aspect='equal')
            ax.streamplot(gx, gy, jx, jy, color=value[0])
            ax.set_xlim(gx[0, 0], gx[-1, -1])
            ax.set_ylim(gy[0, 0], gy[-1, -1])
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            ax.set_title(name)

        def density_plot(subplot_number, value, label, name, labels):
            extent = (gx[0, 0], gx[-1, -1], gy[0, 0], gy[-1, -1])

            ax = fig.add_subplot(120 + subplot_number, aspect='equal')
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            ax.set_title(name)
            ax.imshow(value[0], extent=extent)
            ax.contour(gx, gy, value[1].real, [0.0], colors='red', extent=extent)
            ax.contour(gx, gy, value[1].imag, [0.0], colors='blue', extent=extent)
        
        if 'stream' in kwargs and kwargs['stream']:
            stream_plot(1, (u, angle(self.solution)), 'phase gradient', 'Condensate streams', ('x', 'y'))
            density_plot(2, (u, self.solution), 'density', 'Density distribution of BEC.', ('x', 'y'))
        else:
            helper_plot = contour_plot if 'contour' in kwargs and kwargs['contour'] else surface_plot

            helper_plot(1, p, 'pumping', 'Pumping profile.', ('x', 'y', 'p'))
            helper_plot(2, u, 'density', 'Density distribution of BEC.', ('x', 'y', 'u'))
            helper_plot(3, n, 'reservoir', 'Density distribution of reservoir.', ('x', 'y', 'n'))

        if 'filename' in kwargs:
            fig.savefig(kwargs['filename'])

    def show(self):
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