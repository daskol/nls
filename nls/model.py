#   nls/model.py
#   This module define core abstractions that maps to problem and model definition.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import absolute_import, print_function

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
from .version import *


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
        if 'filename' in kwargs:
            return self.modelFromFile(kwargs['filename'])

        if 'params' in kwargs:
            params = kwargs.pop('params')

        kwargs['model'] = 'default' if 'model' not in kwargs else kwargs['model']
        kwargs['original_params'] = {} if 'original_params' not in kwargs else kwargs['original_params']

        if 'R' not in kwargs['original_params']:
            kwargs['original_params']['R'] = 0.0242057488654
        if 'gamma' not in kwargs['original_params']:
            kwargs['original_params']['gamma'] = 0.0242057488654
        if 'g' not in kwargs['original_params']:
            kwargs['original_params']['g'] = 0.00162178517398
        if 'tilde_g' not in kwargs['original_params']:
            kwargs['original_params']['tilde_g'] = 0.0169440242057
        if 'gamma_R' not in kwargs['original_params']:
            kwargs['original_params']['gamma_R'] = 0.242057488654

        if kwargs.get('model') in ('1d', 'default', str(Model1D)):
            return self.fabricateModel1D(*args, **kwargs)
        elif kwargs.get('model') in ('2d', str(Model1D)):
            return self.fabricateModel2D(*args, **kwargs)
        else:
            raise Exception('Unknown model passed!')

    def modelFromFile(self, filename):
        def modelFromFileLikeObject(filename):
            mat = loadmat(filename)
            if 'model' in mat:
                return self.model(model=mat['model'][0]).restore(filename)

        if isinstance(filename, file):
            return modelFromFileLikeObject(filename)
        else:
            with open(filename) as f:
                return modelFromFileLikeObject(filename)

    def fabricateModel1D(self, *args, **kwargs):
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
        self.dt = kwargs['dt']
        self.dx = kwargs['dx']
        self.order = kwargs['order']
        self.num_nodes = kwargs['num_nodes']
        self.num_iters = kwargs['num_iters']
        self.pumping = kwargs['pumping']
        self.init_sol = kwargs['u0']
        self.originals = kwargs['original_params']
        self.coeffs = zeros(23)
        self.verbose = bool(kwargs.get('verbose'))
        self.solver = None

        hbar = 6.61e-34
        m_e = 9.1e-31
        m_0 = 1.0e-5 * m_e

        phi0 = sqrt(self.originals['gamma'] / (2.0 * self.originals['g']))
        t0 = phi0
        x0 = sqrt(hbar * t0 / (2 * m_0))
        n0 = 2.0 / (self.originals['R'] * t0)

        # NLS equation coeficients
        self.coeffs[0] = 1.0  # \partial_t
        self.coeffs[1] = 1.0  # \nabla^2
        self.coeffs[2] = 1.0  #
        self.coeffs[3] = 1.0  # linear damping
        self.coeffs[4] = 1.0  # self.originals['g'] * phi0 ** 3  # nonlinearity
        self.coeffs[5] = 4.0 * self.originals['tilde_g'] / self.originals['R']  #* phi0 * n0  # interaction to reservoir

        # Reservoir equation coefficients
        self.coeffs[10] = 0.0  # \parital_t
        self.coeffs[11] = 1.0 / (n0 * self.originals['gamma_R'])  # pumping coefficient
        self.coeffs[12] = 1.0  # damping
        self.coeffs[13] = self.originals['R'] * phi0 ** 2 / self.originals['gamma_R']  # interaction term
        self.coeffs[14] = 0.0  # diffusive term

    def __repr__(self):
        from pprint import pformat
        return pformat({
            'dt': self.dt,
            'dx': self.dx,
            'order': self.order,
            'num_nodes': self.num_nodes,
            'num_iters': self.num_iters,
            'pumping': self.pumping,
            'originals': self.originals,
            }) + '\n' + str(self.coeffs)

    def getApproximationOrder(self):
        return self.order

    def getCharacteristicScale(self, scale):
        hbar = 6.61e-34
        m_e = 9.1e-31
        m_0 = 1.0e-5 * m_e

        phi0 = sqrt(self.originals['gamma'] / (2.0 * self.originals['g']))
        t0 = phi0
        x0 = sqrt(hbar * t0 / (2 * m_0))
        n0 = 2.0 / (self.originals['R'] * t0)

        scales = {
            'x': x0,
            't': t0,
            'n': n0,
            'phi': phi0,
        }

        return scales[scale] if scale in scales else None

    def getChemicalPotential(self):
        """Call solver in order to calculate chemical potential.
        """
        self.mu = self.solver.chemicalPotential()
        return self.mu

    def getCoefficients(self):
        return self.coeffs

    def getInitialSolution(self):
        return self.init_sol

    def getModel(self):
        return self.model

    def getNumberOfIterations(self):
        return self.num_iters

    def getNumberOfNodes(self):
        return self.num_nodes

    def getParticleNumber(self, method='simps'):
        return simps((self.solution.conj() * self.solution).real, dx=self.dx)

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

    def getSpatialStep(self):
        return self.dx

    def getSolver(self):
        return self.solver

    def getTimeStep(self):
        return self.dt

    def setNumberOfIterations(self, num_iters):
        self.num_iters = num_iters

    def setPumping(self, pumping):
        self.pumping = pumping

    def setInitialSolution(self, solution):
        self.init_sol = solution

    def solve(self, num_iters=None):
        """Call solver that is aggregated certain child objects.
        """
        return self.solver(num_iters)

    def store(self, filename=None, label=None, desc=None, date=None):
        """Store object to mat-file. TODO: determine format specification
        """
        date = date if date else datetime.now()
        date = date.isoformat()
        filename = filename if filename else date + '.mat'

        matfile = {
            'model': str(type(self)),
            'date': date,
            'dim': len(self.init_sol.shape),
            'dimlesses': self.coeffs,
            'init_solution': self.init_sol,
            'num_iters': self.num_iters,
            'num_nodes': self.num_nodes,
            'order': self.order,
            'originals': self.originals,
            'pumping': self.getPumping(),
            'spatial_step': self.dx,
            'time_step': self.dt,
        }

        if desc:
            matfile['desc'] = desc
        if label:
            matfile['label'] = label

        savemat(filename, matfile)

    def restore(self, filename):
        """Restore object from mat-file. TODO: determine format specification
        """
        matfile = loadmat(filename)
        matfile['originals'] = matfile['originals'][0, 0]

        if matfile['dim'] == 1:
            matfile['init_solution'] = matfile['init_solution'][0, :]
            matfile['pumping'] = matfile['pumping'][0, :]

        self.coeffs = matfile['dimlesses'][0, :]
        self.init_sol = matfile['init_solution']
        self.num_nodes = matfile['num_nodes'][0, 0]
        self.num_iters = matfile['num_iters'][0, 0]
        self.pumping = GridPumping(matfile['pumping'])

        types = matfile['originals'].dtype
        values = matfile['originals']

        self.originals = dict(zip(types.names, (value[0, 0] for value in values)))

        if 'desc' in matfile:
            self.desc = str(matfile['desc'][0])
        if 'label' in matfile:
            self.label = str(matfile['label'][0])

        return self


class Model1D(AbstractModel):
    """Default model that is NLS equation with reservoir in axe symmentic case.
    """

    def __init__(self, *args, **kwargs):
        super(Model1D, self).__init__(*args, **kwargs)

        self.solver = Solver1D(self)


class Model2D(AbstractModel):
    """Model that is NLS equation with reservoir on two dimensional grid.
    """

    def __init__(self, *args, **kwargs):
        super(Model2D, self).__init__(*args, **kwargs)

        self.solver = Solver2D(self)


class Solution(object):
    """Object that represents solution of a given model. Also it contains all model parameters and has ability to store
    and to load solution.

    TODO: improve design.
    """

    def __init__(self, model, solution=None, verbose=False):
        self.elapsed_time = 0.0
        self.model = model
        self.solution = solution
        self.verbose = verbose

    def getElapsedTime(self):
        return self.elapsed_time

    def getSolution(self):
        return self.solution

    def setElapsedTime(self, seconds):
        self.elapsed_time = seconds

    def setSolution(self, solution):
        self.solution = solution

    def visualize(self, *args, **kwargs):
        if len(self.model.init_sol.shape) == 1:
            self.visualize1d(*args, **kwargs)
        else:
            self.visualize2d(*args, **kwargs)

    def visualize1d(self, *args, **kwargs):
        x = arange(0.0, self.model.dx * self.model.num_nodes, self.model.dx)
        p = self.model.pumping(x)  # pumping profile
        u = (self.solution.conj() * self.solution).real  # density profile
        n = self.model.coeffs[11] *  p / (self.model.coeffs[12] + self.model.coeffs[13] * u)

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
        right = self.model.num_nodes * self.model.dx / 2
        left = -right
        x = linspace(left, right, self.model.num_nodes)
        gx, gy = meshgrid(x, x)
        p = self.model.getPumping()
        u = (self.solution.conj() * self.solution).real  # density profile
        n = self.model.coeffs[11] *  p / (self.model.coeffs[12] + self.model.coeffs[13] * u)

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

        if kwargs.get('filename'):
            fig.savefig(kwargs['filename'])

    def show(self):
        show()

    def store(self, filename=None, label=None, desc=None, date=None):
        """Store object to mat-file. TODO: determine format specification
        """
        date = datetime.now() if date is None else date
        filename = filename if filename else date.isoformat() + '.mat'

        def storeWithFileLikeObject(file_like):
            content = {
                'elapsed_time': self.elapsed_time,
                'solution': self.solution,
                'version': version(),
            }
            self.model.store(file_like, label, desc, date)
            savemat(file_like, content, appendmat=True)

        if isinstance(filename, file):
            storeWithFileLikeObject(filename)
        else:
            with open(filename, 'wb') as f:
                storeWithFileLikeObject(f)

    def restore(self, filename):
        """Restore object from mat-file. TODO: determine format specification
        """
        matfile = loadmat(filename)

        if matfile['dim'] == 1:
            matfile['solution'] = matfile['solution'][0, :]

        self.elapsed_time = matfile['elapsed_time'][0, 0]
        self.solution = matfile['solution']

        return self

    def report(self):
        message = 'Elapsed in {0} seconds with {1} iteration on {2} grid nodes.'
        print(message.format(self.elapsed_time, self.model.getNumberOfIterations(), self.model.getNumberOfNodes()))


from .solver import *  # cyclic import fix