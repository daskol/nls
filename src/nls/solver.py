#   nls/solver.py
#   This file contains solver abstraction layer that hides relevant implimentation.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from .native import nls
from .model import *


class AbstractSolver(object):
    """Wrapper(or Facade template) of native fortran subroutine `solve_nls`.
    """
    
    def __init__(self, solution):
        self.solution = solution
        self.elapsed_time = 0.0

    def __call__(self, num_iters=None):
        if num_iters:
            self.solution.setNumberOfIterations(num_iters)

        self.elapsed_time = -time()
        self.solution.setSolution(self.solve(
            self.solution.getTimeStep(),
            self.solution.getSpatialStep(),
            self.solution.getApproximationOrder(),
            self.solution.getNumberOfIterations(),
            self.solution.getPumping(),
            self.solution.getCoefficients(),
            self.solution.getInitialSolution()))
        self.elapsed_time += time()
        self.solution.setElapsedTime(self.elapsed_time)

        return self.solution

    def solve(*argv, **kwargs):
        """This method should be override in child classes that specify solver routine(template method pattern).
        """
        raise Exception('AbstractSolver: native solver routine is not passed!')


class Solver1D(AbstractSolver):

    def __init__(self, solution):
        super(AbstractSolver, self).__init__()

        self.solution = solution
        self.elapsed_time = 0.0


    def solve(self, *args, **kwargs):
        return nls.solve_nls(*args, **kwargs)


class Solver2D(AbstractSolver):

    def __init__(self, solution):
        super(AbstractSolver, self).__init__()

        self.solution = solution
        self.elapsed_time = 0.0

    def solve(self, *args, **kwargs):
        return nls.solve_nls_2d(*args, **kwargs)