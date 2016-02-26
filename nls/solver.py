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

    def solve(*args, **kwargs):
        """This method should be override in child classes that specify solver routine(template method pattern).
        """
        raise Exception('AbstractSolver: native solver routine is not passed!')

    def chemicalPotential(self, *args, **kwargs):
        return self.chemicalPotentialRoutine(
            self.solution.getSpatialStep(),
            self.solution.getPumping(),
            self.solution.getCoefficients(),
            self.solution.getSolution())

    def chemicalPotentialRoutine(self, *args, **kwargs):
        """This method should be override in child classses that specify chemical potencial calculation routine.
        """
        raise Exception('AbstractSolver: native solver routine is not passed!')


class Solver1D(AbstractSolver):
    """One dimensional solver that call native Fortran routine that solves NLS equation in axial symmentry in 2D.
    """

    def __init__(self, solution):
        super(Solver1D, self).__init__(solution)

    def solve(self, *args, **kwargs):
        return nls.solve_nls(*args, **kwargs)

    def chemicalPotentialRoutine(self, *args, **kwargs):
        return nls.chemical_potential_1d(*args, **kwargs)


class Solver2D(AbstractSolver):
    """One dimensional solver that call native Fortran routine that solves NLS equation on a squared grid.
    """


    def __init__(self, solution):
        super(Solver2D, self).__init__(solution)

    def solve(self, *args, **kwargs):
        return nls.solve_nls_2d(*args, **kwargs)

    def chemicalPotentialRoutine(*args, **kwargs):
        return nls.chemical_potential_2d(*args, **kwargs)