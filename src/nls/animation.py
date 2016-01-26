#   nls/animation.py
#   This module define routines and type to record evolution of solution and store it into file.
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function
from time import time
from matplotlib import animation
from matplotlib.pyplot import figure
from .pumping import GaussianPumping


class AbstractAnimation(object):
    """Animation base class that contains common method of initialization and rendering video.
    """

    def __init__(self, model, frames, step=1):
        self.elapsed_time = 0.0
        self.model = model
        self.frames = frames
        self.step = step
        self.writer = animation.writers['ffmpeg'](fps=15, metadata={'title': 'Exciton-polariton condensation.'})

    def getElapsedTime(self):
        return self.elapsed_time

    def animate(self, filename):
        self.elapsed_time = -time()
        dpi = 100
        fig = figure(figsize=(16, 9), dpi=dpi)
        with self.writer.saving(fig, filename, dpi):
            pass
        self.elapsed_time += time()

    def report(self):
        message = 'Elapsed in {0} seconds with {1} frames and {2} step.'
        print(message.format(self.elapsed_time, self.frames, self.step))


class IterationIncreaseAnimation(AbstractAnimation):

    def __init__(self, model, frames, step=1):
        super(IterationIncreaseAnimation, self).__init__(model, frames, step)

    def animate(self, filename):
        self.elapsed_time = -time()
        dpi = 100
        fig = figure(figsize=(16, 9), dpi=dpi)
        with self.writer.saving(fig, filename, dpi):
            for i in xrange(self.frames + 1):
                solution = self.model.solve(self.step)  # Fix references entanglement
                solution.setInitialSolution(solution.getSolution())
                solution.visualize()
                self.writer.grab_frame()
        self.elapsed_time += time()


class PumpingRadiusIncreaseAnimation(AbstractAnimation):

    def __init__(self, model, frames, step=1):
        super(PumpingRadiusIncreaseAnimation, self).__init__(model, frames, step)

    def animate(self, filename):
        self.elapsed_time = -time()
        dpi = 100
        fig = figure(figsize=(16, 9), dpi=dpi)
        with self.writer.saving(fig, filename, dpi):
            for i in xrange(self.frames + 1):
                origin = i * self.step
                pumping = GaussianPumping(power=3.0, x0=+origin, variation=6.84931506849) \
                       + GaussianPumping(power=3.0, x0=-origin, variation=6.84931506849)
                self.model.solution.setPumping(pumping)
                solution = self.model.solve()  # Fix references entanglement
                solution.setInitialSolution(solution.getSolution())
                solution.visualize()
                self.writer.grab_frame()
        self.elapsed_time += time()