#!/usr/bin/env python2
#   bench.py
#   (c) Daniel Bershatsky, 2016
#   See LISENCE for details.

from __future__ import print_function
from sys import path

path.append('..')

from nls.model import Problem
from nls.pumping import GaussianPumping2D

from scipy.io import savemat, loadmat
from matplotlib.pyplot import figure, show


def create_model(num_nodes=200, num_iters=2000, order=5):
    return Problem().model(
        model = '2d',
        dx = 2.0e-1,
        dt = 2.0e-3,
        t0 = 0.0,
        u0 = 0.1,
        order = order,
        num_nodes = num_nodes,
        num_iters = num_iters,
        pumping = GaussianPumping2D(power=15.0, variation=3.14),
        original_params = {
            'R': 0.0242057488654,
            'gamma': 0.0242057488654,
            'g': 0.00162178517398,
            'tilde_g': 0.0169440242057,
            'gamma_R': 0.242057488654
        })

def bench_num_nodes():
    nodes = [50, 100, 200, 500, 1000]
    times = []

    for i, num_nodes in enumerate(nodes):
        model = create_model(num_nodes=num_nodes)
        solution = model.solve()
        solution.report()
        times.append(solution.getElapsedTime())

    savemat('nodes.mat', {'nodes': nodes, 'times': times})

def bench_num_iters():
    iters = range(500, 5001, 500)
    times = []

    for i, num_iters in enumerate(iters):
        model = create_model(num_iters=num_iters)
        solution = model.solve()
        solution.report()
        times.append(solution.getElapsedTime())

    savemat('iters.mat', {'iters': iters, 'times': times})

def bench_orders():
    orders = [3, 5, 7]
    nodes = [200, 400]
    times = [[], []]

    for i, num_nodes in enumerate(nodes):
        for j, order in enumerate(orders):
            model = create_model(num_nodes=num_nodes, order=order)
            solution = model.solve()
            solution.report()
            times[i].append(solution.getElapsedTime())

    savemat('orders.mat', {'orders': orders, 'nodes': nodes, 'times': times})

def visualize(filenames=('before_nodes.mat', 'before_iters.mat', 'before_orders.mat',
                         'after_nodes.mat', 'after_iters.mat', 'after_orders.mat')):
    before_nodes_mat = loadmat(filenames[0])
    before_iters_mat = loadmat(filenames[1])
    before_orders_mat = loadmat(filenames[2])

    after_nodes_mat = loadmat(filenames[3])
    after_iters_mat = loadmat(filenames[4])
    after_orders_mat = loadmat(filenames[5])

    fig = figure()
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(before_nodes_mat['nodes'][0], before_nodes_mat['times'][0], 'x-')
    ax.plot(after_nodes_mat['nodes'][0], after_nodes_mat['times'][0], 'x-')
    ax = fig.add_subplot(1, 3, 2)
    ax.plot(before_iters_mat['iters'][0], before_iters_mat['times'][0], 'x-')
    ax.plot(after_iters_mat['iters'][0], after_iters_mat['times'][0], 'x-')
    ax = fig.add_subplot(1, 3, 3)
    ax.plot(after_orders_mat['orders'][0], after_orders_mat['times'][0], 'x-')
    ax.plot(before_orders_mat['orders'][0], before_orders_mat['times'][1], 'x-')
    show()

def main():
    bench_num_nodes()
    bench_num_iters()
    bench_orders()
    visualize()


if __name__ == '__main__':
    main()