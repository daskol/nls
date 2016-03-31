#!/usr/bin/env python2
#   bench.py
#   (c) Daniel Bershatsky, 2016
#   See LISENCE for details.

from __future__ import print_function
from argparse import ArgumentParser
from sys import path

path.append('..')

from nls.model import Problem
from nls.pumping import GaussianPumping2D

from numpy import polyfit, polyval, linspace
from scipy.io import savemat, loadmat
from seaborn import color_palette
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
    nodes = [50, 100, 200, 300, 400, 500, 700, 1000]
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

def visualize(filenames=('nodes.mat', 'iters.mat', 'orders.mat')):
    nodes_mat = loadmat(filenames[0])
    iters_mat = loadmat(filenames[1])
    orders_mat = loadmat(filenames[2])

    palette = color_palette('muted')

    # curve fitting
    nodes = linspace(min(nodes_mat['nodes'][0]), max(nodes_mat['nodes'][0]), 50)
    nodes_poly = polyfit(list(nodes_mat['nodes'][0]) + [0.0], list(nodes_mat['times'][0]) + [0.0], 2)

    iters = linspace(min(iters_mat['iters'][0]), max(iters_mat['iters'][0]), 50)
    iters_poly = polyfit(iters_mat['iters'][0], iters_mat['times'][0], 1)

    # plotting

    fig = figure(figsize=(17, 6))

    ax = fig.add_subplot(1, 3, 1)
    ax.plot(nodes_mat['nodes'][0], nodes_mat['times'][0], 'o-', label='exp')
    ax.plot(nodes, polyval(nodes_poly, nodes), label='fit')
    ax.set_ylim(0.0)
    ax.set_xlabel('nodes')
    ax.set_ylabel('time, s')
    ax.legend(loc='best')

    ax = fig.add_subplot(1, 3, 2)
    ax.plot(iters_mat['iters'][0], iters_mat['times'][0], 'o-', label='exp')
    ax.plot(iters, polyval(iters_poly, iters), label='fit')
    ax.set_xlabel('iterations')
    ax.set_ylabel('time, s')
    ax.legend(loc='best')

    ax = fig.add_subplot(1, 3, 3)
    ax.bar(orders_mat['orders'][0] + 0.2, orders_mat['times'][0], color=palette[0], label='200x200')
    ax.bar(orders_mat['orders'][0] + 1.0, orders_mat['times'][1], color=palette[1], label='400x400')
    ax.set_xlabel('approximation order')
    ax.set_ylabel('time, s')
    ax.set_xticks((4, 6, 8))
    ax.set_xticklabels((3, 5, 7))
    ax.legend(loc='best')
    
    fig.savefig('benchmark.png')
    show()

def main():
    parser = ArgumentParser(prog='benchmark.py', description='Test perfomance.')
    parser.add_argument('--only-bench', action='store_true')
    parser.add_argument('--only-show', action='store_true')

    args = parser.parse_args()

    if args.only_bench and args.only_show:
        parser.print_help()
        parser.exit()

    if args.only_bench:
        bench_num_nodes()
        bench_num_iters()
        bench_orders()

    if args.only_show:
        visualize()


if __name__ == '__main__':
    main()