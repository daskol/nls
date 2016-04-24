#!/usr/bin/env python2
#   workr.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import absolute_import, print_function
from argparse import ArgumentParser
from os import mkdir
from os.path import exists, join
from sys import path
from requests import Session
from json import loads

path.append('..')

from nls.model import Problem
from nls.pumping import GaussianPumping1D, GaussianPumping2D

import logging


PORT = 5000
HOST = 'localhost'
URL = 'http://{0}:{1}/api/{2}'


def register(session):
    response = session.post(URL.format('register'))
    return response.json().get('worker_id')

def unregister(session, worker_id):
    response = session.post(URL.format('unregister'), json=dict(worker_id=worker_id))
    return response.status_code == 202

def get_job(session, worker_id):
    response = session.get(URL.format('job'), json=dict(worker_id=worker_id))
    return {} if response.status_code == 204 else response.json()

def post_job(session, worker_id, job_id):
    response = session.post(URL.format('job'), json=dict(worker_id=worker_id, job_id=job_id))
    return response.status_code == 204

def solve_problem(job_id, desc, output_dir='output'):
    desc = loads(desc.encode('utf8'))

    model = Problem().model(
        model = '2d',
        dx = 2.0e-1,
        dt = 2.0e-3,
        t0 = 0.0,
        u0 = 0.1,
        order = 5,
        num_nodes = 200,
        num_iters = 2000,
        pumping = GaussianPumping2D(power=desc['power'], variation=5.0),
        original_params = {
            'R': 0.0242057488654,
            'gamma': 0.0242057488654,
            'g': 0.00162178517398,
            'tilde_g': 0.0169440242057,
            'gamma_R': 0.242057488654,
        },
        dimless_params = {
        })

    if not exists(output_dir):
        mkdir(output_dir)

    solution = model.solve()
    solution.report()
    solution.store(join(output_dir, 'solution.mat'))

    return True

def main():
    # Logging settigns
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

    # Initialization and configuration
    parser = ArgumentParser(prog='workr', description='Workr peeks job from Schdulr and runs it.')
    parser.add_argument('--host', '-H', default=HOST, help='set host name')
    parser.add_argument('--port', '-p', default=PORT, help='set port number')

    args = parser.parse_args()

    global URL
    URL = URL.format(args.host, args.port, '{0}')

    # Main pipeline
    session = Session()
    worker_id = register(session)

    if not worker_id:
        exit('Could not register!')

    logging.info('Workr:', worker_id)

    while True:
        job = get_job(session, worker_id);

        if 'job_id' not in job or 'description' not in job:
            break

        job_id = job['job_id']
        description = job['description']

        logging.info('Job:', job_id, description)

        solved = solve_problem(job_id, description)

        if not post_job(session, worker_id, job_id):
            logging.error('Could not post the job!', exc_info=False)
            break

    if not unregister(session, worker_id):
        exit('Could not unregister!')

    logging.info('Done')


if __name__ == '__main__':
    main()