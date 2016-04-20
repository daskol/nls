#!/usr/bin/env python2
#   schedulr.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function, absolute_import
from os import path
from sqlite3 import connect, Row
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify


app = Flask(__name__)
app.config.update(dict(
    DATABASE=path.join(app.root_path, 'schedulr.db'),
    DEBUG=True,
    SECRET_KEY='development key',
    USERNAME='admin',
    PASSWORD='default',
))
app.config.from_envvar('SCHEDULR_SETTINGS', silent=True)


def db_connection():
    connection = connect(app.config['DATABASE'])
    connection.row_factory = Row
    return connection

def db_get():
    if not hasattr(g, 'db'):
        g.db = db_connection()
    return g.db

def update_worker(worker_id, status):
    db = db_get()
    rv = db.execute("UPDATE workers SET stage = ?, last_seen = DATETIME('now','localtime') WHERE id = ?;", (status, worker_id))
    db.commit()

@app.teardown_appcontext
def db_close(error):
    if hasattr(g, 'db'):
        g.db.close()

@app.route('/')
def root():
    return render_template('index.tpl')

@app.route('/api/register', methods=['POST'])
def register():
    db = db_get()
    rv = db.execute("INSERT INTO workers VALUES(NULL, 1, 'registration', DATETIME('now','localtime'));")
    db.commit()
    return jsonify(worker_id=rv.lastrowid), 201

@app.route('/api/unregister', methods=['POST'])
def unregister():
    json = request.get_json()

    if not json:
        return 'Multiform part data is not application/json.', 400

    if 'worker_id' not in json:
        return 'Field `worker_id` not in request data.', 400

    db = db_get()
    rv = db.execute("UPDATE workers SET is_active = 0, stage = ? WHERE id = ?;", ('invalidation', json['worker_id'],))
    db.commit()

    return '', 202

@app.route('/api/job', methods=['GET', 'POST'])
def job():
    if request.method == 'GET':
        return get_job()
    elif request.method == 'POST':
        return post_job()

def get_job():
    json = request.get_json()

    if not json:
        return 'Multiform part data is not application/json.', 400

    if 'worker_id' not in json:
        return 'Field `worker_id` not in request data.', 400

    update_worker(json['worker_id'], 'get_job')

    db = db_get()
    rv = db.execute("UPDATE jobs SET worker_id = ?, is_active = 1, started_at = DATETIME('now','localtime') WHERE id in (SELECT id FROM jobs WHERE is_active = 0 AND worker_id IS NULL LIMIT 1);", (json['worker_id'],))
    db.commit()

    rv = db.execute('SELECT id FROM jobs WHERE worker_id = ? AND is_active = 1 ORDER BY started_at DESC LIMIT 1;', (json['worker_id'],))

    row = rv.fetchone()

    if not row:
        return jsonify(), 204

    job_id = row['id']

    rv = db.execute('SELECT description FROM jobs WHERE id = ?;', (job_id,))
    values = rv.fetchone()

    if values:
        job = dict(job_id=job_id, description=values['description'])
        return jsonify(job), 200
    else:
        return jsonify(), 204

def post_job():
    json = request.get_json()

    if not json:
        return 'Multiform part data is not application/json.', 400

    if 'worker_id' not in json:
        return 'Field `worker_id` not in request data.', 400

    update_worker(json['worker_id'], 'post_job')

    db = db_get()
    rv = db.execute("UPDATE jobs SET is_active = 0, is_done = 1, finished_at = DATETIME('now','localtime') WHERE id = ? AND worker_id = ?;", (json['job_id'], json['worker_id'],))
    db.commit()

    return '', 204

@app.route('/api/jobs', methods=['GET', 'POST'])
def jobs():
    if request.method == 'GET':
        return get_jobs()
    elif request.method == 'POST':
        return post_jobs()

def get_jobs():
    result = db_get().execute('SELECT * FROM jobs WHERE is_active = 1 AND is_done = 0;')
    rows = [dict(zip(row.keys(), row)) for row in result.fetchall()]
    return jsonify(jobs=rows)

def post_jobs():
    json = request.get_json(force=True)

    db = db_get()
    rv = db.executemany('INSERT INTO jobs VALUES(NULL, NULL, ?, 0, 0, NULL, NULL);', zip(json['jobs']))
    db.commit()

    return '', 201

@app.route('/api/workers', methods=['GET'])
def workers():
    result = db_get().execute('SELECT * FROM workers;')
    rows = [dict(zip(row.keys(), row)) for row in result.fetchall()]
    return jsonify(workers=rows)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)