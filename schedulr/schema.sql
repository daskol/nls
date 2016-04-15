DROP TABLE IF EXISTS jobs;
CREATE TABLE jobs (
    id INTEGER  PRIMARY KEY AUTOINCREMENT,
    worker_id INTEGER,
    description TEXT,
    is_active BOOLEAN,
    is_done BOOLEAN,
    started_at DATETIME,
    finished_at DATETIME
);

DROP TABLE IF EXISTS workers;
CREATE TABLE workers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    is_active BOOLEAN,
    stage TEXT,
    last_seen DATETIME
);
