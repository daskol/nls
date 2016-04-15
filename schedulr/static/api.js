//  api.js
//  (c) Daniel Bershatsky, 2016
//  See LICENSE for details

function cleanTable (table) {
    while (table.firstChild) {
        table.removeChild(table.firstChild);
    }
}

function updateActiveJobs () {
    $.get('/api/jobs', function (response) {
        var active_jobs = document.getElementById('active-jobs');
        var jobs = response['jobs'];

        cleanTable(active_jobs);

        for (var job of jobs) {
            var row = document.createElement('tr');

            for (var attr in job) {
                var col = document.createElement('td');
                col.textContent = job[attr];
                row.appendChild(col);
            }

            active_jobs.appendChild(row);
        }
    });
}

function updateActiveWorkers () {
    $.get('/api/workers', function (response) {
        var active_workers = document.getElementById('active-workers');
        var workers = response['workers'];

        cleanTable(active_workers);

        for (var worker of workers) {
            var row = document.createElement('tr');

            for (var attr in worker) {
                var col = document.createElement('td');
                col.textContent = worker[attr];
                row.appendChild(col);
            }

            active_workers.appendChild(row);
        }
    });
}

function uploadJobs () {
    var jobs = document.getElementById('jobs');
    var list = [];

    try {
        var json_jobs = JSON.parse(jobs.value);

        for (var job of json_jobs) {
            list.push(JSON.stringify(job));
        }

        var content = JSON.stringify({
            jobs: list
        });

        $.post('/api/jobs', content, null, "json");
    } catch (e) {
        alert('Job description is not a JSON serializable!');
    }
}

$(function () {
    updateActiveJobs();
    updateActiveWorkers();
})