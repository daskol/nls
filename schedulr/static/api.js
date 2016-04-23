//  api.js
//  (c) Daniel Bershatsky, 2016
//  See LICENSE for details

function cleanTable (table) {
    while (table.firstChild) {
        table.removeChild(table.firstChild);
    }
}

function makeTableHeader (parent, columns) {
    var row = document.createElement('tr');

    for (var column of columns) {
        var tag = document.createElement('th')
        tag.textContent = column;
        row.appendChild(tag)
    }

    parent.appendChild(row)
}

function makeTableBody (parent, columns, rows) {
    for (var row of rows) {
        var tag = document.createElement('tr');

        for (var column of columns) {
            var col = document.createElement('td');

            if (column == "description") {
                col.textContent = JSON.stringify(JSON.parse(row[column]), null, 2);
            } else {
                col.textContent = row[column];
            }

            tag.appendChild(col);
        }

        parent.appendChild(tag);
    }
}

function buildTable (element, header, rows) {
    cleanTable(element);
    makeTableHeader(element, header);
    makeTableBody(element, header, rows);
}

function updateActiveJobs () {
    $.get('/api/jobs', function (response) {
        var active_jobs = document.getElementById('active-jobs');
        var columns = response['columns'];
        var jobs = response['jobs'];

        buildTable(active_jobs, columns, jobs);
    });
}

function updateActiveWorkers () {
    $.get('/api/workers', function (response) {
        var active_workers = document.getElementById('active-workers');
        var columns = response['columns'];
        var workers = response['workers'];

        buildTable(active_workers, columns, workers);
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