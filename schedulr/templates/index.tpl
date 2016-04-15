<!doctype html>
<title>Schedulr</title>
<link rel="stylesheet" type="text/css" href="{{ url_for('static', filename='style.css') }}" media="screen" />
<script type="text/javascript" src="{{ url_for('static', filename='jquery.js') }}"></script>
<script type="text/javascript" src="{{ url_for('static', filename='api.js') }}"></script>
<div class="page" onload="updateOnLoad()">
  <h1>Schedulr</h1>

  <h2>Active Workers</h2>
  <input type="button" value="Update" onclick="updateActiveWorkers()">
  <table id="active-workers" border="1px">
  </table>

  <h2>Active Jobs</h2>
  <input type="button" value="Update" onclick="updateActiveJobs()">
  <table id="active-jobs" border="1px">
  </table>

  <h2>Uploads Jobs</h2>
  <form>
    <input type="button" value="Upload" onclick="uploadJobs()">
    <div></div>
    <textarea id='jobs' rows="20" cols="120"></textarea>
  </form>
</div>
