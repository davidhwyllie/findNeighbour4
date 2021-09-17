# wsgi.py
""" instatiates a findNeighbour4_server application

called by gunicorn to run parallel findNeighbour4 servers in wsgi container
example:
```
pipenv run gunicorn wsgi.py
```

Note: you cannot run unittests against the server using gunicorn using more than one worker using the default configuration file.
This is because each worker thread will start in 'debug' mode, deleting the contents of the database on each startup.
To unit test using gunicorn and the standard unit test suite, use
```
pipenv run gunicorn wsgi:app --workers 1 --bind 127.0.0.1:5020
```

To run findNeighbour4_server.py using the developmental (werkzeug, part of flask) server do 
```
pipenv run python3 findNeighbour4_server.py <config_file_name> 
```
"""

from findNeighbour4_server import create_app
app = create_app()
