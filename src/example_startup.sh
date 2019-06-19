# shell script illustrating testing of restful endpoint.
# note that this option uses the built in flask webserver.
# this is not the recommended means of production
# http://flask.pocoo.org/docs/0.12/deploying/mod_wsgi/#configuring-apache

# from /src directory
python3 findNeighbour3-server.py ../config/default_config.json   & # runs on 5000
