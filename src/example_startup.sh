# shell script illustrating testing of restful endpoint.
# note that this option uses the built in flask webserver.
# this is not the recommended means of production
# http://flask.pocoo.org/docs/0.12/deploying/mod_wsgi/#configuring-apache

# from /src directory
python3 webservice-server.py ../config/default_config.json        &   # uses default configuration; runs on localhost port 8184.
python3 webservice-server-rest.py ../config/default_config.jsom   & # runs on 5000
