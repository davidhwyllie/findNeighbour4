# start the server
pipenv run python3 configure.py demos/AC587/config/config.json --prepare --n_workers 10
nohup pipenv run gunicorn wsgi:app --workers 10 --bind 127.0.0.1:5032 --log-level info &

# load data
pipenv run demo/covidsim_load.py 

