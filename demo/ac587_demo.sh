pipenv run python3 configure.py demos/AC587/config/config.json --prepare --n_workers 2
nohup pipenv run gunicorn wsgi:app --workers 2 --bind 127.0.0.1:5032 --log-level info &

