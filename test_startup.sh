nohup pipenv run gunicorn wsgi:app --workers 4 --bind 127.0.0.1:5020 --log-level info & 
