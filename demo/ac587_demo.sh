# start the server
# make sure you have installed prerequisites 
# pipenv install . --skip-lock

./fn4_startup.sh demos/AC587/config/config.json

# load data
# doesn't find library fn4client if in demo/ [TODO]
pipenv run python3 demo_ac587.py 

