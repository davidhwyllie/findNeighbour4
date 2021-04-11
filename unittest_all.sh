# unit tests and lint testss
# suitable for running as part of CI

pipenv run flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics

# F401 important
pipenv run flake8 . --count --ignore=W293,E266,E302,E251,E225,E265,W291 --exit-zero --max-complexity=10 --max-line-length=127 --statistics

# startup the server
echo "Starting test findNeighbour server to run tests with; waiting 15 seconds to ensure it has started  .."
nohup pipenv run python3 findNeighbour4_server.py &
sleep 15 # wait for it to start

pipenv run pytest test
