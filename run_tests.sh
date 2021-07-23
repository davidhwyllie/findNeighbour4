# unit tests and lint testss
# suitable for running as part of CI

pipenv run flake8 . --count --ignore=W293,E266,E302,E251,E225,E265,W291,E501,W503,C901 --show-source --max-complexity=200 --statistics

# startup the server
echo "Starting test findNeighbour server to run tests with; waiting 15 seconds to ensure it has started  .."
nohup pipenv run python3 findNeighbour4_server.py &
sleep 15 # wait for it to start

pipenv run pytest
