# unit tests and lint testss
# suitable for running as part of CI
# runs with a default config file, using mongo
pipenv run flake8 . --count --ignore=W293,E266,E302,E251,E225,E265,W291,E501,W503,C901 --show-source --max-complexity=200 --statistics

# startup the server
echo "Starting test findNeighbour servers to run tests with; waiting 15 seconds to ensure it has started  .."
nohup pipenv run python3 findNeighbour4_server.py config/default_test_config.json &
nohup pipenv run python3 findNeighbour4_server.py config/default_test_config_rdbms.json &
sleep 15 # wait for them to start

pipenv run pytest

# shut down test servers
kill -9 $(pgrep -f config/default_test_config.json)
kill -9 $(pgrep -f config/default_test_config_rdbms.json)
