# unit tests and lint testss
# suitable for running as part of CI

# optionally: --show-source
pipenv run flake8 . --count --select=E9,F63,F7,F8,F401,E713,W605,F811,F841,W191,E101,F522,F521,F712,E703,E714,F601,E711,F523,F524  --statistics > flake8_severe.txt

cat flake8_severe.txt

# F401 important
pipenv run flake8 . --count --ignore=W293,E266,E302,E251,E225,E265,W291 --exit-zero --max-complexity=10 --max-line-length=127 --statistics  > flake8_less_severe.txt

cat flake8_less_severe.txt

# startup the server
echo "Starting test findNeighbour server to run tests with; waiting 15 seconds to ensure it has started  .."
nohup pipenv run python3 findNeighbour4_server.py &
sleep 15 # wait for it to start

pipenv run coverage test

pipenv run coverage.py report

