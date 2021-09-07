# unit tests and lint testss
# suitable for running as part of CI
# runs with a default config file, using mongo
pipenv run flake8 . --count --ignore=W293,E266,E302,E251,E225,E265,W291,E501,W503,C901 --show-source --max-complexity=200 --statistics

# shut down any running test servers
echo "Terminating any running test catwalk  processes"
for pid in $(pgrep -f CatWalk-PORT-599); do 
echo "$pid"
kill -9 $pid
done
echo "Terminating any running findneighbour4 processes"
for pid in $(pgrep -f config/default_test_config.json); do 
echo "$pid"
kill -9 $pid
done
for pid in $(pgrep -f config/default_test_config_rdbms.json); do 
echo "Terminating $pid"
kill -9 $pid
done

# startup the server
echo "Starting test findNeighbour servers to run tests with; waiting 15 seconds to ensure it has started  .."
echo "starting mongodb server with gunicorn and 4 workers"
rm test_startup.sh -f
pipenv run python3 configure.py config/default_test_config.json --prepare --n_workers 4 > test_startup.sh
chmod +x test_startup.sh
./test_startup.sh

#echo "starting rdbms server with werkzeug"
#nohup pipenv run python3 findNeighbour4_server.py config/default_test_config_rdbms.json &
sleep 15 # wait for them to start

pipenv run pytest test
#pipenv run python3 -m unittest test/test_server_rdbms.py

# shut down any running test servers
echo "Terminating any running findneighbour4 processes"
for pid in $(pgrep -f config/default_test_config.json); do 
echo "Terminating $pid"
kill -9 $pid
done
for pid in $(pgrep -f config/default_test_config_rdbms.json); do 
echo "Terminating $pid"
kill -9 $pid
done
echo "Terminating any running test catwalk  processes"
for pid in $(pgrep -f CatWalk-PORT-599); do 
echo "$pid"
kill -9 $pid
done

# remove test script
rm test_startup.sh

