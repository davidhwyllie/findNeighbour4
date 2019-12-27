# unit tests all scripts

# startup the server
echo "Starting test findNeighbour server to run tests with .."
nohup pipenv run python3 findNeighbour4_server.py &
echo $! > save_pid.txt      # get the PID
sleep 3 # wait for it to start

pipenv run python3 -m unittest discover -p="*.py"


# shut down the server
echo "Shutting down the temporary server for testing"
kill -9 `cat save_pid.txt`

rm save_pid.txt
