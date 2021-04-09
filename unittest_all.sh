# unit tests all scripts
# suitable for running as part of CI, or initial startup


#pipenv run python3 unittest_core.py 

# startup the server
echo "Starting test findNeighbour server to run tests with; waiting 15 seconds to ensure it has started  .."
nohup pipenv run python3 findNeighbour4_server.py &
echo $! > save_pid.txt      # get the PID
sleep 15 # wait for it to start

# test client
echo "Running temporary test server as PID = " `cat save_pid.txt`

#pipenv run python3 unittest_core.py
pipenv run python3 unittest_core.py

# shut down the server
echo "Shutting down the temporary server for testing"
kill -9 `cat save_pid.txt`

rm save_pid.txt
