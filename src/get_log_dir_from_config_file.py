# reads directory of the log file specified in a config file
# expects the config file to be passed as a single argument.

# if the config file cannot be found, or there is no LOGFILE key in the config json,
# returns nothing.
import os
import json
import sys

if __name__ == '__main__':

    if not os.path.exists(sys.argv[1]):
        exit(0)
    
    with open(sys.argv[1]) as jdata:
        data = json.load(jdata)
        try:
            logfile = data['LOGFILE']
        except KeyError:
            exit(0);

        print(os.path.abspath(os.path.dirname(logfile))+"/")
        
