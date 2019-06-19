
# The key dependencies of the main server are as follows:  
Python 3.5+ (tested with 3.5 and 3.7)  
Mongodb  (tested with 3.6.1. and 4.0).  See [here](mongoinstall.md).

The below method installs packages necessary for execution of all parts of the software in a virtual environment.
In you have not already done so,  pipenv.
```pip3 install pipenv```  

If you do not have root/admin access, you can install this locally:  
```pip3 install pipenv --user```

Create a virtual environment
```pipenv shell```

Populate the virtual environment tested as follows:
```
pipenv install biopython blinker urllib3 requests pymongo pandas markdown sentry-sdk matplotlib networkx bokeh flask flask-cors psutil scipy
```
After testing, a working package configuration can be saved with
```
pipenv lock
```
and loaded in a production environment with
```
pipenv install --ignore-pipfile
```

The following are required only for simulated data and are not installed by default:
ete3
pyvolve





