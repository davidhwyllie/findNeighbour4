Testing RDBMS storage systems other than Oracle
===============================================

In theory it should be possible to use RDBMS other than Oracle ADW with findNeighbour4, but this would need to be tested.

If you wish to do this, it is easy to do so.

The script
```
test/testrdbmsstore.py
```

contains unit tests of all interactions the server makes with the database.
It can be invokes as follows
```
pipenv run python3 -m unittest test/testrdbmsstore.py
```

By default, it will test   
* Sqlite
* Any configures Oracle database connections identified in the [relevant configuration file](database_credentials.md) by a key starting with 'unittest_ora'

However, you can other database connections readily.
In the *Test_Database* class, additional connections can be added

```
class Test_Database(unittest.TestCase):
    """establishes database connection strings for cross-database testing.
    Currently tests OCI (if relevant environment variables are set) and Sqlite

    To test other databases, such as MySql, add the relevant connection string &
    database name to the dictionary self.engines"""

    def setUp(self):
        self.engines = {}

        # add additional connection strings to unit test on different databases
        self.engines["Sqlite"] = "sqlite://"  # in memory sqlite
        # self.engines["mysql"] = "mysql+pymysql://root:root@localhost:3306/test_db"  
        # NOTE: accessing a database test_db, with user root and password root; known issues, see above.
```

Thus, if you uncomment the line, which is provided as an example (we have not fully tested mysql with findneighbour4, and cannot guarantee it will work)

```
# self.engines["mysql"] = "mysql+pymysql://root:root@localhost:3306/test_db"
```

The code will run all unit tests against a mysql on localhost, authenticating by username & password both as 'root', on test_db.  You will need to add additional driver modules to the virtual environment, e.g. 

``` 
pipenv install MySQL-python --skip-lock
```
but no other changes shold be needed to start testing.  

Please let us know if you have success with other databases.

