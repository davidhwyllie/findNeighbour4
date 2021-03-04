Running sphinx is not very straightforward.

Issues:
1. Install sphinx within virtual environment; you have then to run the sphinx commands (e.g. sphinx_build) from within the pipenv shell.  
Otherwise, it won't import modules only installed within the pipenv.
2. Modules names like findNeighbour4-server are considered invalid by sphinx and it won't import them.  findNeighbour4_server is OK.  
To fix this, we'll need to do quite a lot of renaming including of findNeighbour4-clustering.py which is called from some unit tests.
3. From within the pipenv shell, 
sphinx build . _build
does appear to work to a limited extent.  Code cross referencing is not working.
4. A lot of unittests get documented.  This is not terribly helpful.
It is probably possible to avoid this


http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
