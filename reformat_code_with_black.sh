# standardise syntax
# suitable for running as pre-commit  CI

pipenv run black *.py -l 127 -t py38 
pipenv run black findn -l 127 -t py38 
pipenv run black test -l 127 -t py38 
pipenv run black demo -l 127 -t py38 
pipenv run black pca -l 127 -t py38
pipenv run black snpclusters -l 127 -t py38 
pipenv run black tree -l 127 -t py38 
pipenv run black utils -l 127 -t py38 