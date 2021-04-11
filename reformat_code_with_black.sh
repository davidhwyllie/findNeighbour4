# standardise syntax
# suitable for running as part of CI

pipenv run black *.py -l 127 -t py38 --check
pipenv run black findn -l 127 -t py38 --check
pipenv run black test -l 127 -t py38 --check
pipenv run black demo -l 127 -t py38 --check
pipenv run black pca -l 127 -t py38 --check
pipenv run black snpclusters -l 127 -t py38 --check
pipenv run black tree -l 127 -t py38 --check
pipenv run black utils -l 127 -t py38 --check