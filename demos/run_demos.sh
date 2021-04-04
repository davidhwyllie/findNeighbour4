#!/bin/bash

# script to run benchmarks

echo "Running benchmark 1"
pipenv run python3 demo_benchmark.py ../demos/benchmark/config/config1.json

echo "Running benchmark 2"
pipenv run python3 demo_benchmark.py ../demos/benchmark/config/config2.json

echo "Running benchmark 3"
pipenv run python3 demo_benchmark.py ../demos/benchmark/config/config3.json

echo "Running benchmark 4"
pipenv run python3 demo_benchmark.py ../demos/benchmark/config/config4.json

echo "Startup complete"


