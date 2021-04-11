from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="findNeighbour4",
    description="A server application for investigating bacterial relatedness using reference-mapped data",
    long_description=Path("./README.md").read_text(),
    url="https://github.com/davidhwyllie/findNeighbour4",
    packages=find_packages(),
)
