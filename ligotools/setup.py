from setuptools import setup, find_packages

VERSION='0.0.1'
setup(
    name="ligotools",
    version=VERSION
    author="Ligo Scientific Collaboration (LSC) and Leif-Martin Sæther Sunde",
    author_email="lmsunde@berkeley.edu",
    install_requires = []
    packages=find_packages()
    classifiers=['Programming language :: Python :: 2']
)