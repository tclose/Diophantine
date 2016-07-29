#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="Diophantine",
    version="0.1",
    py_modules=['diophantine'],
    author="Thomas G. Close",
    author_email="tom.g.close@gmail.com",
    description=(
        "A python package for finding small solutions of systems of "
        "diophantine (integer algebra) equations"),
    long_description=open("README.rst").read(),
    license="MIT License",
    keywords=["mathematics", "diophantine", "algebra", "integer",
              "systems"],
    url="http://github.com/tclose/Diophantine/tarball/0.1",
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2',
                 'Topic :: Scientific/Engineering'],
    install_requires=['sympy'],
    tests_require=['nose']
)
