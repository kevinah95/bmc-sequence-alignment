# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

try:
    long_description = open("README.md").read()
except IOError:
    long_description = ""

setup(
    name="bmc-sequence-alignment",
    version="0.1.0",
    description="Sequence alignments for Computational Molecular Biology",
    license="MIT",
    author="Kevin Hernandez Rostran",
    packages=find_packages(),
    install_requires=[],
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
    ]
)
