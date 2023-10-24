#!/usr/bin/env python

# from distutils.core import setup, Extension
import codecs
import os
import re

# This may cause an error in installation because numpy is still not installed.
from numpy.distutils.misc_util import get_numpy_include_dirs
from setuptools import Extension, dist, setup

# Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(
    os.path.join(os.path.dirname(__file__), "pairlist.py"), encoding="utf8"
) as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))

# bootstrap numpy
dist.Distribution().fetch_build_eggs(["numpy"])

setup(
    ext_modules=[
        Extension(
            "cpairlist",
            ["cpairlist.c", "pairlist.c"],
            extra_compile_args=[
                "-std=c99",
            ],
            include_dirs=get_numpy_include_dirs(),
        )
    ],
    headers=["pairlist.h"],
    name="PairList",
    version=metadata["version"],
    zip_safe=False,
    py_modules=["pairlist"],
    description="Generate neighbor list for the particles in a periodic boundary cell.",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
    ],
    author="Masakazu Matsumoto",
    author_email="vitroid@gmail.com",
    url="https://github.com/vitroid/PairList/",
    keywords=[
        "pairlist",
    ],
    license="MIT",
    install_requires=[
        "numpy",
    ],
    entry_points={"console_scripts": ["pairlist = pairlist:main"]},
)
