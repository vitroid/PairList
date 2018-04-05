#!/usr/bin/env python

# from distutils.core import setup, Extension
from setuptools import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(os.path.dirname(__file__), 'pairlist.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))

setup(ext_modules=[Extension("cpairlist", ["cpairlist.c", "pairlist.c"])],
      include_dirs=get_numpy_include_dirs(),
      name='PairList',
      version=metadata['version'],
      zip_safe=False,
      py_modules=['pairlist'],
      description='Generate neighbor list for the particles in a periodic boundary cell.',
      #long_description=README + '\n\n' +  CHANGES,
      classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/PairList/',
      keywords=['pairlist',],
      license='MIT',
#      packages=['pairlist',],
#                'genice.molecules',
#                'genice.lattices',
#                'genice.formats',
#                ],
      install_requires=['numpy',],
      entry_points = {
              'console_scripts': [
                  'pairlist = pairlist:main'
              ]
          }
)
