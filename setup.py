# -*- coding: utf-8 -*-
from setuptools import setup

modules = \
['pairlist']
install_requires = \
['numpy>=2.0,<3.0']

setup_kwargs = {
    'name': 'pairlist',
    'version': '0.6',
    'description': 'Generate neighbor list for the particles in a periodic boundary cell.',
    'long_description': '# pairlist\n\nGenerates the pair list of atoms that are closer to each other than the\ngiven threshold under the periodic boundary conditions.\n\nversion 0.6\n\n## Usage\n\nSee `pairlist.h` for the function definition and `pairlist-test.c` for usage.\n\nPython API is served in pairlist.py. The API document is [here](https://vitroid.github.io/PairList/pairlist.html).\n\nSee [benchmark.ipynb](https://colab.research.google.com/github/vitroid/PairList/blob/master/benchmark.ipynb) for the comparison with other methods.\n<a href="https://colab.research.google.com/github/vitroid/PairList/blob/master/benchmark.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>\n![benchmark](https://github.com/vitroid/PairList/raw/master/benchmark/benchmark.png)\n\n## Algorithm\n\nA simple cell division algorithm is implemented.\n\n## Demo\n\nIt requires [GenIce](https://github.com/vitroid/GenIce) to make the test data.\n\n```shell\n% make test\n```\n\n## Requirements\n\n* python\n* numpy\n\n\n## Bugs\n',
    'author': 'vitroid',
    'author_email': 'vitroid@gmail.com',
    'maintainer': 'None',
    'maintainer_email': 'None',
    'url': 'None',
    'py_modules': modules,
    'install_requires': install_requires,
    'python_requires': '>=3.9,<4.0',
}
from setup_build import *
build(setup_kwargs)

setup(**setup_kwargs)
