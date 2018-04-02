from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

setup(ext_modules=[Extension("cpairlist", ["cpairlist.c", "pairlist.c"])],
      include_dirs=get_numpy_include_dirs())
