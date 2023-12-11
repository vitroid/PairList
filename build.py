from setuptools.command.build_ext import build_ext
from setuptools import Extension
from distutils.errors import (
    DistutilsPlatformError,
    CCompilerError,
    DistutilsExecError,
    DistutilsPlatformError,
)

import numpy as np

ext_modules = [
    Extension(
        "cpairlist",
        include_dirs=[
            np.get_include(),
        ],
        sources=["csource/cpairlist.c", "csource/pairlist.c"],
        extra_compile_args=[
            "-std=c99",
        ],
    ),
]


class BuildFailed(Exception):
    pass


class ExtBuilder(build_ext):
    def run(self):
        try:
            build_ext.run(self)
        except (DistutilsPlatformError, FileNotFoundError):
            raise BuildFailed("File not found. Could not compile C extension.")

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError, ValueError):
            raise BuildFailed("Could not compile C extension.")


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    setup_kwargs.update(
        {"ext_modules": ext_modules, "cmdclass": {"build_ext": ExtBuilder}}
    )
