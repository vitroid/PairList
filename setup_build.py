from setuptools.command.build_ext import build_ext
from setuptools import Extension
from distutils.errors import (
    CCompilerError,
    DistutilsExecError,
    DistutilsPlatformError,
)
import sys
import platform

import numpy as np

# Windows (MSVC) では -std=c99 フラグは不要（サポートされていない）
# Unix系（GCC/Clang）では -std=c99 を追加
compile_args = []
if platform.system() != "Windows":
    compile_args.append("-std=c99")

ext_modules = [
    Extension(
        "cpairlist",
        include_dirs=[
            np.get_include(),
        ],
        sources=["csource/cpairlist.c", "csource/pairlist.c"],
        extra_compile_args=compile_args,
    ),
]


class BuildFailed(Exception):
    pass


class ExtBuilder(build_ext):
    def run(self):
        try:
            build_ext.run(self)
        except (DistutilsPlatformError, FileNotFoundError):
            self.warn_and_continue()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (
            CCompilerError,
            DistutilsExecError,
            DistutilsPlatformError,
            ValueError,
        ):
            self.warn_and_continue()

    def warn_and_continue(self):
        print("\n" + "=" * 80)
        print("WARNING: C extension could not be compiled.")
        print(
            "This is often because Python development headers are missing (e.g. 'Python.h')."
        )
        print("\nTo fix this, you might need to install the python-dev package:")
        print("- On Ubuntu/Debian: sudo apt-get install python3-dev")
        print("- On RHEL/CentOS: sudo yum install python3-devel")
        print("- On macOS: xcode-select --install")
        print("\nInstallation will proceed using a pure Python fallback.")
        print("=" * 80 + "\n")


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    setup_kwargs.update(
        {"ext_modules": ext_modules, "cmdclass": {"build_ext": ExtBuilder}}
    )
