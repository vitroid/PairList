#!/usr/bin/env python
import sys
import os
from genice.tool import line_replacer
import distutils.core

setup = distutils.core.run_setup("setup.py")

d = {
    "%%sample.py%%"   : "".join(open("sample.py").readlines()),
    "%%sample2.py%%"  : "".join(open("sample2.py").readlines()),
    "%%version%%" : setup.get_version(),
    "%%package%%" : setup.get_name(),
    "%%url%%"     : setup.get_url(),
    "%%requires%%": "\n".join(setup.install_requires),
}


for line in sys.stdin:
    print(line_replacer(line, d), end="")
