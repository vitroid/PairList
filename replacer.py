#!/usr/bin/env python
import sys
import distutils.core
from logging import getLogger
from jinja2 import Template

setup = distutils.core.run_setup("setup.py")


def prefix(L, pre):
    return pre + ("\n" + pre).join(L) + "\n"


d = {
    "sample_py": "".join(open("sample.py").readlines()),
    "sample2_py": "".join(open("sample2.py").readlines()),
    "benchmark_py": "".join(open("benchmark/benchmark.py").readlines()),
    "version": setup.get_version(),
    "package": setup.get_name(),
    "url": setup.get_url(),
    "requires": prefix(setup.install_requires, "* "),
}


# for line in sys.stdin:
#     print(line_replacer(line, d), end="")
t = Template(sys.stdin.read())
print(t.render(**d))
