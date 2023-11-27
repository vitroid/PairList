#!/usr/bin/env python
import sys
from jinja2 import Template
import toml

project = toml.load("pyproject.toml")


def prefix(L, pre):
    return pre + ("\n" + pre).join(L) + "\n"


project["benchmark_py"] = "".join(open("benchmark/benchmark.py").readlines())

# for line in sys.stdin:
#     print(line_replacer(line, d), end="")
t = Template(sys.stdin.read())
print(t.render(**project))
