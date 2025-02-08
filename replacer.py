#!/usr/bin/env python
import sys
from jinja2 import Template
import toml


def prefix(L, pre):
    return pre + ("\n" + pre).join(L) + "\n"


d = toml.load("pyproject.toml")
d |= {
    # "benchmark_py": "".join(open("benchmark/benchmark.py").readlines()),
    # "version": setup.get_version(),
    # "package": setup.get_name(),
    # "url": setup.get_url(),
    # "requires": prefix(setup.install_requires, "* "),
}


# for line in sys.stdin:
#     print(line_replacer(line, d), end="")
t = Template(sys.stdin.read())
print(t.render(**d))
