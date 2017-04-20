#!/usr/bin/env python3


import sys

for filename in sys.argv[1:]:
    with open(filename) as file:
        print("".join(file.readlines()), end="")
        print("")
