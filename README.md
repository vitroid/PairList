# pairlist

Generates the pair list of atoms that are closer to each other than the
given threshold under the periodic boundary conditions.

version 0.6.2

## Usage

See `pairlist.h` for the function definition and `pairlist-test.c` for usage.

Python API is served in pairlist.py. The API document is [here](https://vitroid.github.io/PairList/pairlist.html).

See [benchmark.ipynb](https://colab.research.google.com/github/vitroid/PairList/blob/master/benchmark.ipynb) for the comparison with other methods.
<a href="https://colab.research.google.com/github/vitroid/PairList/blob/master/benchmark.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
![benchmark](https://github.com/vitroid/PairList/raw/master/benchmark/benchmark.png)

## Algorithm

A simple cell division algorithm is implemented.

## Demo

It requires [GenIce](https://github.com/vitroid/GenIce) to make the test data.

```shell
% make test
```

## Requirements

* python
* numpy


## Bugs
