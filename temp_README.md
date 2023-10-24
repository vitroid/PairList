# {{package}}
Generates the pair list of atoms that are closer to each other than the
given threshold under the periodic boundary conditions.

version {{version}}

## Usage

See `pairlist.h` for the function definition and `pairlist-test.c` for usage.

Python API is served in pairlist.py. Here is a sample code to use it.

<!-- ```python
{{sample_py}}
```

```python
{{sample2_py}}
``` -->

<!-- ## Benchmark tests -->

To find the neighbors in a face-centered cubic lattice of size 10x10x10 on a MacBook Air 2021 (Apple Silicon),

```shell
$ python benchmark.py
INFO crude: Neighboring pair list by a crude double loop.
INFO crude: 18024 ms
INFO crude: end.
24000 pairs
INFO numpyish: Neighboring pair list by numpy fancy array.
INFO numpyish: 741 ms
INFO numpyish: end.
24000.0 pairs
INFO pairlist_py: Neighboring pair list by pairlist in pure python.
INFO pairlist_py: 125 ms
INFO pairlist_py: end.
24000 pairs
INFO pairlist_c: Neighboring pair list by pairlist in c.
INFO pairlist_c: end.
INFO pairlist_c: 16 ms
24000 pairs
```

```python
{{benchmark_py}}
```

![benchmark](https://github.com/vitroid/PairList/raw/master/benchmark/benchmark.png)

## Algorithm

A simple cell division algorithm is implemented.

## Demo

It requires [GenIce](https://github.com/vitroid/GenIce) to make the test data.

```shell
% make test
```

## Requirements

{{requires}}

## Bugs

