# {{package}}
Generates the pair list of atoms that are closer to each other than the
given threshold under the periodic boundary conditions.

version {{version}}

## Usage

See `pairlist.h` for the function definition and `pairlist-test.c` for usage.

Python API is served in pairlist.py. Here is a sample code to use it.

```python
{{sample_py}}
```

```python
{{sample2_py}}
```

## Demo

It requires [GenIce](https://github.com/vitroid/GenIce) to make the test data.

```shell
% make test
% make test2
% make test3
```

## Requirements

{{requires}}

## Bugs

