pairlist
========

Generates the pair list of atoms that are closer each other than the
given threshold under the orthorhombic periodic boundary conditions.

Usage
-----

See ``pairlist.h`` for the function definition and ``pairlist-test.c``
for usage.

Demo
----

::

    % make
    % ./pairlist-test | wc -l

The input file ``pairlist-test.data`` contains the positions of 432
atoms, each of which is bonded to the four neighbor atoms with bond
distance 2.76. So the total number of bonds should be 864.

Bugs
----
