# QMatrices

Plain-data matrices and other objects for quantum computation.

The goal is to make this sufficiently simple, complete, and narrow in scope, that there
is no temptation to write the millionth definition of explicit representations of standard matrices for quantum information rather than depend on this package.

There are no `struct`s or types defined here.
Instead, a minimal number of the most important Julia types are used.
Thus far `Array` and `StaticArrays.SArray` are used.

For performance, the exported fixed gates are of type `SArray`.
The equivalent dynamic (heap allocated, of type `Array`) are in submodule `QMatrices.Dynamic`.
