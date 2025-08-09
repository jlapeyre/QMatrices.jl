# Spin 1 matrices
# The scale factors are chosen so that the eigenvalues of each matrix are {-1, 0, 1}
# I think we may want to multiply these by a factor of two. This is
# the same scaling we use for X, Y, Z. But this needs to be double checked.
#
# The eigenvalues of the spin-1/2 ops are {-1/2, 1/2}
# We scale these for the Pauli matrices to {-1, 1}.
# The eigenvalues of the spin-1 ops are {-1, 0, 1}. So the same scaling would give
# {-2, 0, 2}.
# However each Pauli op generates a finite group: Eg, X == X^3.
# The spin-1/2 ops do not generate a finite group. The elements approach zero.
# Each spin-1 op, *does* generate a finite group S1x^3 == S1x.
# This is without scaling. I'm a bit puzzled by this.

const S1x =
   (1/sqrt(2)) * [
        0 1 0;
        1 0 1;
        0 1 0
    ]

const S1y =
      (1/sqrt(2)) * [
        0  -im  0;
        im 0  -im;
        0 im 0
    ]

const S1z =
    float(
    [
        1 0 0;
        0 0 0;
        0 0 -1
    ]
)
