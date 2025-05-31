@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    nothing
    using QMatrices
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        # These don't take very long to compile anyway.
        # Especially for larger values of N, the run time is much larger than compile time.
        RZpi(1)^2
        random_unitary(2)
        random_unitary(4)
        random_special_unitary(4)
        random_unitary_hermitian(4)
        random_orthogonal(8)
    end
end
