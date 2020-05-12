using QMatrices
using LinearAlgebra: LinearAlgebra

@testset "Basic algebraic props" begin
    for m in (X, Y, Z, H)
        @test isapprox(m^2, LinearAlgebra.I) # TODO: (slow) norm is unnecessary
    end
end

@testset "control" begin
    toffoli = [1 0 0 0 0 0 0 0;
               0 1 0 0 0 0 0 0;
               0 0 1 0 0 0 0 0;
               0 0 0 1 0 0 0 0;
               0 0 0 0 1 0 0 0;
               0 0 0 0 0 1 0 0;
               0 0 0 0 0 0 0 1;
               0 0 0 0 0 0 1 0]

    @test control(X, 2) == toffoli
end
