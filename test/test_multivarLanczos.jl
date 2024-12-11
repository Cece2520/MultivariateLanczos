using ITensors, LinearAlgebra

@testset "MultivariateLanczos" begin
    @testset "dummy" begin
        n=5;
        m=2;
        A1 = diagm(randn(n));
        A2 = diagm(randn(n));
        input_A = Dict(Pair("x", A1), Pair("y", A2));
        b = randn(n);

        T = MultivariateLanczos(input_A,b,n,m);
        @show T
    end
    # @testset "real bivariate diagonal" begin
    #     n=5;
    #     A1 = diagm(randn(n));
    #     A2 = diagm(randn(n));
    #     input_A = [Pair("x", A1), Pair("y", A2)];
    #     b = randn(n);

    #     T = MultivariateLanczos(input_A,b, n)
    #     @show T
    # end
end