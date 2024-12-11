using LinearAlgebra

@testset "MultivariateLanczos" begin
    @testset "dummy" begin
        n=5;
        m=2;
        A1 = diagm(randn(n));
        A2 = diagm(randn(n));
        input_A = [A1,A2];
        b = randn(n);

        (T,Q) = MultivariateLanczos(input_A,b,m);
    end
    @testset "Q orthogonal for real bivariate diagonal" begin
        n=5; m=4;
        A1 = diagm(randn(n));
        A2 = diagm(randn(n));
        input_A = [A1,A2];
        b = randn(n);

        (T,Q) = MultivariateLanczos(input_A,b,m)
        
        #test Q orthonormal
        k = length(Q)
        mat_Q = zeros(n,k)
        index = 1
        for vec in values(Q)
            mat_Q[:,index] = vec
            index += 1
        end
        @test isapprox([norm(mat_Q[:,i]) for i in 1:k], ones(k))
        @test isapprox(mat_Q' * mat_Q,I(k))
    end
    @testset "Q orthogonal for real 4-variable diagonal" begin
        n=20; m=8;
        A1 = diagm(randn(n));
        A2 = diagm(randn(n));
        A3 = diagm(randn(n));
        A4 = diagm(randn(n));
        input_A = [A1,A2,A3,A4];
        b = randn(n);

        (T,Q) = MultivariateLanczos(input_A,b,m)
        
        #test Q orthonormal
        k = length(Q)
        mat_Q = zeros(n,k)
        index = 1
        for vec in values(Q)
            mat_Q[:,index] = vec
            index += 1
        end
        @test isapprox(mat_Q' * mat_Q,I(k))
    end
end