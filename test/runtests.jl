using WienerHopfCollocation
using GenericLinearAlgebra
using LinearAlgebra
using StaticArrays
using Test

include("../zoo/planewave.jl")
include("../zoo/spanwise.jl")
include("../zoo/sommerfeld.jl")

@testset "WienerHopfCollocation.jl" begin

    freq = 1
    kvec = @SVector [1, 0.5, 0]

    wave = PlaneWave{3,ComplexF64}(kvec, freq)

    z = randn(ComplexF64)

    tol = sqrt(eps(one(real(z))))

    @testset "Sommerfeld" begin

        Ψ = SommerfeldExact(wave)

        @testset "precompute off" begin

            prob = CommonSolve.init(ProblemWHC, SommerfeldWH(wave); ifsave=false)

            pre = precompute(prob)
            sol = preallocate(prob)

            Φ = CommonSolve.solve!(prob, sol, pre)

            @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
            @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol

        end

        @testset "precompute on" begin

            eqn = SommerfeldWH(wave)
            prob = CommonSolve.init(ProblemWHC, eqn; ifsave=true)

            pre = precompute(prob)
            sol = preallocate(prob)

            @testset "save data" begin

                Φ = CommonSolve.solve!(prob, sol, pre)

                @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
                @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol
            end

            @testset "load data" begin
                Φ = CommonSolve.solve!(prob, sol, pre)

                @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
                @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol
            end


        end
    end




    @testset "spanwise" begin

        μ = [1, -0.2, -0.2, 0.1, 0.1, 0.01, 0.02, 0.01, 0.02 + 0.1im]
        eqn = SpanwiseWH(wave, μ, 5)

        prob = CommonSolve.init(ProblemWHC, eqn; ifsave=true)

        Φ = CommonSolve.solve!(prob)

        ε = sqrt(eps(ApproxFun.prectype(prob)))
        @test norm(Φ(ε, true) - Φ(-ε, true), Inf) < cbrt(ε)
        @test norm(Φ(ε, false) - Φ(-ε, false), Inf) < cbrt(ε)
    end


    #Failing (at least) due to needing specialisations in WienerHopf.jl
    @testset "big" begin

        freq = big(1im)
        kvec = @SVector [big(2+1im), big(1 // 2), big(0)]

        wave = PlaneWave{3,Complex{BigFloat}}(kvec, freq)

        z = big(randn(ComplexF64))

        Ψ = SommerfeldExact(wave)

        prob = CommonSolve.init(ProblemWHC, SommerfeldWH(wave);
            ifsave=false,
            domain = SqrtLine{false}(big(float(0 + 0im))),
            m = 400,
            n = 300)

        @test ApproxFun.prectype(prob) == BigFloat

        @test WienerHopfCollocation.isbig(ApproxFun.prectype(prob))

        Φ = CommonSolve.solve!(prob)

        ε = sqrt(eps(ApproxFun.prectype(prob)))
        @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < ε
        @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < ε

    end

    @testset "scaled" begin

        @testset "Sommerfeld" begin

            Ψ = SommerfeldExact(wave)

            @testset "precompute off" begin

                prob = CommonSolve.init(ProblemWHC, SommerfeldWH(wave);
                        ifsave=false,
                        domain = ScaledLine(SqrtLine{-1/4}(0.0),2.0)
                        )

                pre = precompute(prob)
                sol = preallocate(prob)

                Φ = CommonSolve.solve!(prob, sol, pre)

                @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
                @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol

            end

            @testset "precompute on" begin

                eqn = SommerfeldWH(wave)
                prob = CommonSolve.init(ProblemWHC, SommerfeldWH(wave);
                ifsave=true,
                domain = ScaledLine(SqrtLine{-1/4}(0.0),2.0)
                )


                pre = precompute(prob)
                sol = preallocate(prob)

                @testset "save data" begin

                    Φ = CommonSolve.solve!(prob, sol, pre)

                    @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
                    @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol
                end

                @testset "load data" begin
                    Φ = CommonSolve.solve!(prob, sol, pre)

                    @test (1 - Φ(z, true)[1] / Ψ(z, true) |> norm) < tol
                    @test (1 - Φ(z, false)[1] / Ψ(z, false) |> norm) < tol
                end


            end
        end
    end
    # Clean up
    rm("cauchydata", recursive=true)
end
