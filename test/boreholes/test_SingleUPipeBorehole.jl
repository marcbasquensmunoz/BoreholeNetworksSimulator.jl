import BoreholeNetworksSimulator: get_H, get_D, get_rb, uniform_Tb_coeffs

@testset "test_SingleUPipeBorehole" begin
    H = 100.
    D = 20.
    rb = 0.1
    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb)

    @test get_H(borehole) == H
    @test get_D(borehole) == D
    @test get_rb(borehole) == rb
end

@testset "test_SingleUPipeBorehole_Tb_coeffs_1" begin
    H = 50.
    D = 20.
    rb = 0.075
    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb)

    λ = 5.
    mass_flow = 2.
    Tref = 10.
    fluid = Fluid(cpf=4000.)

    Ei, Eo, Eb = uniform_Tb_coeffs(borehole, λ, mass_flow, Tref, fluid)

    @test Ei == -0.9606199814627016
    @test Eo == 1.0422072067496164
    @test Eb == -0.08158722528691487
end

@testset "test_SingleUPipeBorehole_Tb_coeffs_2" begin
    H = 100.
    D = 20.
    rb = 0.1
    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb)

    λ = 3.
    mass_flow = 10.
    Tref = 10.
    fluid = Fluid(cpf=4000.)

    Ei, Eo, Eb = uniform_Tb_coeffs(borehole, λ, mass_flow, Tref, fluid)

    @test Ei == -0.9881396306590124
    @test Eo == 1.0122556128584963
    @test Eb == -0.024115982199483854
end


