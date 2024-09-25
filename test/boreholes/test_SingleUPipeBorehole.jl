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
