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

function Tb_coeffs(;H, D, rb, λ, mf, T, cpf)
    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb)
    fluid = Fluid(cpf=cpf, name="INCOMP::MEA-20%")
    uniform_Tb_coeffs(borehole, λ, mf, T, fluid)
end

@testset "test_SingleUPipeBorehole_Tb_coeffs" begin
    @test Tb_coeffs(H=100., D=0., rb=0.1, λ=3., mf=5., T=10., cpf=4000.) == (-0.9767097344219599, 1.0248603726039356, -0.048150638181975736)
    @test Tb_coeffs(H=50., D=0., rb=0.1, λ=3., mf=5., T=10., cpf=4000.) == (-0.9881609464332975, 1.012231541811421, -0.024070595378123394)
    @test Tb_coeffs(H=100., D=0., rb=0.075, λ=3., mf=5., T=10., cpf=4000.) == (-0.9681140081263698, 1.0336253102416069, -0.06551130211523704)
    @test Tb_coeffs(H=100., D=0., rb=0.1, λ=10., mf=5., T=10., cpf=4000.) == (-0.9767445614358461, 1.0249574330282976, -0.04821287159245147)
    @test Tb_coeffs(H=100., D=0., rb=0.1, λ=3., mf=1., T=10., cpf=4000.) == (-0.8988859213845516, 1.1395608553237033, -0.24067493393915174)
    @test Tb_coeffs(H=100., D=0., rb=0.1, λ=3., mf=5., T=10., cpf=1000.) == (-0.9159058832119038, 1.1092651648954326, -0.19335928168352878)
end
