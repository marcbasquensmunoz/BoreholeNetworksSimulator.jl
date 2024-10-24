
@testset "test_Approximations" begin
    @test MeanApproximation() isa Any
    @test MidPointApproximation() isa Any
end

@testset "test_Boreholes" begin
    @test SingleUPipeBorehole(H = 100., D = 0.) isa Any
end

@testset "test_Borefields" begin
    borehole = SingleUPipeBorehole(H = 100., D = 0.)
    @test EqualBoreholesBorefield(borehole_prototype=borehole, positions=[(0., 0.)]) isa Any
end

@testset "test_Boundary_Conditions" begin
    @test AdiabaticBoundaryCondition() isa Any
    @test DirichletBoundaryCondition() isa Any
    @test NoBoundary() isa Any
end

@testset "test_Constraints" begin
    Nt = 10
    Nbr = 10

    @test HeatLoadConstraint(rand(Nbr, Nt)) isa Any
    @test constant_HeatLoadConstraint(rand(Nbr), Nt) isa Any
    @test uniform_HeatLoadConstraint(rand(Nt), Nbr) isa Any

    @test InletTempConstraint(rand(Nbr, Nt)) isa Any
    @test constant_InletTempConstraint(rand(Nbr), Nt) isa Any
    @test uniform_InletTempConstraint(rand(Nt), Nbr) isa Any

    @test TotalHeatLoadConstraint(rand(Nt)) isa Any    
end

@testset "test_Fluids" begin
    @test Water() isa Any
    @test EthanolMix() isa Any
end

@testset "test_Mediums" begin
    @test FlowInPorousMedium() isa Any
    @test GroundMedium() isa Any
end

@testset "test_Methods" begin
    @test ConvolutionMethod() isa Any
    @test NonHistoryMethod() isa Any
end