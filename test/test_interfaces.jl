using RequiredInterfaces
const RI = RequiredInterfaces

@testset "test_Medium_interfaces" begin
    @test RI.check_interface_implemented(Medium, GroundMedium)
    @test RI.check_interface_implemented(Medium, FlowInPorousMedium)
end

@testset "test_Constraint_interfaces" begin
    @test RI.check_interface_implemented(Constraint, TotalHeatLoadConstraint)
    @test RI.check_interface_implemented(Constraint, HeatLoadConstraint)
    @test RI.check_interface_implemented(Constraint, InletTempConstraint)
end

@testset "test_Borehole_interfaces" begin
    @test RI.check_interface_implemented(Borehole, SingleUPipeBorehole)
end

@testset "test_Borefield_interfaces" begin
    @test RI.check_interface_implemented(Borefield, EqualBoreholesBorefield)
    @test RI.check_interface_implemented(Borefield, HeterogeneousBorefield)
end

@testset "test_TimeSuperpositionMethod_interfaces" begin
    @test RI.check_interface_implemented(TimeSuperpositionMethod, ConvolutionMethod)
    @test RI.check_interface_implemented(TimeSuperpositionMethod, NonHistoryMethod)
end

@testset "test_Fluid_interfaces" begin
    @test RI.check_interface_implemented(Fluid, Water)
end

@testset "test_Operator_interfaces" begin
    @test RI.check_interface_implemented(Operator, ConstantOperator)
end
