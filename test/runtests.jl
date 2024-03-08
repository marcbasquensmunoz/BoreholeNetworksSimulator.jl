using BoreholeNetworksSimulator
using Test
using DelimitedFiles
using GeometryTypes

const ϵ = 5*10^-14

function initialize(data, tstep, tmax)
    df = readdlm("$(@__DIR__)/$data", ';', Float64, header=true)[1]
    borehole_positions = [Point2(df[i, 1], df[i, 2]) for i in 1:size(df)[1]]
    borefield = EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=50., D=4.), positions=borehole_positions, medium=GroundWaterMedium(), T0 = 10.)
    parameters = compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
    method = ConvolutionMethod(parameters=parameters, borefield=borefield)
    containers = SimulationContainers(parameters)

    return parameters, method, containers, borefield
end

function load_expected(file)
    readdlm("$(@__DIR__)/$file", ';', Float64)
end

@testset "test" begin
    parameters, method, containers, borefield = initialize("test_data.txt", 8760*3600/12., 8760*3600*10.)
    constraint = InletTempConstraint([i%12 in 1:6 ? 90. : 55. for i = 1:120])
    networks = 
    [
        [
            [1],
            [2]
        ]
    ]

    function operator(i, Tin, Tout, Tb, Δq, Q)
        BoreholeOperation(networks[1], 0.5 .* ones(8), 4182.)
    end

    simulate(parameters=parameters, containers=containers, operator=operator, borefield=borefield, constraint=constraint, method=method)

    expected_result = load_expected("test.csv")

    @test containers.X ≈ expected_result atol = ϵ
end
