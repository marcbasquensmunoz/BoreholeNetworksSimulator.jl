import BoreholeNetworksSimulator: n_boreholes, get_H, get_rb, segment_coordinates, internal_model_coeffs!, internal_model_b!

@testset "test_HeterogeneousBorefield" begin
    H = 100.
    D = 20.
    rb = 0.1

    x_pos = 0.
    y_pos = 0.
    σ = 5.
    positions = [(x_pos, y_pos), (x_pos, y_pos+σ), (x_pos, y_pos+2σ), (x_pos+σ, y_pos)]
    boreholes = [BoreholeMock(H=H, D=D, rb=rb) for _ in 1:length(positions)]
    borefield = HeterogeneousBorefield(boreholes=boreholes, positions=positions)

    @test n_boreholes(borefield) == length(positions)
    for i in eachindex(positions)
        @test get_H(borefield, i) == H
        @test get_rb(borefield, i) == rb
    end
    @test segment_coordinates(borefield, 1) == (x_pos, y_pos, D, H)
    @test segment_coordinates(borefield, 2) == (x_pos, y_pos+σ, D, H)
    @test segment_coordinates(borefield, 3) == (x_pos, y_pos+2σ, D, H)
    @test segment_coordinates(borefield, 4) == (x_pos+σ, y_pos, D, H)
end

@testset "test_HeterogeneousBorefield_M" begin
    Ci = 5.
    Co = 2.
    Cb = 3.
    σ = 5.
    positions = [(0., 0.), (0., σ), (0., 2σ)]
    Nb = length(positions)
    boreholes = [BoreholeMock(Ci=Ci, Co=Co, Cb=Cb) for _ in 1:Nb]
    borefield = HeterogeneousBorefield(boreholes=boreholes, positions=positions)

    medium = MediumMock()
    fluid = FluidMock(cpf = 1000.)
    T_fluid = 10 * ones(2*Nb)

    network = all_parallel_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nb))
    compute_mass_flows!(mass_flows, network, operation)

    M = zeros(Nb, 4*Nb)
    internal_model_coeffs!(M, borefield, medium, mass_flows, T_fluid, fluid)

    expected = [
        (1, 1, Ci), (1, 2, Co), (1, 2Nb+1, Cb),
        (2, 3, Ci), (2, 4, Co), (2, 2Nb+2, Cb),
        (3, 5, Ci), (3, 6, Co), (3, 2Nb+3, Cb)
        ]
    @test test_sparse_matrix(M, expected)
end


@testset "test_HeterogeneousBorefield_b" begin
    H = 100.
    D = 20.
    rb = 0.1
    Ci = 5.
    Co = 2.
    Cb = 3.
    σ = 5.
    positions = [(0., 0.), (0., σ), (0., 2σ)]
    Nb = length(positions)
    boreholes = [BoreholeMock(H=H, D=D, rb=rb, Ci=Ci, Co=Co, Cb=Cb) for _ in 1:Nb]
    borefield = HeterogeneousBorefield(boreholes=boreholes, positions=positions)

    b = zeros(Nb)
    internal_model_b!(b, borefield)
    
    @test b == zeros(Nb)
end
