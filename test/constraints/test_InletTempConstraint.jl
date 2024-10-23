import BoreholeNetworksSimulator: constraints_coeffs!, constraints_b!

@testset "test_InletTempConstraint_M_parallel" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    T_in = ones(Nbr, Nt)
    constraint = InletTempConstraint(T_in)

    M = zeros(Nbr, 4*Nb)

    borefield = BorefieldMock()

    network = all_parallel_network(Nb)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 1, 1.), (2, 3, 1.), (3, 5, 1.)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_InletTempConstraint_M_series" begin
    Nbr = 1
    Nb = 3
    Nt = 1

    T_in = ones(Nbr, Nt)
    constraint = InletTempConstraint(T_in)

    M = zeros(Nbr, 4*Nb)

    borefield = BorefieldMock()

    network = all_series_network(Nb)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 1, 1.)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_InletTempConstraint_M_mixed" begin
    Nbr = 3
    Nb = 5
    Nt = 1

    T_in = ones(Nbr, Nt)
    constraint = InletTempConstraint(T_in)

    M = zeros(Nbr, 4*Nb)

    network = BoreholeNetwork(5)
    connect_to_source!(network, [1, 3, 5])
    connect!(network, 1, 2)
    connect!(network, 3, 4)
    connect_to_sink!(network, [2, 4, 5])
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    compute_mass_flows!(mass_flows, network, operation)
    
    borefield = BorefieldMock()

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 1, 1.), (2, 2*2+1, 1.), (3, 4*2+1, 1.)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_InletTempConstraint_b" begin
    Nbr = 3
    Nb = 3
    Nt = 10

    T_in = ones(Nbr, Nt)
    for i in axes(T_in, 2)
        T_in[:, i] .= Float64(i)
    end
    constraint = InletTempConstraint(T_in)

    b = zeros(Nbr)

    network = all_parallel_network(Nb)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))

    b_time = zeros(Nbr, Nt)
    for step in 1:Nt
        constraints_b!(b, constraint, operation, step)
        b_time[:, step] .= b
    end

    @test b_time == T_in
end