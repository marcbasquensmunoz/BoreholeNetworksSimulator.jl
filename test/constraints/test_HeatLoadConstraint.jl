import BoreholeNetworksSimulator: constraints_coeffs!, constraints_b!

@testset "test_HeatLoadConstraint_M_parallel" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    H = 23.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = all_parallel_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    compute_mass_flows!(mass_flows, network, operation)

    borefield = BorefieldMock(H = H * ones(3))

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 3Nb+1, H), (2, 3Nb+2, H), (3, 3Nb+3, H)]
    @test test_sparse_matrix(M, expected)
end


@testset "test_HeatLoadConstraint_M_parallel_no_flow" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    H = 23.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = all_parallel_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=zeros(Nbr))
    compute_mass_flows!(mass_flows, network, operation)

    borefield = BorefieldMock(H = H * ones(3))

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [
                (1, 1, 1.), (1, 2, -1.), 
                (2, 3, 1.), (2, 4, -1.), 
                (3, 5, 1.), (3, 6, -1.), 
                ]
    @test test_sparse_matrix(M, expected)
end


@testset "test_HeatLoadConstraint_M_parallel_no_flow_partial" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    H = 23.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = all_parallel_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=vcat(zeros(1), ones(1), zeros(1)))
    compute_mass_flows!(mass_flows, network, operation)

    borefield = BorefieldMock(H = H * ones(3))

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [
                (1, 1, 1.), (1, 2, -1.), 
                (2, 3Nb+2, H),
                (3, 5, 1.), (3, 6, -1.), 
                ]
    @test test_sparse_matrix(M, expected)
end

@testset "test_HeatLoadConstraint_M_series" begin
    Nbr = 1
    Nb = 3
    Nt = 1

    H = 32.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    borefield = BorefieldMock(H = H * ones(3))

    network = all_series_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    compute_mass_flows!(mass_flows, network, operation)

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 3Nb+1, H), (1, 3Nb+2, H), (1, 3Nb+3, H)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_HeatLoadConstraint_M_series_no_flow" begin
    Nbr = 1
    Nb = 3
    Nt = 1

    H = 32.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    borefield = BorefieldMock(H = H * ones(3))

    network = all_series_network(Nb)
    mass_flows = initialize_mass_flows(network)
    operation = BoreholeOperation(network=network, mass_flows=zeros(Nbr))
    compute_mass_flows!(mass_flows, network, operation)

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 1, 1.), (1, 2, -1.)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_HeatLoadConstraint_M_mixed" begin
    Nbr = 3
    Nb = 5
    Nt = 1

    H = 51.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = BoreholeNetwork(5)
    connect_to_source!(network, [1, 3, 5])
    connect!(network, 1, 2)
    connect!(network, 3, 4)
    connect_to_sink!(network, [2, 4, 5])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    borefield = BorefieldMock(H = H * ones(5))

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [(1, 3Nb+1, H), (1, 3Nb+2, H), (2, 3Nb+3, H), (2, 3Nb+4, H), (3, 3Nb+5, H)]
    @test test_sparse_matrix(M, expected)
end

@testset "test_HeatLoadConstraint_M_mixed_no_flow" begin
    Nbr = 3
    Nb = 5
    Nt = 1

    H = 51.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = BoreholeNetwork(5)
    connect_to_source!(network, [1, 3, 5])
    connect!(network, 1, 2)
    connect!(network, 3, 4)
    connect_to_sink!(network, [2, 4, 5])
    operation = BoreholeOperation(network=network, mass_flows=zeros(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    borefield = BorefieldMock(H = H * ones(5))

    constraints_coeffs!(M, constraint, borefield, network, mass_flows)

    expected = [
                (1, 1, 1.), (1, 2, -1.), 
                (2, 5, 1.), (2, 6, -1.), 
                (3, 9, 1.), (3, 10, -1.)
                ]
    @test test_sparse_matrix(M, expected)
end

@testset "test_HeatLoadConstraint_b" begin
    Nbr = 3
    Nb = 3
    Nt = 10

    Q_tot = ones(Nbr, Nt)
    for i in axes(Q_tot, 2)
        Q_tot[:, i] .= Float64(i)
    end
    constraint = HeatLoadConstraint(Q_tot)

    b = zeros(Nbr)

    network = all_parallel_network(Nb)
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    b_time = zeros(Nbr, Nt)
    for step in 1:Nt
        constraints_b!(b, constraint, network, mass_flows, step)
        b_time[:, step] .= b
    end

    @test b_time == Q_tot
end

@testset "test_HeatLoadConstraint_b_no_flow" begin
    Nbr = 3
    Nb = 3
    Nt = 10

    Q_tot = ones(Nbr, Nt)
    for i in axes(Q_tot, 2)
        Q_tot[:, i] .= Float64(i)
    end
    constraint = HeatLoadConstraint(Q_tot)

    b = zeros(Nbr)

    network = all_parallel_network(Nb)
    operation = BoreholeOperation(network=network, mass_flows=zeros(Nbr))
    mass_flows = initialize_mass_flows(network)
    compute_mass_flows!(mass_flows, network, operation)

    b_time = zeros(Nbr, Nt)
    for step in 1:Nt
        constraints_b!(b, constraint, network, mass_flows, step)
        b_time[:, step] .= b
    end

    @test b_time == zeros(Nbr, Nt)
end