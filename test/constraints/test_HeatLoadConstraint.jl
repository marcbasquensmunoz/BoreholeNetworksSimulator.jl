import BoreholeNetworksSimulator: constraints_coeffs!, constraints_b!, BorefieldMock

@testset "test_HeatLoadConstraint_M_parallel" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    H = 23.

    Q_tot = ones(Nbr, Nt)
    constraint = HeatLoadConstraint(Q_tot)

    M = zeros(Nbr, 4*Nb)

    network = BoreholeNetwork([[1], [2], [3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock(H = H * ones(3))

    constraints_coeffs!(M, constraint, operation, borefield)

    expected = [(1, 3Nb+1, H), (2, 3Nb+2, H), (3, 3Nb+3, H)]
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

    network = BoreholeNetwork([[1, 2, 3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock(H = H * ones(3))

    constraints_coeffs!(M, constraint, operation, borefield)

    expected = [(1, 3Nb+1, H), (1, 3Nb+2, H), (1, 3Nb+3, H)]
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

    network = BoreholeNetwork([[1, 2], [3, 4], [5]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock(H = H * ones(5))

    constraints_coeffs!(M, constraint, operation, borefield)

    expected = [(1, 3Nb+1, H), (1, 3Nb+2, H), (2, 3Nb+3, H), (2, 3Nb+4, H), (3, 3Nb+5, H)]
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

    network = BoreholeNetwork([[1], [2], [3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))

    b_time = zeros(Nbr, Nt)
    for step in 1:Nt
        constraints_b!(b, constraint, operation, step)
        b_time[:, step] .= b
    end

    @test b_time == Q_tot
end