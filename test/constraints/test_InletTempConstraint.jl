import BoreholeNetworksSimulator: constraints_coeffs!, constraints_b!

@testset "test_InletTempConstraint_M_parallel" begin
    Nbr = 3
    Nb = 3
    Nt = 1

    T_in = ones(Nbr, Nt)
    constraint = InletTempConstraint(T_in)

    M = zeros(Nbr, 4*Nb)

    network = BoreholeNetwork([[1], [2], [3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock()

    constraints_coeffs!(M, constraint, operation, borefield)

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

    network = BoreholeNetwork([[1, 2, 3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock()

    constraints_coeffs!(M, constraint, operation, borefield)

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

    network = BoreholeNetwork([[1, 2], [3, 4], [5]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))
    borefield = BorefieldMock()

    constraints_coeffs!(M, constraint, operation, borefield)

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

    network = BoreholeNetwork([[1], [2], [3]])
    operation = BoreholeOperation(network=network, mass_flows=ones(Nbr))

    b_time = zeros(Nbr, Nt)
    for step in 1:Nt
        constraints_b!(b, constraint, operation, step)
        b_time[:, step] .= b
    end

    @test b_time == T_in
end