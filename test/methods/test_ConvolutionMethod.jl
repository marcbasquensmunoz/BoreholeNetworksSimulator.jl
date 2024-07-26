import BoreholeNetworksSimulator: method_coeffs!, method_b!, precompute_auxiliaries!, update_auxiliaries!


@testset "test_ConvolutionMethod" begin

end


#=
    method_coeffs!(M, ::TimeSuperpositionMethod, borefield, medium, boundary_condition)
    method_b!(b, ::TimeSuperpositionMethod, borefield, medium, step)
    precompute_auxiliaries!(method::TimeSuperpositionMethod, options)
    update_auxiliaries!(method::TimeSuperpositionMethod, X, borefield, step)
=#