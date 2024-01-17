module BTESGroundWaterSimulator

using LinearAlgebra, GeometryTypes
using Distances
using Parameters
using BoreholeResponseFunctions

export GroundWaterFlow, evaluate_relevant_distances, map_unique_pairs
include("mls_simmetries.jl")

export build_matrix!, build_giventerm!, update_b!, solve_full_convolution_step!
include("model_builder.jl")

export rotation, rotation_z
include("geometrical_transformation.jl")

export BoreholePara, resistance_network, coefficient_matrix, deltacircuit, effective_borehole_resistance, uniformTb_koeff
include("innerborehole_model.jl")

export BorefieldProperties, BoreholeProperties, sim1
include("../examples/example1/sim1_encapsulated.jl")

end # module
