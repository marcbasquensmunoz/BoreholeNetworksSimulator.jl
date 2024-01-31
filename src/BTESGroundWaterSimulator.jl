module BTESGroundWaterSimulator

using LinearAlgebra, GeometryTypes
using Distances
using Parameters
using BoreholeResponseFunctions
using CoolProp

include("utils.jl")

export GroundWaterFlow, evaluate_relevant_distances, map_unique_pairs
include("mls_simmetries.jl")

export build_matrix!, build_giventerm!, update_b!, solve_full_convolution_step!
include("model_builder.jl")

export rotation, rotation_z
include("geometrical_transformation.jl")

export BoreholePara, resistance_network, coefficient_matrix, deltacircuit, effective_borehole_resistance, uniformTb_koeff
export heat_transfer_coefficient
include("innerborehole_model.jl")

export Borehole, SingleUPipeBorehole
export Borefield, EqualBoreholesBorefield
export Medium, GroundMedium, GroundWaterMedium
export Constraint, HeatLoadConstraint, InletTempConstraint
export Method, ConvolutionMethod
export SimulationParameters, SimulationContainers, BoreholeOperation, compute_parameters, load_cache!, save_cache
export simulate

modular = get_all_julia_files_in_dir("$(@__DIR__)/modular")
sort_dependencies!(modular, ["interfaces/"])
include.(modular)

end # module
