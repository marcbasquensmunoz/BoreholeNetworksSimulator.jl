module BoreholeNetworksSimulator

using LinearAlgebra
using GeometryTypes
using Parameters
using CoolProp
using StaticArrays
using SparseArrays
using LinearSolve
using JLD2
using CSV
using DataFrames
using DataStructures
using QuadGK
using FastGaussQuadrature
using LegendrePolynomials
using SpecialFunctions
using RequiredInterfaces
using Graphs

using FiniteLineSource

include("utils.jl")

modular = get_all_julia_files_in_dir(joinpath(@__DIR__, "modular"))
sort_dependencies!(modular, ["interfaces", "network", "core"])
include.(modular)
include("CoolProp/load_properties.jl")
include("CoolProp/interpolate.jl")

export Fluid, Water, EthanolMix, GlycolMix
export Borehole, SingleUPipeBorehole
export Borefield, EqualBoreholesBorefield, RectangularBorefield, HeterogeneousBorefield
export Medium, GroundMedium, FlowInPorousMedium
export Constraint, TotalHeatLoadConstraint, HeatLoadConstraint, InletTempConstraint, constant_HeatLoadConstraint, uniform_HeatLoadConstraint, constant_InletTempConstraint, uniform_InletTempConstraint
export TimeSuperpositionMethod, ConvolutionMethod, NonHistoryMethod
export BoundaryCondition, NoBoundary, DirichletBoundaryCondition, NeumannBoundaryCondition
export Approximation, MidPointApproximation, MeanApproximation
export SimulationOptions, SimulationContainers
export BoreholeNetwork, Valve, BoreholeOperation
export connect_to_source!, connect_to_sink!, connect!, all_series_network, all_parallel_network, compute_mass_flows!, source, sink, connect_in_series!, equal_valve, valve, absolute_valve, initialize_mass_flows, boreholes_in_branch, first_bhs_in_branch, connect_in_parallel!
export Operator, SimpleOperator, operate
export ConstantOperator
export simulate!, simulate_steps!, compute_parameters, load_cache!, save_cache, initialize, reset!
export n_branches, all_parallel_network, all_series_network
export ThermophysicalProperties
export extract_Tfin, extract_Tfout, extract_Tb, extract_q

end # module