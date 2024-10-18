module BoreholeNetworksSimulator

using LinearAlgebra
using GeometryTypes
using Parameters
using CoolProp
using ExponentialUtilities
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

using FiniteLineSource

include("utils.jl")

modular = get_all_julia_files_in_dir(joinpath(@__DIR__, "modular"))
sort_dependencies!(modular, ["interfaces", "core"])
include.(modular)
include("CoolProp/load_properties.jl")
include("CoolProp/interpolate.jl")

export Fluid, Water, EthanolMix
export Borehole, SingleUPipeBorehole
export Borefield, EqualBoreholesBorefield, RectangularBorefield
export Medium, GroundMedium, FlowInPorousMedium
export Constraint, TotalHeatLoadConstraint, HeatLoadConstraint, InletTempConstraint, constant_HeatLoadConstraint, uniform_HeatLoadConstraint, constant_InletTempConstraint, uniform_InletTempConstraint
export TimeSuperpositionMethod, ConvolutionMethod, NonHistoryMethod
export BoundaryCondition, NoBoundary, DirichletBoundaryCondition, AdiabaticBoundaryCondition
export Approximation, MidPointApproximation, MeanApproximation
export SimulationOptions, SimulationContainers, BoreholeNetwork
export BoreholeOperation, Operator, SimpleOperator, operate
export simulate!, compute_parameters, load_cache!, save_cache, initialize
export n_branches, all_parallel_network, all_series_network
export ThermophysicalProperties

end # module