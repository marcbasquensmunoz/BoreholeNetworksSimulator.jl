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
using PythonCall
using DataStructures
using QuadGK
using FastGaussQuadrature
using LegendrePolynomials
using SpecialFunctions
using RequiredInterfaces

using BoreholeResponseFunctions
using FiniteLineSource

include("utils.jl")

modular = get_all_julia_files_in_dir("$(@__DIR__)/modular")
sort_dependencies!(modular, ["interfaces/", "core/"])
include.(modular)

export Borehole, SingleUPipeBorehole
export Borefield, EqualBoreholesBorefield
export Medium, GroundMedium, FlowInPorousMedium
export Constraint, HeatLoadConstraint, InletTempConstraint, constant_HeatLoadConstraint, uniform_HeatLoadConstraint, constant_InletTempConstraint, uniform_InletTempConstraint
export TimeSuperpositionMethod, ConvolutionMethod, NonHistoryMethod
export BoundaryCondition, NoBoundary, DirichletBoundaryCondition
export SimulationOptions, SimulationContainers, BoreholeOperation, BoreholeNetwork, Fluid
export simulate!, compute_parameters, load_cache!, save_cache, initialize

end # module
