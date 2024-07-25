using BoreholeNetworksSimulator, Test

dir = @__DIR__
include("test_interfaces.jl")
include("mediums/test_GroundMedium.jl")
include("mediums/test_FlowInPorousMedium.jl")