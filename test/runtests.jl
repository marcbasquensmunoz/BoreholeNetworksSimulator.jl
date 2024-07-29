using BoreholeNetworksSimulator, Test

include("utils.jl")

dir = @__DIR__
include("test_interfaces.jl")

include("mediums/test_GroundMedium.jl")
include("mediums/test_FlowInPorousMedium.jl")

include("constraints/test_HeatLoadConstraint.jl")
include("constraints/test_InletTempConstraint.jl")

include("borefields/test_EqualBoreholesBorefield.jl")

include("boreholes/test_SingleUPipeBorehole.jl")

include("methods/test_ConvolutionMethod.jl")
include("methods/test_NonHistoryMethod.jl")