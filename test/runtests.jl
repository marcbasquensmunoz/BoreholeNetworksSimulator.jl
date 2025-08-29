using BoreholeNetworksSimulator, Test
using Parameters, DataStructures

include("Aqua.jl")

include("utils.jl")

interface_mocks = BoreholeNetworksSimulator.get_all_julia_files_in_dir("mocks/interfaces/")
include.(interface_mocks)
include("mocks/buffer.jl")

# Check interface implementations
include("test_interfaces.jl")

# Run smoke tests
include("api/api.jl")

# Run unit tests
include("mediums/test_GroundMedium.jl")
include("mediums/test_FlowInPorousMedium.jl")

include("constraints/test_TotalHeatLoadConstraint.jl")
include("constraints/test_HeatLoadConstraint.jl")
include("constraints/test_InletTempConstraint.jl")

include("borefields/test_EqualBoreholesBorefield.jl")

include("boreholes/test_SingleUPipeBorehole.jl")

include("methods/test_ConvolutionMethod.jl")
include("methods/test_NonHistoryMethod.jl")
