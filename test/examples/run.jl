project_dir = dirname(pwd())
examples_dir = "$project_dir/BoreholeNetworksSimulator.jl/examples"

@info "Running Braedstrup"
include("$examples_dir/Braedstrup/main.jl")
@info "Running g-function"
include("$examples_dir/g-function/main.jl")
@info "Running Tekniska: constant m"
include("$examples_dir/tekniska/constant_m.jl")
@info "Running Tekniska: proportional m"
#include("$examples_dir/tekniska/prop_m.jl")
@info "Running Tekniska: toggle"
#include("$examples_dir/tekniska/toggle.jl")
@info "Running partialization"
include("$examples_dir/partialization/interval.jl")
include("$examples_dir/partialization/temp_diff.jl")
@info "Running reversal"
include("$examples_dir/reversal/main.jl")
@info "Running turbulence"
include("$examples_dir/turbulence/main.jl")
