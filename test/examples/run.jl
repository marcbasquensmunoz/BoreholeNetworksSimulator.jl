
@info "Running Braedstrup"
include("../../examples/Braedstrup/main.jl")
@info "Running g-function"
include("../../examples/g-function/main.jl")
@info "Running Tekniska: constant m"
include("../../examples/tekniska/constant_m.jl")
@info "Running Tekniska: proportional m"
include("../../examples/tekniska/prop_m.jl")
@info "Running Tekniska: toggle"
include("../../examples/tekniska/toggle.jl")
@info "Running partialization"
include("../../examples/partialization/interval.jl")
include("../../examples/partialization/temp_diff.jl")
@info "Running reversal"
include("../../examples/partialization/main.jl")
@info "Running turbulence"
include("../../examples/turbulence/main.jl")
