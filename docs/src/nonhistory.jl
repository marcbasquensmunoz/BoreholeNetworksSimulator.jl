# # Non-history method
#
# The non-history method is a time superposition method introduced in [1] whose computational complexity
# in the number of time steps ``N_t`` is ``\mathcal{O}\left( N_t \right)``. Recall that the standard
# way to do time superposition is via the convolution of the load with the response, which implemented
# via the Fast Fourier Transform, yields a computational complexity of ``\mathcal{O}\left( N_t \log{N_t} \right)``.
# This means that using the non-history method in simulations allows for finer time steps.

# To show this, let us run a simulation with hourly time steps, with a duration of 20 years (so `175200` time steps),
# with both the convolution and the non-history time superposition methods.
# Let us define an example, very similar to 
using BoreholeNetworksSimulator

Δt = 8760*3600/12.     # Hourly time step
Nt = 2

medium = GroundMedium(α=1e-6, λ=3., T0=10.)
positions = [(0., 0.), (0., 5.)]
borehole = SingleUPipeBorehole(H=100., D=10.)
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
configurations = [BoreholeNetwork([[1], [2]])]
constraint = constant_HeatLoadConstraint(5 .* ones(BoreholeNetworksSimulator.n_boreholes(borefield)), Nt)

function operator(i, Tin, Tout, Tb, q, configurations)
    BoreholeOperation(configurations[1], 2 .* ones(2))
end

# Now, we define two different options using different `method` parameters, 
# one with `ConvolutionMethod` corresponding to the convolution,
# and the other with `NonHistoryMethod`, corresponding with the non-history method.
options_convolution = SimulationOptions(
    method = ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

options_nonhistory = SimulationOptions(
    method = NonHistoryMethod(),
    constraint = constraint,
    borefield = borefield,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

# Let us run the convolution
containers_convolution = @time initialize(options_convolution)
@time simulate!(operator=operator, options=options_convolution, containers=containers_convolution)

# And now let us run the non-history
containers_nonhistory = @time initialize(options_nonhistory)
@time simulate!(operator=operator, options=options_nonhistory, containers=containers_nonhistory)

# There is a massive speed up! We can also check that the two solutions are identical
sum(abs.(containers_convolution.X - containers_nonhistory.X))

# ## References
# [1] [Lazzarotto, Alberto; Basquens, Marc; Cimmino, Massimo; 
# _Non-history dependent temporal superposition algorithm for the pint source solution_,
# Research Conference Proceedings of the IGSHPA (2024).](https://doi.org/10.22488/okstate.24.000021)