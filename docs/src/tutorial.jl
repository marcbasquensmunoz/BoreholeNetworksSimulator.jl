# # Basic tutorial
# In this tutorial we will learn how to simulate the expected temperature and heat extraction in a borehole 
# field.
# 
# First, we need to specify the parameters and constraints of our system, as well as some options for the simulation. 
# This is done through a [`SimulationOptions`](@ref) object. Its variables are modular components 
# with which we can describe many scenarios by leveraging Julia's multiple dispatch. 
# Let us see the components one by one with an example.
#
# We start by specifying the simulation time step and the simulation duration. For our example, 
# we will take monthly time steps during 10 years:

using BoreholeNetworksSimulator

Δt = 8760*3600/12.
Nt = 10*12

# Suppose the ground we are interested in simulating is made of solid rock. This means the heat
# transfer will occur by pure conduction. Assume that the thermal diffusivity of the rock is ``α = 10^{-6} \frac{m^2}{s}``
# and its thermal conductivity is ``λ = 3 \frac{W}{m \ K}``. 
# The undisturbed temperature of the ground is ``T_0=10 \ ^{\circ}C``.
# We model the ground with a subtype of [`Medium`](@ref), in 
# our case, as per our assumptions, we are particularly interested in [`GroundMedium`](@ref):

@show GroundMedium()
α = 1e-6
λ = 3.
T0 = 10.
medium = GroundMedium(α=α, λ=λ, T0=T0)


D = 10.
H = 100.
borehole = SingleUPipeBorehole(H=H, D=D)

# Next, we need to know how and where the boreholes in our system are. Suppose we are interested
# in simulating the evolution of two identical vertical boreholes of length ``H=100m`` and 
# buried at depth ``D=10m``, separated by a distance of ``σ=5m``. The boreholes are of the U-type.
# We need to enconde this information
# in a subtype of [`Borefield`](@ref). Since the boreholes in our example have the same properties,
# we can use [`EqualBoreholesBorefield`](@ref), which instantiates several identical boreholes 
# from a prototype. The prototype is specified by a subtype of [`Borehole`](@ref). 
# In our case, we can use [`SingleUPipeBorehole`](@ref) to model a borehole with a single U-pipe.

D = 10.
H = 100.

borehole = SingleUPipeBorehole(H=H, D=D)

# There are more parameters that are also relevant, such as borehole radius, grout properties,
# pipe resistance, etc., but for the moment we will use their default values. 
#
# Next, we need to specify where the borehole are located.

σ = 5.
positions = [(0., 0.), (0., σ)]
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)

# Note that we didn't specify yet how the boreholes are connected. This is because the simulation 
# allows for dynamical changes in the network hydraulic configuration. Using different configurations,
# we can simulate the reverse flow, or even simulate valves opening and closing.
# The reason behind this design is to allow for decision making on the operation of the borefield
# _during_ the simulation, based on the inputs that we are interested in.
# 
# For now, let us state the possible hydraulic configurations that we will allows our borefield to be in.
# This is done with an vector of [`BoreholeNetwork`](@ref). Borefields usually have several branches
# of boreholes. Each borehole in a branch is connected in series, while branches may o may not be connected 
# in parallel between themselves. To create a [`BoreholeNetwork`](@ref), we need to provide
# all of its branches in a vector. Each branch is, in turn, a vector containing the identifiers, in order, 
# of the boreholes present in that branch. Each identifier `i::Int` refers to the borehole located at position
# `positions[i]`.
#
# In our example, we want to simulate two independent boreholes, so each of them must be in a separate branch.
# Also, for the moment, we are only interested in this configuration, so we define:

configurations = [BoreholeNetwork([[1], [2]])]

# Even with all these specifications, the evolution of the system is still not fully determined.
# The missing conditions are referred to as constraints, and are modeled by subtypes of [`Constraint`](@ref).
# For instance, if we would like the two boreholes to be connected in parallel, we would still need to 
# impose that their inlet temperatures be equal. In our example, since we want out boreholes to be independent, 
# we will impose the total amount of heat that we want to extract from each borehole. We will impose a constant
# load, equal for both boreholes. This is specified by
q1 = 5.
q2 = 5.
constraint = constant_HeatLoadConstraint([q1, q2], Nt)

# We can finally create the object with all the options:

options = SimulationOptions(
    method = ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

# As we have mentioned, the simulation is designed to allow for a controllable opeartion during its duration.
# We do this by defining a function that takes as an input the current state of the borefield and outputs
# a [`BoreholeOperation`](@ref) object. This object has two variables: the first specifies which 
# configuration will be used for the next time step. In our case, we only want a single configuration.
# The second specifies the mass flow rate through each branch of the selected configuration, provided as 
# a vector. In our example, we will keep this constant through the simulation:

function operator(i, Tin, Tout, Tb, q, configurations)
    network = configurations[1]
    BoreholeOperation(network, 2 .* ones(n_branches(network)))
end

# Before simulating, we first need to call [`initialize`](@ref) to run some precomputations
# that will be used throught the simulation and to instantiate containers where the result will be written.
containers = initialize(options)

# And finally, we can start the simulation.
@time simulate!(operator=operator, options=options, containers=containers)

# The simulation is over! Note that the bulk of the time is spent in the precoputation, while the 
# simulation itself is quite fast. This is especially good if we want to test different operation 
# strategies.
#
# The result is saved in
containers.X
