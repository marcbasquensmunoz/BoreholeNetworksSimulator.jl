
# BoreholeNetworksSimulator.jl

BoreholeNetworksSimulator.jl is a pure [Julia](https://julialang.org/), performant, and modular framework for simulations of fields of interconnected borehole heat exchangers.
Features:
- Computes fluid temperatures, borehole wall temperatures and heat extracted. 
- Supports many different configurations and settings by being highly modular: boreholes, hydraulic configurations, ground properties, load demand or temperature  constraints, ground boundary conditions can be seamlessly used with minimal effort. 
- Allows for design of operation strategies via an operator callback at each time step.
- Implements the "non-history" time superposition method, reducing the computational complexity in the number of time steps to linear. This allows for simulations with fine time steps.
- Python interoperability.

# Getting started

BoreholeNetworksSimulator.jl is currently not in Julia's General registry (as well as some of its dependencies), however, it is easily available through the local registry [geothermal_registry](https://github.com/marcbasquensmunoz/geothermal_registry). 
In order to install it, start Julia and run the command:

````
using Pkg; pkg"registry add https://github.com/marcbasquensmunoz/geothermal_registry"; Pkg.add("BoreholeNetworksSimulator")
````

In order to learn how to use this package, please visit [Basic tutorial](@ref) first.

For users interested in running BoreholeNetworksSimulator.jl from Python, see also [Running BoreholeNetworksSimulator from Python](@ref).

# Design philosophy

The driving motivation when developing BoreholeNetworksSimulator.jl was to create an easy to use tool that is also highly flexible and performant.

This goal has been accomplished by leveraging Julia's Multiple Dispatch. Through a simple, common interface, the user can run a plethora of different simulations by modularly specifying options.
Some options have physical meaning or interpretation, while others are purely algorithmic, yet they can all be seamlessly changed by changing a parameter.
Please, refer to [Public API](@ref) for the full list of options offered by this package.
