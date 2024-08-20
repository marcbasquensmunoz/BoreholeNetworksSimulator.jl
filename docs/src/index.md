BoreholeNetworksSimulator.jl is a pure Julia, performant, and modular framework for simulations of fields of boreholes.
Features:
- Computes fluid temperatures, borehole wall temperatures and heat extracted. 
- Highly modular: different boreholes, hydraulic configurations, ground properties, load demand or temperature  constraints, ground boundary conditions can be seamlessly used with minimal effort. 
- Allows for design of operation strategies via an operator callback at each time step.
- Implements the "non-history" time superposition method, reducing the computational complexity in the number of time steps to linear. This allows for simulations with fine time steps.
