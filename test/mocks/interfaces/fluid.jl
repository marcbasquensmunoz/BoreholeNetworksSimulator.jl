
"""
    FluidMock <: Fluid 

Mock for testing purposes.
"""
@with_kw struct FluidMock <: Fluid 
    cpf = 0.
    properties = (0., 0., 0., 0.)
end
BoreholeNetworksSimulator.cpf(fluid::FluidMock) = fluid.cpf
BoreholeNetworksSimulator.thermophysical_properties(fluid::FluidMock, T) = fluid.properties
