from juliacall import Main as jl
from juliacall import Pkg as jlPkg
import pathlib

# path = jl.Base.find_package("BoreholeNetworksSimulator")

project_dir = str(pathlib.Path(__file__).parent.parent.resolve().absolute())

jlPkg.activate(jl.pyconvert(jl.String, project_dir))
# jlPkg.activate("C:/BoreholeNetworksSimulator.jl") 

# jlPkg.Registry.add("https://github.com/marcbasquensmunoz/geothermal_registry")
# jlPkg.update()
jlPkg.resolve()
jlPkg.instantiate()
jl.seval("using BNSPythonAdapter")
jl.seval("using BoreholeNetworksSimulator")

import numpy as np

class BoreholeNetwork:
    branches: np.array

    def __init__(self, branches):
        self.branches = branches


class BoreholeOperation:
    network: BoreholeNetwork
    mass_flows: np.array

    def __init__(self, network, mass_flows):
        self.network = network
        self.mass_flows = mass_flows



jl.seval('''
    function convertBoreholeNetwork(::Type{BoreholeNetwork}, x)
        branches = pyconvert(Vector{Int}, PyArray(x.branches))
        PythonCall.pyconvert_return(BoreholeNetwork(branches=branches))
    end
''')
jl.seval('''
    function convertBoreholeOperation(::Type{BoreholeOperation}, x)
        network = pyconvert(BoreholeNetwork, PyArray(x.network))
        mass_flows = pyconvert(Vector{Float64}, PyArray(x.mass_flows))
        PythonCall.pyconvert_return(BoreholeOperation(network=network, mass_flows=mass_flows))
    end
''')
jl.PythonCall.pyconvert_add_rule("__main__:BoreholeNetwork", jl.BoreholeNetwork, jl.convertBoreholeNetwork)
jl.PythonCall.pyconvert_add_rule("__main__:BoreholeOperation", jl.BoreholeOperation, jl.convertBoreholeOperation)

jl.seval("@info \"Load successful\"")