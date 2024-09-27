# Running BoreholeNetworksSimulator from Python

It is also possible to run BoreholeNetworksSimulator from Python by using its Julia interoperability (achieved with the package `PythonCall.jl` and the module `juliacall`) while keeping roughly the same syntax. For a more detailed explanation of how this works, please visit the [PythonCall.jl documentation](https://juliapy.github.io/PythonCall.jl/stable/).
This requires to have installed in your python environment the modules `juliacall` and `numpy`.

BoreholeNetworksSimulator.jl contains a second Julia module called `BNSPythonAdapter` whose purpose is to make the necessary conversions between Python and Julia. This module is not intended to be used from Julia, but from within Python.
In this tutorial, we will repeat the example in [Basic tutorial](@ref) but run from Python to show how it works. The full example as a Python script in at `BNSPythonAdapter/example.py`.

BNSPythonAdapter contains a python script `src/adapter.py` that executes the necessary bridging code.
You should start your python program by importing this module. Assuming that the directory `BNSPythonAdapter` is added to the `sys.path`, this should work:
````
import src.adapter
````
The script takes care of activating the Julia package `BNSPythonAdapter` which has `BoreholeNetworksSimulator` as a dependency, and loads both namespaces for later use.
It also defines conversions from Python objects to their Julia counterparts. 
Next, we need to import `juliacall` to be able to reference objects and funcions from `BoreholeNetworksSimulator`. It is useful to give it an alias for conciseness.
````
from juliacall import Main as jl
````
Now, the python variable `jl` represents the `Main` module in Julia, and it has as methods all the Julia functions available in it. Since we imported the namespace `BoreholeNetworksSimulator`, this also includes its objects and functions.

Then, we define the same variables as in [Basic tutorial](@ref):
````
Δt = 8760*3600/12.
Nt = 10*12

α = 1e-6
λ = 3.
T0 = 10.
medium = jl.GroundMedium(α=α, λ=λ, T0=T0)

D = 10.
H = 100.

borehole = jl.SingleUPipeBorehole(H=H, D=D)

σ = 5.
positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(0., 0.), (0., σ)])
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)

configurations = [jl.BoreholeNetwork(jl.Vector[jl.Vector[jl.Int]]([[1], [2]]))]

q1 = 5.
q2 = 5.
loads = jl.Array[jl.Float64]([q1, q2])
constraint = jl.constant_HeatLoadConstraint(loads, Nt)

options = jl.SimulationOptions(
    method = jl.ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    medium = medium,
    fluid = jl.Water(),
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)
````
Note that the code itself is not very different from its Julia version, but there are two remarks worth making. First, we need to call any object defined in Julia by typing `jl.` in front. This creates a python object with the same fields that `PythonCall.jl` knows how to convert back into a Julia object.
Second, note that we have defined the arrays by explicitly declaring their generic type. If we don't do this, they will be converted into `Vector{Any}` in the Julia code, which is not desirable.

Another difference is the definition of our `[Operator](@ref)` object. Since we are writing python, it should now be a python object with a method `operate` that returns an object of type `jl.BoreholeOperation`, which `PythonCall.jl` knows how to convert to `BoreholeOperation`. For the ssake of the example, we will implement again the `[SimpleOperator](@ref)` operator that is already implemented in Julia. For further examples, check the python version of the code for `[Example: Braedstrup borefield ](@ref)`.
````
class SimpleOperator():
    def __init__(self, mass_flow, branches):
        self.mass_flows = jl.Vector[jl.Float64](mass_flow * np.ones(branches))

    def operate(self, i, options, Tfin, Tfout, Tb, q):
        return jl.BoreholeOperation(options.configurations[0], operator.mass_flows)


operator = SimpleOperator(2., 2)
````

Note that this is equivalent to calling the same code from Julia:
````
operator = jl.SimpleOperator(mass_flow = 2., branches = 2)
````


Finally we can run the simulation
````
containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)
````
Note that the bang in `simulate!` is replaced by a `_b` in the juliacall version of the function. 
````
containers.X
````
Now, one can compare this result with its Julia version. If everything went right, they should agree!
