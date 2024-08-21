# Public API

## General 

These are the functions used to run the simulation. Note that an object of type `SimulationOptions` is needed to call them.

```@docs
initialize
```

```@docs
simulate!
```

Note that an array of `BoreholeNetwork` containing all the configurations allowed during the simulation must be specified in `SimulationOptions`.
On the other hand, `BoreholeOperation` is an object that needs to be returns by an `operator` function at each time step , representing the dynamical changes in the operation. See [Basic tutorial](@ref) for more details.

```@docs
BoreholeNetwork
```

```@docs
BoreholeOperation
```


## Simulation Options

The available options for the simulation are specified through a `SimulationOptions` object. This needs to be passed to the `initialize` and `simulate!` functions.
The options are modular: each particular option can be chosen independently of the others, allowing for a wide range of possible simulations. 
Note that in some particular cases there might be some incompatibilities between options or non-implemented interactions. In those cases, an error will appear explaining the reason.

```@docs
SimulationOptions
```

The several options are listed below:

### Medium

Models the underground medium through which the heat will transfer between boreholes.

```@docs
Medium
```

#### Options

```@docs
GroundMedium
```

```@docs
FlowInPorousMedium
```

### Borefield

Models the geometry of the borefield.

```@docs
Borefield
```

#### Options

```@docs
EqualBoreholesBorefield
```

### Borehole

Models the internal heat transfer in the borehole.

```@docs
Borehole
```

#### Options

```@docs
SingleUPipeBorehole
```

### Constraint

Imposes the working conditions and demands of the whole system.

```@docs
Constraint
```

#### Options

```@docs
HeatLoadConstraint
```

```@docs
InletTempConstraint
```

### Time Superposition Method

Applies methods for time superposition.

```@docs
TimeSuperpositionMethod
```

#### Options

```@docs
ConvolutionMethod
```

```@docs
NonHistoryMethod
```

### Boundary Condition

Models the ground surface.

```@docs
BoundaryCondition
```

#### Options
```@docs
NoBoundary
```

```@docs
DirichletBoundaryCondition
```
