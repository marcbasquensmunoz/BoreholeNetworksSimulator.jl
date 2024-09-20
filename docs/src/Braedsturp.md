# Example: Braedsturp borefield 

The example considers a borehole field of 48 boreholes connected according to scheme utilized in the installation in Braedsturp, Denmark. 
The source code of this example can be found in `examples/Braedsturp/main.jl`.
The field is in a porous medium with an underground water flow of ``0.01 \frac{m}{\text{day}}``, which we assume to be in the ``x`` direction. 
The borefield consists of 8 branches, each branch in parallel, and within each branch there are ``6``boreholes in series. 
The positions and connections are represented in the following plot: 

![](./examples/braedsturp/configuration.png)

## Run the example
```
include("examples/Braedsturp/main.jl")
```
In this example, the simulation is executed and the results are saved in the cache file named 'cache_3.1536e8.jld2'. The result data is also available in memory:
```
julia> containers.X
192×120 Matrix{Float64}:
  70.8847   75.0666   77.993    79.8948   81.1608  …  55.0      55.0      55.0      55.0
  67.5464   71.9367   75.0863   77.1524   78.5424     54.0451   53.9479   53.8829   53.8375
  74.3405   78.1643   80.7556   82.4294   83.5353     54.0451   53.9479   53.8829   53.8375
  70.8847   75.0666   77.993    79.8948   81.1608     53.3751   53.162    53.0222   52.9267
  78.015    81.4475   83.6781   85.1081   86.0436     53.3751   53.162    53.0222   52.9267
   ⋮                                               ⋱                                
 166.095   -45.855   -26.2117  -19.1312  -14.2904      8.85064   6.64457   5.00616   3.79414
 176.965   -50.0092  -28.9111  -21.2066  -15.7431      8.59637   6.57639   4.9994    3.79898
 187.948   -55.3887  -30.2105  -21.656   -15.8568  …   7.91946   6.12587   4.71477   3.6243
 179.911   -47.192   -24.6851  -18.3117  -13.8542      6.52704   5.39342   4.38165   3.54095
```

## Plots
Plots showing the results of the simulation as also available
```
julia > include("examples/old/sim1.jl") 
```

Inlet borehole temperatures and heat flows for boreholes along two branches in the borehole field. The time series are color coded according to the previous configuration plot above. In addition to the inlet temperature, the output temperature from the branch (grey dot), and the mean output temperature from the field (black dot) are displayed.

![](./examples/braedsturp/branch1_test1.png)
![](./examples/braedsturp/branch2_test1.png)

Finally we can display the heatmap of the temperature field in the borehole region during the 10th year of operation

![](./examples/braedsturp/heatmap_test1.png)


## Running the code in Python
The example is also avaiable from Python. For details of how this is done refer to [Running BoreholeNetworksSimulator from Python](@ref). 
The Python version is in:
```
examples/Braedsturp/main.py
```
