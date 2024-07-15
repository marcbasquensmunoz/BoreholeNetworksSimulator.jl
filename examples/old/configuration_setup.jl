using DataFrames, CSV
using Makie, Colors
using WGLMakie

cdir = @__DIR__
df = CSV.File("$cdir/data/Braedstrup_borehole_coordinates.txt"; decimal=',', delim = ';') |> DataFrame
#geometry of the field
borehole_positions =  [(x, y) for (x,y) in zip(df.X,df.Y) ]

# connection pattern for external_to_internal flow 
external_to_internal = [ [22,30,37,29,36,35],
                        [43,48,42,41,40,34],
                        [47,46,45,39,32,33],
                        [44,38,31,24,25,26],
                        [23,16,17,11,12,18],
                        [10,5,6,7,13,19],
                        [1,2,3,8,14,20],
                        [4,9,15,21,28,27]                                        
                        ]

# connection pattern for internal_to_external flow
internal_to_external = reverse.(external_to_internal)



# VISUALIZATION OF THE BOREHOLE FIELD

# draw boreholes
# scene1, layout = fig(resolution = (1000, 1000))

ff1 = Figure(size = (500, 500))
ax1 = ff1[1, 1] = Axis(ff1, ylabel = "y [m]", xlabel = "x[m]",aspect = DataAspect())
scatter!(ax1,  Makie.Point2.(borehole_positions) , markersize = 15.)

# numbering boreholes
for i =1:length(borehole_positions)
    text!(ax1, "$i", fontsize = 10., position  = (borehole_positions[i] .+  (0.7 , -0.3)) , color = :blue )
end

# visualize connections
for b in internal_to_external
    for k=2:length(b)
        s = Makie.Point.([ borehole_positions[b][k-1], borehole_positions[b][k] ] )
        linesegments!(ax1, s, color = :red, linewidth = 1.5)
    end
end

# color code boreholes
nbranch1, nbranch2 = 8, 3 
branch1 = external_to_internal[nbranch1]
branch2 = external_to_internal[nbranch2]
color_branch1 = range(colorant"blue", stop=colorant"darkorange", length=length(branch1));
color_branch2 = range(colorant"red", stop=colorant"green", length=length(branch2));

for b in zip(branch1, color_branch1)
    local_p =  [ Makie.Point(borehole_positions[b[1]]) ]
    scatter!(ax1,  local_p , color = b[2], markersize = 24)
end

for b in zip(branch2, color_branch2)
    local_p =  [ Makie.Point(borehole_positions[b[1]]) ]
    scatter!(ax1,  local_p , color = b[2], markersize = 24)
end


# create arrows to visualize ground water velocity
p1 = [Makie.Point(-15.2,y) for y = -8:3:8]
v1 = fill(Makie.Vec(2.,0), length(p1))
text!(ax1, "ux", fontsize = 1.5, position  = Makie.Point(-15. , -1.2) , color = :dodgerblue )

arrows!(ax1, p1, v1, linecolor =:dodgerblue, arrowcolor =:dodgerblue, linewidth = 4., arrowsize = 15)

limits!(ax1, -17.5, 13.5, -15.5, 13.5)

# Makie.save("$cdir/results/configuration.png",ff1)

ff1