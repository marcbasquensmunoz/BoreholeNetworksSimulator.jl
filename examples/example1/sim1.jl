    
using CSV, DataFrames, JLD2, Parameters
using WGLMakie, Colors, GeometryTypes
using Distances,LinearAlgebra, Dierckx
using DSP,Statistics

using BoreholeResponseFunctions

using BTESGroundWaterSimulator


cdir = @__DIR__

# 0. NAMES OF FOLDER FOR SAVING RESULTS FILES
symtitle = "1"
results_directory = "$cdir/results/sym$(symtitle)"
!isdir(results_directory) && mkdir(results_directory)


# 1. PROPERTIES OF BOREHOLE FIELD
ux_in_meterperday  = 1. * 1e-2      # groundwater speed along the flow coordinate
ux =  ux_in_meterperday/(3600*24)
λw = 0.6                            # water thermal conductivity
λs = 2.                             # ground thermal conductivity
Cw = 4.18*1e6                       # water thermal capacity
Cs = 1.7*1e6                        # ground thermal capacity
θ = 0.                              # angle of Darcy velocity
Φ = 0.2                             # porosity
λ = λs *(1-Φ) + λw*Φ                # porous medium conductivity
C = Cs *(1-Φ) + Cw*Φ                # porous medium capacity
α  = λ/C                            # porus medium thermal diffusivity
vt = ux * Cw/C                      # porous medium darcy velocity


# 2. GEOMETRICAL CONFIGURATION OF THE BOREHOLES ALONG THE VERTICAL DIRECTION

# Geometry on the vertical direction 

# H length of the borehole 
# h length of segment 
# D groundwater level
H,h,D = 50., 50., 4.


# 3. U-PIPE MODEL 
params = BoreholePara(λg=2.5,λs = λ)                                 # Define borehole u-pipe cross-section geometry 
upipe_crossection = GeometryTypes.Point2.([[0.03,0.0],[-0.03,.0]])   # (x,y) pipe positions on the borehole cross-section
@unpack rb,rpo = params                                              # borehole radius (rb) and pipe outer radius (rpo)

mf  =  0.5          # mass flow kg/s
ρf  =  1000.        # density 
cpf =  4182.        # specific heat capacity
Cf  =  ρf*cpf       # capacity  
Vf  =  mf/ρf        # volume flow

R = resistance_network(params, upipe_crossection)   # resistance network in the cross-section
A = coefficient_matrix(R, Cf, Vf)

# determine coefficients k_in, k_out and k_b for matrix build (Cimmino (2016))
k_in, k_out,  k_b = uniformTb_koeff(A,H)    # COEFFICIENTS OF THE BOREHOLE MODEL

# 4. BOREHOLE FIELD CONFIGURATION
#import configuration
cdir = @__DIR__
df = CSV.File("$cdir/data/Braedstrup_borehole_coordinates.txt"; decimal=',', delim = ';') |> DataFrame
#geometry of the field
borehole_positions =  [(x, y) for (x,y) in zip(df.X,df.Y) ]

# borehole connections map for extraction mode
external_to_internal = [ [22,30,37,29,36,35], 
                        [34,40,41,42,48,43],  
                        [33,32,39,45,46,47],                        
                        [26,25,24,31,38,44],  
                        [23,16,17,11,12,18],  
                        [10,5,6,7,13,19],     
                        [1,2,3,8,14,20],      
                        [4,9,15,21,28,27]                                         
                        ]

# borehole connections map for injection mode
internal_to_external = reverse.(external_to_internal)

# 5. DISCRETIZATION
z_ref = collect(D:h:D+H-h)          # line source reference point
z_eval = collect(D+h/2:h:D+H-h/2)   # evaluation points (evaluate at the mid point of the segment)

Nb = length(borehole_positions)     # number of boreholes
Nsb = length(z_eval)                # number of segments per borehole !
Ns  = Nb*Nsb                        # total number of segmensts

bh_map     =  reshape(ones(Int64, Nsb)*collect(1:Nb)',(Ns,1))
coord_source = [(x[1],x[2],p) for x in borehole_positions for p in z_ref] # position of sources 
coord_eval = [(x[1],x[2],p) for x in borehole_positions for p in z_eval]  # position of evaluation points

p =  GeometryTypes.Point3{Float64}.(coord_source)
tp = GeometryTypes.Point3{Float64}.(coord_eval)  

# # rotation of points in new coordinate system where Darcy velocity is parallel to x axis
p_rot  = rotation_z(p,-θ) 
tp_rot = rotation_z(tp,-θ) 

d = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp_rot) 
d = [d[1] == 0. && d[2] == 0. ?  (0.,params.rb, d[3],d[4]) : d for d in d]

tstep, tmax = 8760*3600/12., 8760*3600*10.
t = tstep:tstep:tmax
Nt = length(t) # number of time steps


# 6. THERMAL RESPONSES
# mutual response function between pairs of segments (adiabatic surface boundary condition)
g = [[1/(2π*params.λs)*mfls_adiabatic_surface(tt, α, coord[1:3]... ,vt, h, coord[4]; atol =1e-9) for tt in t] for coord in d]
# Matrix containing response function for each pair of segments at time-step 1
G = [g[1] for g in g]


# 7. LOADING CONDITION: for this particular simulation we impose temperature as boundary condition
T0 = 10.                                                 # undisturbed temperature
Tfin_constraint =  [i%12 in 1:6 ? 90. : 55. for i=1:Nt]  # input temperature
# Tfin_constraint = 90*ones(Nt)


# 8. MATRIX ASSEMBLY EXAMPLE
M_injection = zeros(3Nb+Ns,3Nb+Ns)      # matrix describing topology of injection problem
M_extraction = zeros(3Nb+Ns,3Nb+Ns)     # matrix describing topology of extraction problem
b = zeros(3Nb+Ns)                       # given term 
X = zeros(Nt,3Nb+Ns)                    # vector of unknowns
qprime = zeros(Nt, Ns)                  # heat injection per meter
Δqbcurrentsum = zeros(Ns)               # net heat injection on a given segment (this variable is needed by the solver)


build_matrix!(M_injection,Nb,Ns,
              k_in, k_out, k_b, 
              G,
              mf, cpf, bh_map, h,             
              internal_to_external
            )

build_matrix!(M_extraction,Nb,Ns,
            k_in, k_out, k_b, 
            G,
            mf, cpf, bh_map, h,             
            external_to_internal
          )

build_giventerm!(b,Nb,Ns,Tfin_constraint[1],T0,internal_to_external)



function solve_problem!(X,M_injection,M_extraction,
                        b,
                        Nb, Ns, Nt,
                        Tfin_constraint,
                        internal_to_external,external_to_internal,
                        qprime, g, T0,
                        Δqbcurrentsum, h                  
                        )

    for i=1:Nt
        M = i%12 in 1:6 ? M_injection : M_extraction
        branches = i%12 in 1:6 ? internal_to_external : external_to_internal
        # b = i%12 in 1:6 ? b_injection : b_extraction        

        solve_full_convolution_step!(X,M,b, i, Nb, Ns,
                    Tfin_constraint,branches,
                    qprime, g, T0,
                    Δqbcurrentsum, h                  
                    )
    end
end



@time solve_problem!(X,M_injection,M_extraction,
                b,
                Nb, Ns, Nt,
                Tfin_constraint,
                internal_to_external, external_to_internal,
                qprime, g, T0,
                Δqbcurrentsum, h                  
                )

# solve_problem!(X,M_injection,M_injection,
#                 b, 
#                 Nb, Ns, Nt,
#                 Tfin_constraint,
#                 internal_to_external, internal_to_external,
#                 qprime, g, T0,
#                 Δqbcurrentsum, h                  
#                 )

# solve_problem!(X,M_extraction,M_extraction,
#                 b,
#                 Nb, Ns, Nt,
#                 Tfin_constraint,
#                 external_to_internal,external_to_internal,
#                 qprime, g, T0,
#                 Δqbcurrentsum, h                  
#                 )

# 9. EXTRACT FROM SOLUTION VECTOR X
Tfin  = X[:, 1:2:2Nb]
Tfout = X[:, 2:2:2Nb]
Tb    = X[:, 2Nb+1:3Nb] 
q     =  cumsum(qprime, dims = 1) 

# output temperature to compute  energy and exergy echanged 
last_borehole_in_branch = [[i%12 in 1:6 ? x[1][end] : x[2][end] for i=1:Nt] for x in zip(internal_to_external, external_to_internal)]
Tfos  = hcat([[Tfout[idx,i] for (idx,i) in enumerate(ll)] for ll in last_borehole_in_branch]...)
Tfo   = mean(Tfos, dims =2)
Tfo   = reshape(Tfo, length(Tfo))


# 10. KPIs COMPUTATION
nbranches = 8
Energy_exchanged   = mf * cpf .* nbranches .* (Tfin_constraint .- Tfo) * (8760/12.)/1e6
Energy_exchanged_by_season = collect(reshape(Energy_exchanged, 6, div(Nt,6))')
Energy_exchanged_by_season2 = sum(Energy_exchanged_by_season , dims =2)

E_injected   =  Energy_exchanged_by_season2[isodd.(1:length(Energy_exchanged_by_season2))]
E_extracted  =  Energy_exchanged_by_season2[iseven.(1:length(Energy_exchanged_by_season2))]

eta = - E_extracted ./ E_injected

Exergy_exchanged = mf * cpf .* nbranches .* (Tfin_constraint .- Tfo  .- (T0 +273.15) .* log.( (Tfin_constraint .+ 273.15) ./ (Tfo .+ 273.15) )  ) * (8760/12.)/1e6
Exergy_exchanged_by_season = collect(reshape(Exergy_exchanged, 6, div(Nt,6))')
Exergy_exchanged_by_season2 = sum(Exergy_exchanged_by_season , dims =2)

Ex_injected   =  Exergy_exchanged_by_season2[isodd.(1:length(Exergy_exchanged_by_season2))]
Ex_extracted  =  Exergy_exchanged_by_season2[iseven.(1:length(Exergy_exchanged_by_season2))]

psi = - Ex_extracted ./ Ex_injected

##.

save("$(results_directory)/data_test$(symtitle).jld2" , 
        Dict( "ux" => ux,              
              "Tfin" => Tfin,
              "Tfin_constraint" => Tfin_constraint,
              "Tb" => Tb,
              "Tfout" => Tfout,
              "q"  => q,
              "qprime" => qprime,
              "borehole_positions" => borehole_positions,
              "E_injected" => E_injected,
              "E_extracted" => E_extracted,
              "eta" => eta,
              "psi" => psi
                )
)

##.
# 11. PLOTS 
set_theme!(
    font = "Arial", # inherited by layoutables if not overridden
    fontsize = 15, # inherited by layoutables if not overridden
    colgap = 20,
)


# FIGURE 1: borehole configuration
nbranch1, nbranch2 = 8, 3 
branch1 = external_to_internal[nbranch1]
branch2 = external_to_internal[nbranch2]
color_branch1 = range(colorant"blue", stop=colorant"darkorange", length=length(branch1));
color_branch2 = range(colorant"red", stop=colorant"green", length=length(branch2));


scene1 = Figure(resolution = (400, 400))
ax1 = scene1[1, 1] = Axis(scene1, ylabel = "y [m]", xlabel = "x[m]",aspect = DataAspect())
scatter!(ax1,  Makie.Point2.(borehole_positions) , markersize = 15.)

for i =1:length(borehole_positions)
    text!(ax1, "$i", fontsize = 9, position  = (borehole_positions[i] .+  (0.8 , -0.7)) , color = :blue )
end

for b in internal_to_external
    for k=2:length(b)
        s = Makie.Point.([ borehole_positions[b][k-1], borehole_positions[b][k] ] )
        linesegments!(ax1, s, color = :red, linewidth = 1.5 )
    end
end


for b in zip(branch1, color_branch1)
    local_p =  [ Makie.Point(borehole_positions[b[1]]) ]
    scatter!(ax1,  local_p , color = b[2], markersize = 24)
end

for b in zip(branch2, color_branch2)
    local_p =  [ Makie.Point(borehole_positions[b[1]]) ]
    scatter!(ax1,  local_p , color = b[2], markersize = 24)
end

scene1
##.

# FIGURE 2: temperature and heat flow along selected boreholes

scene12 = Figure(resolution = (450, 450))
ax2 = scene12[1, 1] = Axis(scene12, ylabel = "T [°C]") #, xlabel = "time [months]")
for b in zip(branch1, color_branch1)
    lines!(ax2, collect(t ./ 12tstep), Tfin[:,b[1]], color = b[2], linewidth = 2.)
    # scatter!(ax2, collect(t ./ tstep), Tfin[:,b[1]], color = b[2], markersize = 12., marker = '▲')
end
scatter!(ax2, collect(t ./ 12tstep), Tfos[:,nbranch1],  color = :black, markersize = 8.)
scatter!(ax2, collect(t ./ 12tstep), Tfo,  color = :grey, markersize = 12.)

ax3 = scene12[2, 1] = Axis(scene12, ylabel = "q [W/m]", xlabel = "time [years]")
for b in zip(branch1, color_branch1)    
    lines!(ax3, collect(t ./ 12tstep), q[:,b[1]], color = b[2], linewidth = 2.)
end

hidexdecorations!(ax2, grid = false)

scene12
##.
# FIGURE 3: temperature and heat flow along selected boreholes
scene13 = Figure(resolution = (450, 450))
ax4 = scene13[1, 1] = Axis(scene13, ylabel = "T [°C]") #, xlabel = "time [months]")
for b in zip(branch2, color_branch2)
    lines!(ax4, collect(t ./ 12tstep), Tfin[:,b[1]],  color = b[2], linewidth = 2.)    
end

scatter!(ax4, collect(t ./ 12tstep), Tfos[:,nbranch2],  color = :black, markersize = 12.)
scatter!(ax4, collect(t ./ 12tstep), Tfo,  color = :grey, markersize = 8.)

ax5 = scene13[2, 1] = Axis(scene13, ylabel = "q (W/m)", xlabel = "time [years]") #, xlabel = "time [months]")
for b in zip(branch2, color_branch2)
    lines!(ax5, collect(t ./ 12tstep), q[:,b[1]],  color = b[2], linewidth = 2.)
    
end

hidexdecorations!(ax4, grid = false)

scene13
##.
#############################
# COLOR MAP VISUALIZATION

# in order to visualize the simulation on a color map we need to use the obtained load distribution 
# and calculate the effect of such load on the grid

# Define the grid and compute the effect of the sources on the grid

x = -18:0.5:30
y = -18:0.5:18
test_points2 = [(x,y) for x in x for y in y]
tp2 = GeometryTypes.Point2{Float64}.(test_points2)  
coord2_eval = [(x[1],x[2],p) for x in test_points2 for p in z_eval[1]]
# rotation of points in new coordinate system where Darcy velocity is parallel to x axis 
tp2 = GeometryTypes.Point3{Float64}.(coord2_eval)  

tp2_rot = rotation_z(tp2,-θ) 
# evaluate unique tuples (Δx, |Δy|) among mutual points 
d2 = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp2_rot)
d2 = [d[1] == 0. && d[2] == 0. ?  (0.,params.rb, d[3],d[4]) : d for d in d2]


t_eval =  10 .^(collect(range( log10(tstep), stop = log10(tmax), length = 15)))

gg = [[1/(2π*params.λs)*mfls_adiabatic_surface(tt, α, coord[1:3]... ,vt, h, coord[4]; atol =1e-8) for tt in t_eval] for coord in d2]
gspline = [Spline1D(t_eval, g) for g in gg]     # interpolant
ggg = [g(t) for g in gspline]                   # evaluate the pairwise response function for all time-step in the simulation

# # matrix of vectors containing mutual contributions to the overall temperature build-up 
Ns, Ne = size(ggg)
ΔTji = [conv(qprime[:,j] , ggg[j,i])[1:Nt] for i=1:Ne,j=1:Ns]
# # summing the contribution of all sources on the segment i we can calculate its temperature
Ts_verification = sum(ΔTji, dims = 2) 

Tmax = maximum(vcat(Ts_verification...) .+ T0)


# FIGURE 4: heat map of the temperature field in time
# use the slider to visualize the evolution of the temperature field 

##.
# Slider figure
scene3 = Figure(size = (400, 400),  rowgap = 0.);

sl1 = scene3[2, 1] = Slider(scene3, range = 2:2:120, startvalue = 3, tellheight = true)

month_idx = @lift("month "*string($(sl1.value)))
Tt = @lift( [T[$(sl1.value)] .+ T0 for T in  Ts_verification] )
Z  = @lift( collect(reshape($(Tt), (length(y),length(x) ) )')  )

ax31 = scene3[1, 1] = Axis(scene3; aspect = DataAspect());
heatmap1 = heatmap!(ax31, x, y, Z, interpolate = true, colorrange = (10.,Tmax) )
scatter!(ax31, p, markersize = 10)

text!(ax31, month_idx,
      position = GeometryTypes.Point2(-15, 12.5),
      color= :white,
      align = (:left, :bottom), 
      fontsize = 15
     )

limits!(ax31, -17.5, 29.5, -17.5, 17.5)
hidedecorations!(ax31)

cbar = scene3[1, 2] = Colorbar(scene3, heatmap1, label = "T [°C]", width = 10)
cbar.height = Relative(2/3)
     
supertitle = scene3[0, :] = Label(scene3, "ground water speed $ux_in_meterperday m/day", fontsize = 18 ,tellheight = true)

scene3

##.


# FIGURE 5: snapshots of the temperature field
# six figures stacked
steps_of_interest = collect(2:2:12) .+ 12*9
Tts = [ [T[idx] .+ T0 for T in  Ts_verification] for idx in steps_of_interest]
Zs  = [collect(reshape(Tt, (length(y),length(x) ) )')   for Tt in Tts]

scene3

scene4 = Figure(resolution = (2700, 1300),  rowgap = 15, colgap =15 );
grid_structure = [(i,j) for i=1:2 for j =1:3]

for i=1:6
    ax = scene4[grid_structure[i]...] = Axis(scene4; aspect = DataAspect());
    heatmap = heatmap!(ax, x, y, Zs[i], interpolate = true, colorrange = (10.,Tmax) )
    scatter!(ax, p, markersize = 15)
    
    month_string = "month $(steps_of_interest[i])"
    text!(ax, month_string,
          position = GeometryTypes.Point2(-15, 12.5),
          color= :white,
          align = (:left, :bottom), 
          fontsize = 2.5
         )

    limits!(ax, -17.5, 29.5, -17.5, 17.5)
    hidedecorations!(ax)
end

cbar = scene4[1:end, end+1 ] = Colorbar(scene4, heatmap1, label = "T [°C]", width = 50)
cbar.height = Relative(2/3)
     
supertitle = scene4[0, :] = Label(scene4, "ground water speed $ux_in_meterperday m/day", fontsize = 50)
scene4

# Makie.save("$(results_directory)/configuration.png",scene1)
Makie.save("$(results_directory)/branch1_test$symtitle.png",scene12)
Makie.save("$(results_directory)/branch2_test$symtitle.png",scene13)
Makie.save("$(results_directory)/heatmap_test$symtitle.png",scene4)

