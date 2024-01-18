using CSV, DataFrames, JLD2, Parameters
    
using WGLMakie, Colors, GeometryTypes
using Distances,LinearAlgebra, Dierckx
using DSP,Statistics

using BoreholeResponseFunctions
using BTESGroundWaterSimulator

using PythonCall

@with_kw struct BorefieldProperties
    ux_in_meterperday      # groundwater speed along the flow coordinate
    λw                     # water thermal conductivity
    λs                     # ground thermal conductivity
    Cw                     # water thermal capacity
    Cs                     # ground thermal capacity
    θ                      # angle of Darcy velocity
    Φ  
end

@with_kw struct BoreholeProperties
    H       # length of the borehole   
    h       # length of segment 
    D       # groundwater level
    λg
    mf      # mass flow kg/s
    ρf      # density 
    cpf     # specific heat capacity
end

function sim1(;operator, borefield::BorefieldProperties, borehole::BoreholeProperties, configurations, tstep, tmax)

    cdir = @__DIR__

    # 0. NAMES OF FOLDER FOR SAVING RESULTS FILES
    symtitle = "1"
    results_directory = "$cdir/results/sym$(symtitle)"
    !isdir(results_directory) && mkdir(results_directory)

    # 1. PROPERTIES OF BOREHOLE FIELD
    @unpack D, h, H, λg, mf, ρf, cpf = borehole
    @unpack ux_in_meterperday, λw, λs, Cw, Cs ,θ, Φ = borefield

    ux = ux_in_meterperday/(3600*24)
    λ = λs *(1-Φ) + λw*Φ        # porous medium conductivity
    C = Cs *(1-Φ) + Cw*Φ        # porous medium capacity
    α  = λ/C                                                            # porus medium thermal diffusivity
    vt = ux * Cw/C                                            # porous medium darcy velocity

    # 3. U-PIPE MODEL 
    params = BoreholePara(λg = λg, λs = λ)                               # Define borehole u-pipe cross-section geometry 
    upipe_crossection = GeometryTypes.Point2.([[0.03,0.0],[-0.03,.0]])   # (x,y) pipe positions on the borehole cross-section
    @unpack rb,rpo = params                                              # borehole radius (rb) and pipe outer radius (rpo)

    Cf = ρf * cpf       # capacity  
    Vf = mf / ρf        # volume flow

    R = resistance_network(params, upipe_crossection)   # resistance network in the cross-section
    A = coefficient_matrix(R, Cf, Vf)

    # determine coefficients k_in, k_out and k_b for matrix build (Cimmino (2016))
    k_in, k_out, k_b = uniformTb_koeff(A, H)    # COEFFICIENTS OF THE BOREHOLE MODEL

    # 4. BOREHOLE FIELD CONFIGURATION
    #import configuration
    cdir = @__DIR__
    df = CSV.File("$cdir/data/Braedstrup_borehole_coordinates.txt"; decimal=',', delim = ';') |> DataFrame
    #geometry of the field
    borehole_positions = [(x,y) for (x,y) in zip(df.X,df.Y)]

    # 5. DISCRETIZATION
    z_ref  = collect(D:h:D+H-h)         # line source reference point
    z_eval = collect(D+h/2:h:D+H-h/2)   # evaluation points (evaluate at the mid point of the segment)

    Nb  = length(borehole_positions)    # number of boreholes
    Nsb = length(z_eval)                # number of segments per borehole !
    Ns  = Nb*Nsb                        # total number of segmensts
    Nc  = size(configurations)[1]       # number of configurations

    bh_map       = reshape(ones(Int64, Nsb)*collect(1:Nb)',(Ns,1))
    coord_source = [(x[1],x[2],p) for x in borehole_positions for p in z_ref]  # position of sources 
    coord_eval   = [(x[1],x[2],p) for x in borehole_positions for p in z_eval] # position of evaluation points

    p  = GeometryTypes.Point3{Float64}.(coord_source)
    tp = GeometryTypes.Point3{Float64}.(coord_eval)  

    # # rotation of points in new coordinate system where Darcy velocity is parallel to x axis
    p_rot  = rotation_z(p,-θ) 
    tp_rot = rotation_z(tp,-θ) 

    d = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp_rot) 
    d = [d[1] == 0. && d[2] == 0. ?  (0.,params.rb, d[3],d[4]) : d for d in d]

    t = tstep:tstep:tmax
    Nt = length(t) # number of time steps

    # 6. THERMAL RESPONSES
    # mutual response function between pairs of segments (adiabatic surface boundary condition)
    g = [1/(2π*params.λs) * mfls_adiabatic_surface(tt, α, coord[1:3]..., vt, h, coord[4], atol = 1e-9) for coord in d, tt in t]
    # Matrix containing response function for each pair of segments at time-step 1
    G = g[:, :, 1]

    # 7. LOADING CONDITION: for this particular simulation we impose temperature as boundary condition
    T0 = 10.                                                 # undisturbed temperature
    Tfin_constraint =  [i%12 in 1:6 ? 90. : 55. for i=1:Nt]  # input temperature

    # 8. MATRIX ASSEMBLY EXAMPLE
    M = zeros(3Nb+Ns, 3Nb+Ns, Nc)           # matrix describing topology
    b = zeros(3Nb+Ns)                       # given term 
    X = zeros(Nt,3Nb+Ns)                    # vector of unknowns
    qprime = zeros(Nt, Ns)                  # heat injection per meter
    Δqbcurrentsum = zeros(Ns)               # net heat injection on a given segment (this variable is needed by the solver)

    # Build matrix M
    for i in 1:Nc
        @views build_matrix!(
            M[:,:,i], Nb, Ns,
            k_in, k_out, k_b, 
            G,
            mf, cpf, bh_map, h,             
            configurations[i]
        )
    end

    build_giventerm!(b, Nb, Ns, Tfin_constraint[1], T0, configurations[1])
    last_borehole_in_branch = [zeros(Int, Nt) for i in eachindex(configurations[1])]

    # Solve problem by iterating
    for i=1:Nt
        topology = pyconvert(Int, operator(i, X, qprime))
        println(topology)
        @views solve_full_convolution_step!(X, 
                    M[:, :, topology], b, i, Nb, Ns,
                    Tfin_constraint, configurations[topology],
                    qprime, g, T0,
                    Δqbcurrentsum, h                 
                    )
        for j in eachindex(configurations[1])
            last_borehole_in_branch[j][i] = last(configurations[topology][j])
        end
    end

    # 9. EXTRACT FROM SOLUTION VECTOR X
    Tfin  = X[:, 1:2:2Nb]
    Tfout = X[:, 2:2:2Nb]
    Tb    = X[:, 2Nb+1:3Nb] 
    q     = cumsum(qprime, dims = 1) 

    # output temperature to compute energy and exergy echanged 
    Tfos = hcat([[Tfout[idx,i] for (idx,i) in enumerate(ll)] for ll in last_borehole_in_branch]...)
    Tfo  = vec(mean(Tfos, dims = 2))

    # 10. KPIs COMPUTATION
    to_K = 273.15

    nbranches = size(configurations[1], 1)
    Energy_exchanged = @. mf * cpf * nbranches * (Tfin_constraint - Tfo) * (8760/12.)/1e6
    Energy_exchanged_by_season = collect(reshape(Energy_exchanged, 6, div(Nt,6))')
    Energy_exchanged_by_season2 = sum(Energy_exchanged_by_season , dims = 2)

    E_injected   =  Energy_exchanged_by_season2[1:2:length(Energy_exchanged_by_season2)]
    E_extracted  =  Energy_exchanged_by_season2[2:2:length(Energy_exchanged_by_season2)]

    eta = - E_extracted ./ E_injected

    Exergy_exchanged = @. mf * cpf * nbranches * (Tfin_constraint - Tfo  - (T0 + to_K) * log( (Tfin_constraint + to_K) / (Tfo + to_K) )  ) * (8760/12.)/1e6
    Exergy_exchanged_by_season = collect(reshape(Exergy_exchanged, 6, div(Nt,6))')
    Exergy_exchanged_by_season2 = sum(Exergy_exchanged_by_season , dims = 2)

    Ex_injected   =  Exergy_exchanged_by_season2[1:2:length(Exergy_exchanged_by_season2)]
    Ex_extracted  =  Exergy_exchanged_by_season2[2:2:length(Exergy_exchanged_by_season2)]

    psi = - Ex_extracted ./ Ex_injected

    save("$(results_directory)/data_test$(symtitle).jld2" , 
        Dict( 
            "ux" => ux,              
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
    return eta
end


function seasonal_operator(i, X, q)
    i%12 in 1:6 ? 2 : 1
end

borefield = BorefieldProperties(ux_in_meterperday = 1e-2, λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2)
borehole = BoreholeProperties(H = 50., h = 50., D = 4., λg = 2.5, mf = 0.5, ρf = 1000., cpf = 4182.)
configurations = 
[
    [   
        [22,30,37,29,36,35], 
        [34,40,41,42,48,43],  
        [33,32,39,45,46,47],                        
        [26,25,24,31,38,44],  
        [23,16,17,11,12,18],  
        [10,5,6,7,13,19],     
        [1,2,3,8,14,20],      
        [4,9,15,21,28,27]                                         
    ]
    ,
    [
        [35,36,29,37,30,22],
        [43,48,42,41,40,34],
        [47,46,45,39,32,33],
        [44,38,31,24,25,26],
        [18,12,11,17,16,23],
        [19,13,7,6,5,10],
        [20,14,8,3,2,1],
        [27,28,21,15,9,4]
    ]
]
tstep = 8760*3600/12.
tmax  = 8760*3600*10.

sim1(operator=seasonal_operator, borefield=borefield, borehole=borehole, configurations=configurations, tstep=tstep, tmax=tmax)
