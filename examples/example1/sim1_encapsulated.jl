    
using CSV, DataFrames, JLD2, Parameters
using WGLMakie, Colors, GeometryTypes
using Distances,LinearAlgebra, Dierckx
using DSP,Statistics

using BoreholeResponseFunctions

using BTESGroundWaterSimulator

struct BorefieldProperties
    ux_in_meterperday      # groundwater speed along the flow coordinate
    λw                     # water thermal conductivity
    λs                     # ground thermal conductivity
    Cw                     # water thermal capacity
    Cs                     # ground thermal capacity
    θ                      # angle of Darcy velocity
    Φ  
end

struct BoreholeProperties
    H   # length of the borehole   
    h   # length of segment 
    D   # groundwater level
    λg
    mf      # mass flow kg/s
    ρf      # density 
    cpf     # specific heat capacity
end

function sim1(borefield::BorefieldProperties, borehole::BoreholeProperties, external_to_internal; tstep, tmax)

    cdir = @__DIR__

    # 0. NAMES OF FOLDER FOR SAVING RESULTS FILES
    symtitle = "1"
    results_directory = "$cdir/results/sym$(symtitle)"
    !isdir(results_directory) && mkdir(results_directory)


    # 1. PROPERTIES OF BOREHOLE FIELD
    ux = borefield.ux_in_meterperday/(3600*24)
    λ = borefield.λs *(1-borefield.Φ) + borefield.λw*borefield.Φ        # porous medium conductivity
    C = borefield.Cs *(1-borefield.Φ) + borefield.Cw*borefield.Φ        # porous medium capacity
    α  = λ/C                                                            # porus medium thermal diffusivity
    vt = ux * borefield.Cw/C                                            # porous medium darcy velocity


    # 2. GEOMETRICAL CONFIGURATION OF THE BOREHOLES ALONG THE VERTICAL DIRECTION

    # Geometry on the vertical direction 

    # 3. U-PIPE MODEL 
    params = BoreholePara(λg=borehole.λg,λs = λ)                                 # Define borehole u-pipe cross-section geometry 
    upipe_crossection = GeometryTypes.Point2.([[0.03,0.0],[-0.03,.0]])   # (x,y) pipe positions on the borehole cross-section
    @unpack rb,rpo = params                                              # borehole radius (rb) and pipe outer radius (rpo)

   
    Cf  =  borehole.ρf*borehole.cpf       # capacity  
    Vf  =  borehole.mf/borehole.ρf        # volume flow

    R = resistance_network(params, upipe_crossection)   # resistance network in the cross-section
    A = coefficient_matrix(R, Cf, Vf)

    # determine coefficients k_in, k_out and k_b for matrix build (Cimmino (2016))
    k_in, k_out, k_b = uniformTb_koeff(A,borehole.H)    # COEFFICIENTS OF THE BOREHOLE MODEL

    # 4. BOREHOLE FIELD CONFIGURATION
    #import configuration
    cdir = @__DIR__
    df = CSV.File("$cdir/data/Braedstrup_borehole_coordinates.txt"; decimal=',', delim = ';') |> DataFrame
    #geometry of the field
    borehole_positions =  [(x, y) for (x,y) in zip(df.X,df.Y) ]

    # borehole connections map for injection mode
    internal_to_external = reverse.(external_to_internal)

    # 5. DISCRETIZATION
    z_ref = collect(borehole.D:borehole.h:borehole.D+borehole.H-borehole.h)          # line source reference point
    z_eval = collect(borehole.D+borehole.h/2:borehole.h:borehole.D+borehole.H-borehole.h/2)   # evaluation points (evaluate at the mid point of the segment)

    Nb = length(borehole_positions)     # number of boreholes
    Nsb = length(z_eval)                # number of segments per borehole !
    Ns  = Nb*Nsb                        # total number of segmensts

    bh_map     =  reshape(ones(Int64, Nsb)*collect(1:Nb)',(Ns,1))
    coord_source = [(x[1],x[2],p) for x in borehole_positions for p in z_ref] # position of sources 
    coord_eval = [(x[1],x[2],p) for x in borehole_positions for p in z_eval]  # position of evaluation points

    p =  GeometryTypes.Point3{Float64}.(coord_source)
    tp = GeometryTypes.Point3{Float64}.(coord_eval)  

    # # rotation of points in new coordinate system where Darcy velocity is parallel to x axis
    p_rot  = rotation_z(p,-borefield.θ) 
    tp_rot = rotation_z(tp,-borefield.θ) 

    d = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp_rot) 
    d = [d[1] == 0. && d[2] == 0. ?  (0.,params.rb, d[3],d[4]) : d for d in d]

    t = tstep:tstep:tmax
    Nt = length(t) # number of time steps


    # 6. THERMAL RESPONSES
    # mutual response function between pairs of segments (adiabatic surface boundary condition)
    g = [1/(2π*params.λs)*mfls_adiabatic_surface(tt, α, coord[1:3]... ,vt, borehole.h, coord[4]; atol =1e-9) for coord in d, tt in t]
    # Matrix containing response function for each pair of segments at time-step 1
    G = g[:, :, 1]


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
                borehole.mf, borehole.cpf, bh_map, borehole.h,             
                internal_to_external
                )

    build_matrix!(M_extraction,Nb,Ns,
                k_in, k_out, k_b, 
                G,
                borehole.mf, borehole.cpf, bh_map, borehole.h,             
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
                    Δqbcurrentsum, borehole.h                  
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
    Energy_exchanged   = borehole.mf * borehole.cpf .* nbranches .* (Tfin_constraint .- Tfo) * (8760/12.)/1e6
    Energy_exchanged_by_season = collect(reshape(Energy_exchanged, 6, div(Nt,6))')
    Energy_exchanged_by_season2 = sum(Energy_exchanged_by_season , dims =2)

    E_injected   =  Energy_exchanged_by_season2[isodd.(1:length(Energy_exchanged_by_season2))]
    E_extracted  =  Energy_exchanged_by_season2[iseven.(1:length(Energy_exchanged_by_season2))]

    eta = - E_extracted ./ E_injected

    Exergy_exchanged = borehole.mf * borehole.cpf .* nbranches .* (Tfin_constraint .- Tfo  .- (T0 +273.15) .* log.( (Tfin_constraint .+ 273.15) ./ (Tfo .+ 273.15) )  ) * (8760/12.)/1e6
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
end
