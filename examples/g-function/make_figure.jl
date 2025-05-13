
function scientific_notation(x)
    exponent = floor(Int, log10(x))
    mantissa = round(x / 10. ^exponent, digits=2)
    return mantissa, exponent
end

function make_figure(configurations, spacings, compute, title; ts, H, D, rb)
    fig = Figure()

    title = Label(fig[0, 1], title, fontsize = 26)
    grid = fig[1, 1] = GridLayout()
    axes_gfunc = []
    axes_error = []

    for (k, (n, m)) in enumerate(configurations)        
        axis_gfunc = Axis(grid[1, k], ylabel = k == 1 ? L" g_{BNS}" : "", title = "$(n)x$(m) grid")
        axis_error = Axis(grid[2, k], ylabel = k == 1 ? L"\log_{10} \mid g_{pyg} - g_{BNS}\mid " : "", xlabel = L"\text{ln} \, \frac{t}{t_s}")    
        for B in spacings
            positions = [(B*(i-1), B*(j-1)) for i in 1:n for j in 1:m]

            early = compute(30, 3600*24., positions)
            mid = compute(12*10, 3600*24*30., positions)
            late = compute(100, 3600*8760*10., positions)
            far_late = compute(20, 3600*8760*1000., positions)

            tts = vcat(early[1], mid[1], late[1], far_late[1])
            Tbm = vcat(early[2], mid[2], late[2], far_late[2])
            error_gfunc = vcat(early[3], mid[3], late[3], far_late[3])

            zero_error = findall(x -> x==0., error_gfunc)
            error_gfunc[zero_error] .= eps()

            legend_value = B > 1e8 ? "\\infty" : B/H
            lines!(axis_gfunc, tts, Tbm, label = L"%$(legend_value)")
            lines!(axis_error, tts, log10.(error_gfunc))

            push!(axes_gfunc, axis_gfunc)
            push!(axes_error, axis_error)
        end
    end

    for axis_error in axes_error
        ylims!(axis_error, -18, 0)
    end

    for (i, (axis_gfunc, axis_error)) in enumerate(zip(axes_gfunc, axes_error))
        linkxaxes!(axis_gfunc, axis_error)
        if i < length(axes_gfunc)
            linkyaxes!(axes_gfunc[i], axes_gfunc[i+1])
            linkyaxes!(axes_error[i], axes_error[i+1])
        end
        hidexdecorations!(axis_gfunc, grid = false)
    end
    hideydecorations!(axes_gfunc[end], grid = false)
    hideydecorations!(axes_error[end], grid = false)
    rowgap!(grid, 5)
    colgap!(grid, 5)

    rbH_mantissa, rbH_exponent = scientific_notation(rb/H)
    DH_mantissa, DH_exponent = scientific_notation(D/H)
    ts_mantissa, ts_exponent = scientific_notation(ts)

    legend = fig[1, 2] = GridLayout()

    Label(legend[1, 1], L"t_s = %$(ts_mantissa) \times 10^{%$(ts_exponent)} \ \text{s}")
    Label(legend[2, 1], L"\frac{r_b}{H} = %$(rbH_mantissa) \times 10^{%$(rbH_exponent)}")
    Label(legend[3, 1], L"\frac{D}{H} = %$(DH_mantissa) \times 10^{%$(DH_exponent)}")
    Legend(legend[4, 1], axes_gfunc[1], L"B/H")

    fig
end