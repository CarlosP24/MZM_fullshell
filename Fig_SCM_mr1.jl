#region ## Functions 
include("plot_functions_solid.jl")

## Figure
rast = 5

function Fig1st(files, zscales; )
    fig = Figure(resolution = (0.9 * 600, 5/4 * 0.9 * 800), )
    row1 = fig[1, 1] = GridLayout()
    row2 = fig[2, 1] = GridLayout()
    row3 = fig[3, 1] = GridLayout()

    # Row 1 
    gac = row1[1, 1]  = GridLayout()
    gbl = row1[1, 2] = GridLayout()

    ga = gac[1, 1]
    gb = gbl[1, 1]
    gc = gac[2, 1]

    # Panel a - Schematics 

    ax = CairoMakie.Axis(ga; aspect = DataAspect())
    schematics(ax)

    
    # Panel b - Bands
    xlabel = L"k_z"
    ylabel = L"$\epsilon$ (meV)"
    xticks = [0]
    yticks = ([ 0, ], [ L"\mu",]) 
    limits = ((0, 0.7), (-15, 5))

    axB = CairoMakie.Axis(gb; limits, xticks, yticks, xlabel, ylabel, xlabelpadding = -18,  ylabelpadding = -20, yaxisposition = :right)
    hlines!(axB, 0; linestyle = :dash, color = :black)
    hideydecorations!(axB, grid = false, ticklabels = false)

    ax2 = CairoMakie.Axis(gb; limits,  yticks = ([0.5], [L"\epsilon"]), tellheight = false, tellwidth = false, yticklabelrotation = π/2 , yticklabelpad = 10)
    hideydecorations!(ax2, ticklabels = false)
    hidexdecorations!(ax2)
    hidespines!(ax2)

    dict = load(files[:b])["bands"]
    plot_b(dict, axB)

    # Panel c - Dome profile 
    dict = load(files[:c])
    xticks = ([0, 70], ["0", L"R"])
    yticks = ([last(dict["ϵs"]), first(dict["ϵs"])], [L"U_{min}", L"U_{max}"])
    ylabel = L"\epsilon"
    ax = CairoMakie.Axis(gc; xticks, yticks, ylabel)
    text!(ax, 10, first(dict["ϵs"]) + 3; text = L"\mu", color = :black, )
    #hidexdecorations!(ax;  ticks = false,)
    
    plot_dome(dict, ax)

    Legend(gb[1, 2], axB, L"m_J", nbanks = 1, halign =:right, valign = :bottom, orientation = :vertical; framevisible = false, titleposition = :left, titlegap = 5, rowgap = 0.1, labelsize = 10, )


    rowsize!(gbl, 1, Auto(1))
    colgap!(row1, 1, -10)
    rowgap!(gac, 1, 5)

    # Row 2
    # Panel d - LDOS 
    xlabel = L"\Phi/\Phi_0" 
    ylabel = L"$\omega$ (meV)" 
    label =  L"$$LDOS (arb. units)"
    xticks = [0, 1, 2]
    yticks = [-0.2, 0, 0.2]

    ax = CairoMakie.Axis(row2[1, 1]; xlabel, ylabel, xticks, yticks)

    dict = load(files[:d])
    zscale = zscales[:d]
    plot_LDOS(dict, ax, zscale; zlims = (0.01, 1), lowclip = :black,)
    vlines!(ax, 1.49, ymin = 0, ymax = 0.2, color = :white, linestyle = :dash)

    text!(ax, [1.5], [0.2]; text = "●", color = :white, align = (:center, :center))


    # Panel e - LDOS radial 
    xlabel = L"$r$ (nm)"
    xticks = [0, 60]

    ax = CairoMakie.Axis(row2[1, 2]; xlabel, ylabel, xticks, yticks)
    hideydecorations!(ax, ticks = false)
    
    dict = load(files[:e])
    zscale = zscales[:e]
    hmap = plot_radial(dict, ax, zscale; zlims = (0.01, 1), lowclip = :black,)

    Colorbar(row2[1, 3], limits = (0, 1), colormap = :thermal, label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

    colsize!(row2, 1, Auto(6))
    colsize!(row2, 2, Auto(1))

    colgap!(row2, 1, 10)
    colgap!(row2, 2, 5)

    # Row 3 
    # Panel f - WF 
    xlabel = L"$r$ (nm)"
    ylabel = L"$|\Psi|$ (arb. units)"
    xticks = ([0, 70], ["0", L"R"])
    yticks = ([0, 0.035], ["0", "1"])
    limits = ((0, 80), (0, 0.035))

    axW = CairoMakie.Axis(row3[1,1]; limits, xticks, yticks, xlabel, ylabel, ylabelpadding = 25)
    hidexdecorations!(axW, ticks = false, grid = false, ticklabels = false,)

    dict = load(files[:f])["WF"]
    plot_wf_v2(dict, axW)

    # Panel g - LDOS under the carpet
    xlabel = L"\Phi/\Phi_0" 
    ylabel = L"$\omega$ (meV)" 
    xticks = [0.6, 1, 1.4]
    yticks = [-0.1, 0, 0.1]

    ax = CairoMakie.Axis(row3[1, 2]; xlabel, ylabel, xticks, yticks)

    dict = load(files[:g])
    zscale = zscales[:g]
    hmap = plot_LDOS(dict, ax, zscale)

    row3[1, 3] = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
    labelpadding = -10,  width = 10, ticksize = 2,
    ticklabelpad = 5,) 

    Label(row3[1, 2, Top()], L"$m_J = 0$", padding = (-100, 0, -30, 0), color = :white)
    
    colgap!(row3, 1, 5)

    #Adjustments
    rowsize!(fig.layout, 1, Auto(1))
    rowsize!(fig.layout, 2, Auto(1))
    rowsize!(fig.layout, 3, Auto(0.5))
    #

    rowgap!(fig.layout, 1, 10)
    rowgap!(fig.layout, 2, 10)

    colgap!(row3, 2, 5)


    Label(fig.layout[1, 1, TopLeft()], "a",  padding = (-40, 0, -22, 0), font = "CMU Serif Bold")
    Label(fig.layout[1, 1, Top()], "c",  padding = (0, 0, -22, 0), font = "CMU Serif Bold")
    Label(fig.layout[1, 1, Left()], "b",  padding = (-40, 0, -22, 0), font = "CMU Serif Bold")
   
    Label(fig.layout[1, 1, Top()], L"$R = 70$ nm, $d = 10$ nm",  padding = (-230, 0, 0, 0), font = "CMU Serif Bold")


    Label(fig.layout[2, 1, TopLeft()], "d",  padding = (-40, 0, -22, 0), font = "CMU Serif Bold")
    Label(fig.layout[2, 1, Top()], "e",  padding = (265, 0, -22, 0), font = "CMU Serif Bold")

    Label(fig.layout[3, 1, TopLeft()], "f",  padding = (-40, 0, -22, 0), font = "CMU Serif Bold")
    Label(fig.layout[3, 1, Top()], "g",  padding = (-40, 0, -22, 0), font = "CMU Serif Bold")
    
    return fig
end

## Generate Figure
p = "Data/Fig_SCM_mr1/Panel_"
files = Dict( :b => "$(p)bands.jld2",  :c => "$(p)dome.jld2", :d => "$(p)LDOS.jld2", :e => "$(p)LDOS_radial.jld2", :f => "$(p)WF_v2.jld2", :g => "$(p)minigap.jld2")
zscales = Dict(:d => 8e-2, :e => 8e-2, :g => 2e-4)
f = Fig1st(files, zscales)
save("Figures/FigSCM_pos_mr1.pdf", f)
f


