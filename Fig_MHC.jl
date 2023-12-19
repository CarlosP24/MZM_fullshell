
using JLD2, CSV, DataFrames, CairoMakie, ColorSchemes, Colors, Makie.GeometryBasics

## Functions
my_dcolors = [:purple, :blue,  :red]


function plot_LDOS(dict, ax, zscale)
    local hmap 
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    m = zscale * maximum(data)
    zlims = (0.037, 1)
    data ./= m
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, lowclip = :black, rasterize = rast)
    ax.backgroundcolor = :black
    return hmap
end

function plot_mj0(dict, ax, zscale)
    local hmap 
    data = dict["mjdict"][0]
    xs, ys = dict["xs"], real.(dict["ys"])
    m = zscale * maximum(data)
    zlims = (0, 1)
    data ./= m
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = rast)
    ax.backgroundcolor = :black
    return hmap
end

function plot_dIdV(dict, ax, zmax; colormap = :thermal)
    local hmap
    data = sum(values(dict["mjdict"]))
    data = clamp.(data, 0, zmax)
    xs, ys = dict["xs"], real.(dict["ys"])
    zlims = (0, zmax)
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap, colorrange = zlims, rasterize = rast, lowclip = :black)
    ax.backgroundcolor = :black
    return hmap
end

function plot_PD(data, ax, points, μrange, αrange)
    colormap = :imola
    replace!(x -> x < 0.5 || x > 1 ? NaN : x, data)
    xlength, ylength = size(data)
    xs = range(μrange..., xlength)
    ys = range(αrange..., ylength) 
    hmap = heatmap!(ax, xs, ys, data; colormap, rasterize = rast)
    for (i, point) in enumerate(points)
        scatter!(ax, point...;)
    end
    return hmap
end

function n_radial(data, ax)
    rs = data["rs"]
    ns = data["n"]

    Mn = maximum(filter(!isnan, ns))
    mn = minimum(filter(!isnan, ns))

    xs = -80:0.1:80
    ys = -80:0.1:80 

    function tocart(x, y)
        r = sqrt(x^2 + y^2)
        ind = findmin(x -> abs(x-r), rs)[2]
        return r > 70  ? NaN : r < 50 ? 0 : ns[ind] 
    end

    n_cartesian = [tocart(x, y) for x in xs, y in ys]

    hmap = heatmap!(ax, xs, ys, n_cartesian; colormap = :thermal, rasterize = rast, colorrange = (mn*1.1, Mn), lowclip = :black )

    return hmap
end

cmap = Makie.ColorSchemes.colorschemes[:Set1_9]
ncolors = length(cmap)

function plot_wf_v2(dict, ax)
    WF_array =  vcat(collect(values(dict))...)
    ravs= map( x -> sum(x.psi .* x.rs), WF_array)
    bins = minimum(ravs):2e-3:maximum(ravs)
    c_Dict = Dict([bin => findall( x -> abs(x - bin) < 2e-3, ravs) for bin in bins])
    filter!(p -> isempty(p.second) ? false : true, c_Dict )
    color_Dict = Dict([ϵ => cmap[mod1.(unique(floor.(Int, map(x -> x.mJ, WF_array[c]))).+1, ncolors)] for (ϵ, c) in c_Dict])
    mJs = [map(x -> x.mJ, WF_array[c]) for (ϵ, c) in c_Dict][1]
    p = sortperm(mJs)
    for (ϵ, c) in c_Dict
        rs = WF_array[c[1]].rs 
        psi = WF_array[c[1]].psi
        psi = psi .- psi[1]
        nrs = length(rs)
        colD = color_Dict[ϵ]
        colD = colD[p]
        mJs = mJs[p]
        lcol_abs = length(colD)
        colors = repeat(colD, 10)
        ncol = length(colors)
        l = ceil(Int, nrs/ncol) 
        n1 = 2-l
        for (i, col) in enumerate(colors)
            n2 = minimum([n1 + l, nrs])
            label = i <= lcol_abs ? string(round(Int,mJs[i])) : nothing
            lines!(ax, rs[n1:n2], psi[n1:n2], color = col, linewidth = 3, label = label)
            n1 = n1 + l  
        end
    end
end



## Figure
rast = 5
function FigMHC(files, zscales; points = ([0.65, 50],))
    fig = Figure(resolution = (0.9 * 600, 4/3 * 650), font = "CMU Serif Roman")

    grid_up = fig[1, 1] = GridLayout() 
    grid_down = fig[2, 1] = GridLayout()

    # Row 1 
    g_radial = grid_up[1, 1] = GridLayout()
    ax = CairoMakie.Axis(g_radial[1, 1], aspect = 1)
    
    poly!(ax, Circle(Point2f(0,0), 85f0), color = colorant"#B9C8F7", strokewidth = 2)

    poly!(ax, Circle(Point2f(0,0), 77f0), color = colorant"#B9C8F7", strokewidth = 2, linestyle = :dash)


    data = load(files[:a])
    hmap = n_radial(data, ax)
    arrows!(ax, [0], [0], [55], [0], color = :white, linewidth = 3)
    arrows!(ax, [0], [0], [70 * cos(pi + pi/4)], [70 * sin(pi + pi/4)], color = :white, linewidth = 3)

    text!(ax, 10, 5; text = L"R_{av}", color = :white)
    text!(ax, -45, -15; text = L"R_{LP}", color = :white)
    hidedecorations!(ax)
    hidespines!(ax) 

    Colorbar(g_radial[2, 1], colormap = :thermal, vertical = false, flipaxis = false, height = 10, ticksize = 2, ticks = [0, 1], label = "nₑ (arb. units)", labelpadding = -20,  )

    rowgap!(g_radial, 1, 1)

    #Wave functions
    file = files[:b]
    xlabel = L"$r$ (nm)"
    ylabel = L"$|\Psi|$ (arb. units)"
    xticks = ([40, 50, 59.5, 70, 75], ["0", "", L"R_{av}", "", L"R_{LP}"])
    yticks = ([0, 0.03], ["0", "1"])
    limits = ((40, 80), (0, 0.03))
    ax = CairoMakie.Axis(grid_up[1, 2]; xlabel, ylabel, xticks, yticks,  limits, ylabelpadding = -10)
    hidexdecorations!(ax, ticks = false, ticklabels = false)

    dict = load(file)["WF"] 

    band!(ax, [50, 70], 0, 1, color = colorant"#FBF5CD",)
    band!(ax, [70, 80], 0, 1, color = colorant"#B9C8F7",)
    
    plot_wf_v2(dict, ax)

    vlines!(ax, 59.5, linestyle = :dash, color = :black, linewidth = 1)
    vlines!(ax, 75, linestyle = :dash, color = :black, linewidth = 1)

    #text!(ax, 45, 0.01; text = "Empty", align = (:center, :center))
    text!(ax, 55, 0.028; text = "SM", align = (:center, :center))
    text!(ax, 72.2, 0.028; text = "SC", align = (:center, :center))

    Legend(grid_up[1, 2], ax, L"m_J", nbanks = 1, halign =:left, valign = :center, orientation = :vertical; framevisible = false, titleposition = :top, titlegap = 5, rowgap = 0.1, labelsize = 10, tellheight = false, tellwidth = false )

    colsize!(grid_up, 1, Auto(1))
    colsize!(grid_up, 2, Auto(1)) 

    # Row 2
    #Total LDOS
    xlabel = L"\Phi/\Phi_0" 
    ylabel = L"$\omega$ (meV)" 
    label =  L"$$LDOS (arb. units)"
    xticks = [0, 1, 2]
    yticks = [-0.2, 0, 0.2]

    ax = CairoMakie.Axis(grid_down[1, 1]; xlabel, ylabel, xticks, yticks)
    dict = load(files[:c])
    zscale = zscales[:c]
    hmap = plot_LDOS(dict, ax, zscale)

    text!(ax, 0.4, -0.22; text = L"$d = 0$ nm, $g = 0$", color = :white, align = (:center, :center), font = "CMU Serif Roman",)

    hidexdecorations!(ax, ticks = false)

   
  
    text!(ax, 0.62, -0.08; text = L"L_{\Phi}", aling = (:center, :center), color = :white, fontsize = 18)
    arrows!(ax, [0.81], [-0.02], [-0.3], [0], color = :white, arrowsize = 5)
    arrows!(ax, [0.5], [-0.02], [0.3], [0], color = :white, arrowsize = 5)

    style = (color = :white,  align = (:center, :center))
    text!(ax, 1, -0.04; text = "0", style...)
    text!(ax, 1.1, -0.09; text = "-2", style...)
    text!(ax, 0.75, -0.11; text = "-1", style...)
    text!(ax, 0.86, -0.16; text = "+1, 0", style...)
    text!(ax, 0.9, -0.21; text = "+2", style...)

    text!(ax, 1, 0.04; text = "0", style...)
    text!(ax, 1.1, 0.09; text = "+2", style...)
    text!(ax, 0.75, 0.11; text = "+1", style...)
    text!(ax, 0.86, 0.16; text = "-1, 0", style...)
    text!(ax, 0.9, 0.21; text = "-2", style...)
    
    text!(ax, 1.8, 0.135; text = "-¹⁄₂", style...)
    text!(ax, 1.6, 0.05; text = "+³⁄₂", style...)
    text!(ax, 2.4, 0.05; text = "+⁵⁄₂", style...)
    text!(ax, 1.95, 0.21; text = "-³⁄₂", style...)

    text!(ax, 1.8, -0.135; text = "+¹⁄₂", style...)
    text!(ax, 1.6, -0.05; text = "-³⁄₂", style...)
    text!(ax, 2.4, -0.05; text = "-⁵⁄₂", style...)
    text!(ax, 1.95, -0.21; text = "+³⁄₂", style...)
    
    Colorbar( grid_down[1, 2] , colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

    # Row 4
    #LDOS d = 10nm, g = 10

    ax = CairoMakie.Axis(grid_down[2, 1]; xlabel, ylabel, xticks, yticks)
    dict = load(files[:d])
    zscale = zscales[:d]
    hmap = plot_LDOS(dict, ax, zscale)

    hidexdecorations!(ax, ticks = false)

    arrows!(ax, [0.88], [0], [0], [0.048], color = :white, linewidth = 2, 
    arrowsize = 12, arrowhead = :hline)
    arrows!(ax, [0.88], [0.048], [0], [-0.048], color = :white, linewidth = 2, 
    arrowsize = 12, arrowhead = :hline)

    text!(ax, 0.4, -0.22; text = L"$d = 10$ nm, $g = 10$", color = :white, align = (:center, :center), font = "CMU Serif Roman", )
    arrows!(ax, [0.5],[0.15], [0.25], [-0.11], color = :white, linewidth = 1, arrowsize = 12)

    text!(ax, 0.4, 0.2; text = "Topological \nminigap", color = :white, align = (:center, :center), font = "CMU Serif Roman", )
    Colorbar( grid_down[2, 2] , colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)


    # Row 5
    #dI/dV
    xticks = [0, 1, 2]
    xlabel = L"\Phi/\Phi_0" 
    ylabel = L"$V$ (mV)"
    label =  L"$dI/dV$ (G$_0$)"
    ax = CairoMakie.Axis(grid_down[3, 1]; xlabel, ylabel, xticks, yticks)
    dict = load(files[:e])
    zmax = zscales[:e]
    hmap = plot_dIdV(dict, ax, zmax; colormap = :thermal)

    text!(ax, 0.4, -0.22; text = L"$d = 10$ nm, $g = 10$", color = :white, align = (:center, :center), font = "CMU Serif Roman",)

    Colorbar( grid_down[3, 2] , colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)


    colgap!(grid_down, 1, 5)

    rowgap!(fig.layout, 1, 5)
    rowgap!(grid_down, 1, 5)
    rowgap!(grid_down, 2, 5) 

    rowsize!(fig.layout, 1, Auto(0.8))
    rowsize!(fig.layout, 2, Auto(3))

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -22, 0); style...)
    Label(grid_up[1, 2, TopLeft()], "b",  padding = (-30, 0, -22, 0); style...)
    Label(grid_down[1, 1, TopLeft()], "c",  padding = (-40, 0, -22, 0); style...)
    Label(grid_down[2, 1, TopLeft()], "d",  padding = (-40, 0, -22, 0); style...)
    Label(grid_down[3, 1, TopLeft()], "e",  padding = (-40, 0, -22, 0); style...)

    return fig
end

## Generate Figure

p = "Data/Fig_MHC/Panel_"
files = Dict(:a => "$(p)nrad.jld2", :b => "$(p)WF.jld2", :c => "$(p)LDOS.jld2",  :d => "$(p)LDOS_g.jld2", :e => "$(p)dIdV.jld2")
zscales = Dict(:c => 1.5e-1,  :d => 1.1e-1, :e => 1)
f = FigMHC(files, zscales)
save("Figures/FigMHC.pdf", f)
f


