using JLD2, CSV, DataFrames, CairoMakie, ColorSchemes, Images, ImageMagick
## Functions

# Plotting
function plot_Total(dict, ax, zscale)
    local hmap
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    xmin = findmin( x -> abs(x - 0), xs)[2]
    xmax = findmin( x -> abs(x - 1.5), xs)[2]
    data = data[xmin:xmax, :]
    xs = xs[xmin:xmax]
    m = zscale * maximum(data)
    zlims = (0.04, 1)
    data ./= m
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = rast, lowclip = :black)
    ax.backgroundcolor = :black
    return hmap
end

function plot_mJ0(dict, ax, zscale)
    local hmap 
    data = dict["mjdict"][0]
    xs, ys = dict["xs"], real.(dict["ys"])
    m = zscale * maximum(data)
    zlims = (0.01, 1)
    data ./= m
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = rast, lowclip = :black)
    ax.backgroundcolor = :black
    text!(ax, 1, 0; text = L"m_J = 0", color = :white, align = (:center, :center))
    return hmap
end

function plot_2PD(ax, data, dict,μrange, αrange, δx, δy; points = nothing)
    colormap = :imola
    replace!(x -> x < 0 || x > 0.5 ? NaN : x, data)
    xlength, ylength = size(data)
    xs = range(μrange..., xlength)
    ys = range(αrange..., ylength) 
    hmap = heatmap!(ax, xs, ys, data; colormap, colorrange = (0, 0.4), rasterize = rast, )

    xs, ys = dict["xs"], dict["ys"]
    data = dict["PD"]
    border, dborder = find_border(data, xs, ys, x -> x == 0)
    xborder = getindex.(dborder, 1)
    yborder = getindex.(dborder, 2)
    scatter!(ax, xborder .- δx, yborder .* δy, color = :red, markersize = 3)

    if points !== nothing
        for (i, point) in enumerate(points)
            scatter!(ax, point[1], point[2]; color = point[3],)
        end
    end

    return hmap
end

function plot_minigap(dict, ax, zscale)
    local hmap
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    xmin = findmin( x -> abs(x - (-10)), xs)[2]
    xmax = findmin( x -> abs(x - 1), xs)[2]
    #ymin = findmin( x -> abs(x - (-0.05)), ys)[2]
    #ymax = findmin( x -> abs(x - 0.23), ys)[2]
    data = data[xmin:xmax, :]
    xs = xs[xmin:xmax]
    # ys = ys[ymin:ymax]
    m = zscale * maximum(data)
    zlims = (0.04, 1)
    data ./= m
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, lowclip = :black, rasterize = rast)
    ax.backgroundcolor = :black
    return hmap
end

function plot_PD(data, ax, xrange, yrange; points = nothing)
    colormap = :imola
    #replace!(x -> x < 0 || x > 0.5 ? NaN : x, data)
    xlength, ylength = size(data)
    xs = range(xrange..., xlength)
    ys = range(yrange..., ylength) 

    for (i, R) in enumerate(xs)
        n = find_nRborder(R)
        for (j, α) in enumerate(ys)
            data[i, j] = data[i, j] - n
            data[i, j] = ifelse(data[i, j] < 0 || data[i, j] > 0.5, NaN, data[i, j])
        end
    end

    hmap = heatmap!(ax, xs, ys, data; colormap, colorrange = (0, 0.4), rasterize = rast,)

    if points !== nothing
        for (i, point) in enumerate(points)
            scatter!(ax, point[1], point[2]; color = point[3],)
        end
    end
    return hmap
end

# Utilities

function get_Z(Z)
    if Z > 0
        return L"$+%$(Z)$"
    else
        return L"$%$(Z)$"
    end
end


function find_border(data, xs, ys, cond::Function)
    mesh = collect(Iterators.product(xs, ys))

    replace!(x -> cond(x) ? 0 : 1, data)

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end

# R vs. α PD. Finding destructive LP border.
function pairbreaking(n, Δ0, ξd, Rcore, Rshell)
    RLP = (Rcore + Rshell) / 2
    dshell = Rshell - Rcore
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (4 * (n - round(n))^2 + dshell^2 / RLP^2 * (n^2 + (round(n)^2)/3))
    return Λ
end

function uUsadel(Δ0, Λ, ω)
    Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usa = 1/(2 * Δd) *(ω + sign(ω) * rai - sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
    return usa
end

n = 0.5
Δ0 = 0.23
ξd = 70

u0(n, R) = uUsadel(Δ0, pairbreaking(n, Δ0, ξd, R, R), 0 + 1e-3*im)

function find_nRborder(R; ns = range(0.501, 1, length = 1000))
    indn = findfirst(imag.(u0.(ns, R)) .< -1e1)
    indn !== nothing && return ns[indn]
    return 0.501
end

## Figure

rast = 5
function FigHollow(files, zscales, TTs, points_α, points_R, points_Γ, panel_markers; kmax = 0.5, bands_cmap = [ColorSchemes.gist_rainbow.colors[x] for x in 1:9:99], ylims_bands = (-2.2, 5))
    fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")

    # Grids 
    grid_left = fig[1, 1] = GridLayout()

    grid_a = grid_left[1, 1] = GridLayout() 
    grid_b = grid_left[1, 2] = GridLayout() 
    grid_c = grid_left[2, 1] = GridLayout() 
    grid_d = grid_left[2, 2] = GridLayout()
    grid_e = grid_left[3, 1] = GridLayout()
    grid_f = grid_left[3, 2] = GridLayout(alignmode = Mixed(left = 0, ),)

    grid_right = fig[1, 2] = GridLayout()


    # LDOS block 
    xlabel = L"\Phi/\Phi_0"
    ylabel = L"$ω$ (meV)"
    label =  L"$$LDOS (arb. units)"
    xticks0 = [0.6, 1, 1.4]
    xticksT = [0, 0.5, 1]

    # Panel a: no SOC 
    file = files["no_SOC"]*".jld2"
    zscale = zscales["noSOC"]
    dict = load(file)
    ax = CairoMakie.Axis(grid_a[1, 1]; xlabel, ylabel, xticks = xticksT)

    hmap = plot_Total(dict, ax, zscale)
    hidexdecorations!(ax, ticks = false)

    file = files["no_SOC"]*"_0.jld2"
    dict = load(file)
    zscale = zscales["noSOC_0"]

    ax = CairoMakie.Axis(grid_a[1, 2]; xlabel, xticks = xticks0)
    hmap = plot_mJ0(dict, ax, zscale)
    hidexdecorations!(ax, ticks = false)
    hideydecorations!(ax, )

    # Panel b: SOC trivial 
    file = files["SOCtrivial"]*".jld2"
    zscale = zscales["SOCtrivial"]
    dict = load(file)

    ax = CairoMakie.Axis(grid_b[1, 1]; xlabel, ylabel, xticks = xticksT)

    hmap = plot_Total(dict, ax, zscale)
    hidexdecorations!(ax, ticks = false)
    hideydecorations!(ax, ticks = false)

    text!(ax, 1, 0; text = "Degeneracy \npoints", align = (:center, :center), color = :white, fontsize = 14, )
    arrows!(ax, [1], [0.045], [0], [0.05], color = :white)
    arrows!(ax, [1], [-0.045], [0], [-0.05], color = :white)

    file = files["SOCtrivial"]*"_0.jld2"
    zscale = zscales["SOCtrivial_0"]
    dict = load(file)
    ax = CairoMakie.Axis(grid_b[1, 2]; xlabel, xticks = xticks0)

    hmap = plot_mJ0(dict, ax, zscale)
    hidexdecorations!(ax, ticks = false)
    hideydecorations!(ax)

    Colorbar(grid_left[1, 3], colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)


    # Panel c: TT 
    file = files["big_SOC_trivial"]*".jld2"
    dict = load(file)
    zscale = zscales["big_SOC_trivial"]
    ax = CairoMakie.Axis(grid_c[1, 1]; xlabel, ylabel, xticks = xticksT)

    hmap = plot_Total(dict, ax, zscale)

    file = files["big_SOC_trivial"]*"_0.jld2"
    dict = load(file)
    zscale = zscales["big_SOC_trivial_0"]

    ax = CairoMakie.Axis(grid_c[1, 2]; xlabel, xticks = xticks0,)
    hmap = plot_mJ0(dict, ax, zscale)

    hideydecorations!(ax,)

    # Panel d: Topological
    file = files["topo"]*".jld2"
    zscale = zscales["topo"]
    dict = load(file)

    ax = CairoMakie.Axis(grid_d[1, 1]; xlabel, ylabel, xticks = xticksT)

    hmap = plot_Total(dict, ax, zscale)
    text!(ax, 0.54, 0.03; text = L"L_{\Phi}", aling = (:center, :center), color = :white, fontsize = 15)
    arrows!(ax, [0.67], [0.02], [-0.15], [0], color = :white, arrowsize = 5)
    arrows!(ax, [0.5], [0.02], [0.15], [0], color = :white, arrowsize = 5)
    hideydecorations!(ax, ticks = false)
    vlines!(ax, 0.51; ymax = 0.2, color = :white, linestyle = :dash)

    
    file = files["topo"]*"_0.jld2"
    zscale = zscales["topo_0"]
    dict = load(file)

    ax = CairoMakie.Axis(grid_d[1, 2]; xlabel, xticks = xticks0)
    hmap = plot_mJ0(dict, ax, zscale)
    hideydecorations!(ax,)

    Colorbar(grid_left[2, 3],colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

    grids_LDOS = [grid_a, grid_b, grid_c, grid_d]
    for (grid, panel_marker) in zip(grids_LDOS, panel_markers)
        colsize!(grid, 1, Auto(3/2))
        colgap!(grid, 1, 5)
        grid != grid_a && Label(grid[1, 1, Top()], "●", padding = (-50, 0, -30, 0), font = "CMU Serif Bold", color = panel_marker)
        grid == grid_a && Label(grid[1, 1, Top()], L"\alpha = 0", padding = (-50, 0, -30, 0), font = "CMU Serif Bold", color = :white)
    end 

    # Panel e: bands
    file = files["bands"]
    bands_dict = load(file)["bands"]
    Zs = sort(collect(keys(bands_dict)))

    xlabel = L"k_z" 
    ylabel = L"$ε$"
    xticks = [0]
    yticks = ([ 0, ], [ L"μ",]) 
    ax = CairoMakie.Axis(grid_e[1, 1]; xlabel, xticks, yticks, yaxisposition = :right)
    
    ylims!(ax, ylims_bands)
    xlims!(ax, (0, kmax))

    cmap = bands_cmap
    ncolors = length(cmap)
    
    for Z in Zs
        color = cmap[mod1(Z+1, ncolors)]
        pts = bands_dict[Z]
        foreach(p -> lines!(p; color, label = (p == first(pts)) & (Z in -4:4) ? get_Z(Z) : nothing, linewidth = 3), pts)
    end

    hlines!(ax, 0; linestyle = :dash, color = :black)

    text!(ax, 0.3, -0.9; text = "●", color = :brown, align = (:center, :center))
    text!(ax, 0.33, -1; text = L"$\Phi = 0.51 \Phi_0$",  align = (:left, :center), fontsize = 15)

    

    grid_e[1, 2] = Legend(fig, ax, L"m_J",  framevisible = false, nbanks = 1, tellheight = false, tellwidth = false, labelsize = 14, rowgap = 1, patchsize = (15, 5) )
    
    yticks = ([-0.9], [L"\tilde{\mu}",])
    ax2 = CairoMakie.Axis(grid_e[1, 1]; ylabel, tellwidth = false, width = 260, yticks, ylabelpadding = -8)

    hideydecorations!(ax2; label = false, ticklabels = false)
    hidexdecorations!(ax2)
    hidespines!(ax2)
    ylims!(ax2, ylims_bands)
    xlims!(ax2, (-.02, kmax))
    arrows!(ax2, [0], [0, -1.8], [0, 0], [-1.8, 1.8]; arrowhead = '-', arrowsize = 40); 
    colsize!(grid_e, 1, Auto(3))
    colsize!(grid_e, 2, Auto(1))
    colgap!(grid_e, 1, 5)



    # Panel f: under the carpet

    file = files["topo_minigap"]
    zscale = zscales["topo_minigap"]
    dict = load(file)

    ax = CairoMakie.Axis(grid_f[1, 1]; xlabel = L"\Phi / \Phi_0", ylabel = L"$\omega$ (meV)")

    hmap = plot_minigap(dict, ax, zscale)

    ax.backgroundcolor = :black
    
    text!(ax, -6.9, -0.14; text = L"\Phi_{TT}^{(1)}", aling = (:center, :center), color = :white, fontsize = 15)
    arrows!(ax, [-7], [-0.1], [-0.9], [0.06], color = :white)

    text!(ax, -2, 0.05; text = L"\Phi_{TT}^{(2)}", aling = (:center, :center), color = :white, fontsize = 15)
    arrows!(ax, [-0.6], [0.08], [0.8], [-0.04], color = :white)


    Colorbar(grid_left[3, 3], colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)
    
    colgap!(grid_left, 1, 10)
    colgap!(grid_left, 2, 10)
    rowgap!(grid_left, 1, 5)

    Label(grid_f[1, 1, Top()], "●", color = last(panel_markers), padding = (-100, 0, -40, 0))

    Label(grid_f[1, 1, Top()], L"m_J = 0", color = :white, padding = (0, 0, -40, 0))

    # Panel g: α v μ 
    file = files["PD_μα"]

    xlabel = L"$\tilde{\mu}$ (meV)"
    ylabel = L"$α$ (meVnm)"

    μrange = (-1.5, 1.5)
    αrange = (-200, 200)
    xticks = [0, 0.75, 1.5]
    yticks = [0, 50, 100, 150, 200]

    ax = CairoMakie.Axis(grid_right[1, 1]; xlabel, ylabel, xticks, yticks)

    data = Matrix(CSV.read(file[1], DataFrame, header = false))
    dict = load(file[2])

    data = data .- 0.5 # Change Φ_TT to L_Φ
    δx = 132.53
    hmap = plot_2PD(ax, data, dict,  μrange, αrange, δx, 1; points = points_α)

    xlims!(ax, (0, 1.5))
    ylims!(ax, (0, 200))
    text!(ax, 0.1, 175; text = L"$R = 70$ nm", align = (:left, :center), fontsize = 14)
    text!(ax, 0.1, 155; text = L"$\Gamma_S = \Delta_0$", align = (:left, :center), fontsize = 14)


    hidexdecorations!(ax, ticks = false)
    ax.backgroundcolor = :white 

    Colorbar(grid_right[1, 2], colormap = :imola, limits = (0, 0.4), label = L"$L_{\Phi}$", ticklabelsvisible = true, ticks = [0,0.3], labelpadding = -25,  width = 10, ticksize = 2, ticklabelpad = 5)


    # Panel h: Γ v μ 

    file = files["PD_μΓ"]


    ylabel = L"$\Gamma_S / \Delta_0$"

    μrange = (0, 1.5)
    Γrange = (0, 2)
    xticks = range(μrange..., 3)
    yticks = range(Γrange..., 5)

    ax = CairoMakie.Axis(grid_right[2, 1]; xlabel, ylabel, xticks, yticks)

    data = Matrix(CSV.read(file[1], DataFrame, header = false))
    dict = load(file[2])

    data = data .- 0.5 # Change Φ_TT to L_Φ
    δx = 132.4
    δy = 1.23
    hmap = plot_2PD(ax, data, dict, μrange, Γrange, δx, δy; points = [[0.75, 1, :brown]])
    hlines!(ax, 0, color = :white, linewidth = 3)

    text!(ax, 0.1, 1.8; text = L"$R = 70$ nm", align = (:left, :center), fontsize = 14)
    text!(ax, 0.75, 1.8; text = L"$\alpha = 85$ meVnm", align = (:left, :center), fontsize = 14)




    ax.backgroundcolor = :white 

    Colorbar(grid_right[2, 2], colormap = :imola, limits = (0, 0.4), label = L"$L_{\Phi}$", ticklabelsvisible = true, ticks = [0,0.3], labelpadding = -25,  width = 10, ticksize = 2, ticklabelpad = 5)

    # Panel i: α v R 

    file = files["PD_Rα"]


    xlabel = L"$R$ (nm)"
    ylabel = L"$\alpha$ (meVnm)$"
    Rrange = (0, 140)
    αrange = (0, 200)


    xticks = range(Rrange..., 5)
    yticks = range(αrange..., 5)

    ax = CairoMakie.Axis(grid_right[3, 1]; xlabel, ylabel, xticks, yticks)

    data = Matrix(CSV.read(file[1], DataFrame, header = false))
    
    hmap = plot_PD(data, ax,  Rrange, αrange;)
    hlines!(ax, 0, color = :white, linewidth = 3)
    text!(ax, 5, 180; text = L"$\mu = 0.75$ meV", align = (:left, :center), fontsize = 14)
    text!(ax, 5, 160; text = L"$\Gamma_S = \Delta_0$", align = (:left, :center),fontsize = 14)

    vlines!(ax, 33.3;  ymax = 0.6, color = :black, linestyle = :dash, linewidth = 3)

    dict_islands = load(file[2])
    αs, Rs = dict_islands["xs"], dict_islands["ys"]
    islands = dict_islands["PD"]
    border, dborder = find_border(islands, αs, Rs, x -> x == 0)
    xborder = getindex.(dborder, 1)
    yborder = getindex.(dborder, 2)
    scatter!(ax, xborder, yborder, color = :red, markersize = 3)

    for point in points_R
        scatter!(ax, point[1], point[2]; color = point[3],)
    end

    text!(ax, 15, 60; text = "Destructive LP", align = (:center, :center), rotation = π/2, fontsize = 12)


    ax.backgroundcolor = :white 
    Colorbar(grid_right[3, 2], colormap = :imola, limits = (0, 0.4), label = L"$L_{\Phi}$", ticklabelsvisible = true, ticks = [0,0.3], labelpadding = -25,  width = 10, ticksize = 2, ticklabelpad = 5)


    colgap!(grid_right, 1, 5)
    rowgap!(grid_right, 1, 5)
    rowgap!(grid_right, 2, 5)
    rowgap!(grid_left, 2, 5)

    colsize!(fig.layout, 2, Relative(0.7/3))


    Label(grid_left[1, 1, TopLeft()], "a", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_left[1, 2, TopLeft()], "b", padding = (-10, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_left[2, 1, TopLeft()], "c", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_left[2, 2, TopLeft()], "d", padding = (-10, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_left[3, 1, TopLeft()], "e", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_left[3, 2, TopLeft()], "f", padding = (-10, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_right[1, 1, TopLeft()], "g", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_right[2, 1, TopLeft()], "h", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(grid_right[3, 1, TopLeft()], "i", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", fontsize = 20)



    return fig 
end

## Generate Figure

path = "Data/Fig_HCA/"
p = "$(path)Panel_"
try run(`mkdir Figures`) catch end


files = Dict("PD_μα" => ("$(path)PhaseDiagram_μα.csv", "$(path)TCM_w=0_μα_islands.jld2"), "PD_Rα" => ("$(path)PhaseDiagram_Rα.csv", "$(path)TCM_w=0_islands_Rα.jld2"), "PD_μΓ" => ("$(path)PhaseDiagram_μΓ.csv", "$(path)TCM_w=0_μΓ_islands.jld2"),
"bands" => "$(p)bands.jld2", "no_SOC" => "$(p)noSOC", "SOCtrivial" => "$(p)SOCtrivial", "big_SOC_trivial" => "$(p)big_SOC_trivial", "smalltopo" => "$(p)SOC_topo_small", "topo" => "$(p)SOCtopo", "topo_minigap" => "$(p)minigap.jld2")


zscales = Dict("noSOC" => 5e-1, "noSOC_0" => 2e-1, "SOCtrivial" => 4.5e-1, "SOCtrivial_0" => 4e-2, "big_SOC_trivial" => 4.5e-1, "big_SOC_trivial_0" => 1e-1, "smalltopo" => 3e-1, "smalltopo_0" => 4e-3, "topo" => 2e-1, "topo_0" => 2.5e-3,  "topo_minigap" => 1e-3)

TTs = [0.696599, -8.25365] 
points_α = [ [0.75, 45, :fuchsia], [0.75, 62, :orange], [0.75, 85,:brown]] |> reverse
points_R = [ [70, 45, :fuchsia], [70, 62, :orange], [70, 85,:brown]] |> reverse
panel_markers = [:white, :fuchsia, :orange,  :brown]

f = FigHollow(files, zscales, TTs, points_α, points_R, 0, panel_markers)
save("Figures/FigHCA.pdf", f, overwrite = true)
f
