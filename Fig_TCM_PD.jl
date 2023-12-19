using CairoMakie, JLD2, Images, ImageMagick, CSV, DataFrames

## Functions
rast = 5

function plot_PD(ax, dict; log = false, colorrange = (2, 4), dlim = 6e-4, colormap = :haline, μlims = nothing,   kw...)

    xs, ys = dict["xs"], dict["ys"] 
    data = dict["PD"] 

    if μlims !== nothing
        x1 = findmin(abs.(xs .- μlims[1]))[2]
        x2 = findmin(abs.(xs .- μlims[2]))[2]
        xs = xs[x1:x2]
        data = data[x1:x2, :]
    end

    log && map!(log10, ys, ys)
    data_masked = getindex.(data, 1) .* (getindex.(data, 2) .> dlim)
    replace!(x -> x == 0 ? NaN : x, data_masked)

    data_plot = log10.(data_masked .* 5) 

    heatmap!(ax, xs, ys, data_plot; colormap, colorrange, rasterize = rast)

    return xs[end]
end

function find_border(data, xrng, yrng, cond::Function)
    xlength, ylength = size(data)
    xs = range(first(xrng), last(xrng), xlength)
    ys = range(first(yrng), last(yrng), ylength)

    mesh = collect(Iterators.product(xs, ys))

    replace!(x -> cond(x) ? 0 : 1, data)

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end

function plot_MHC(ax, file, xrng, yrng; log = false)
    data = Matrix(CSV.read(file, DataFrame, header = false))
    border, dborder = find_border(data, xrng, yrng, x -> x > 0.5)
    xborder = getindex.(dborder, 1)
    yborder = getindex.(dborder, 2)
    if log 
        yborder = log10.(yborder)
    end
    return scatter!(ax, xborder, yborder; color = :orange, markersize = 3)
end

function plot_PD_LDOS(ax, dict; log = false)
    xs, ys = dict["xs"], dict["ys"]
    log && map!(log10, ys, ys)
    ldos = getindex.(dict["PD"], 2)
    return heatmap!(ax, xs, ys, ldos; colormap = :thermal, rasterize = rast, colorrange = (1e-7, 6e-4), highclip = :yellow)
end

function find_border_num(data, xs, ys; bounds = (0.1, 0.2), sigma = 1.4)
    mesh = collect(Iterators.product(xs, ys))

    edges = canny(data, bounds, sigma)
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end


function plot_islands(ax, dict; log = false, lims = nothing)
    xs, ys = dict["xs"], dict["ys"]
    log && map!(log10, ys, ys)
    ldos = dict["PD"]
    border, dborder = find_border_num(ldos, xs, ys,) 
    if lims !== nothing 
        dborder = filter(x -> (getindex(x, 2) > lims[1]) & (getindex(x, 2) < lims[2]), dborder)
    end
    xborder = getindex.(dborder, 1)
    yborder = getindex.(dborder, 2)
    return scatter!(ax, xborder, yborder; color = :red, markersize = 3)
end

## Figure

function FigTCMPD(files, files_MHC, xtickss, points, ylims)
    fig = Figure(resolution = (1100, 2/3 * 650), font = "CMU Serif Roman")

    row = 1
    xlabel = L"$\mu$ (meV)"
    ylabel = L"$\alpha$ (meVnm)"


    for (col, file) in enumerate(files[row])
        ax = CairoMakie.Axis(fig[row, col]; xlabel, ylabel, xticks = xtickss[row][col], xlabelpadding = -10)
        dict_PD, dict_islands = load.(file)
        plot_PD(ax, dict_PD)

        col == 3 && plot_MHC(ax, files_MHC[row], dict_PD["xs"], dict_PD["ys"])


        col == 5 && hlines!(ax, 0; color = :white, linewidth = 3)
        hlines!(ax, 0; color = (:gray, ifelse(col == 5, 0.5, 0.3)))

        w = ws[col]
        R = w == 0 ? 68 : 70 - w
        scatter!(ax, points[R].μ, points[R].α; color = :white, markersize = 8)
        
        lims = col == 4 ?  (-200, 130) : nothing
        plot_islands(ax, dict_islands; lims)

        col != 1 && hideydecorations!(ax, grid = false, ticks = false)

        label = col == 1 ? L"$w \rightarrow 0$" : L"$w = %$(w)$ nm"

        Label(fig[row, col, Top()], label; padding = (0, 0, 10, 0), fontsize = 20)
    end

    ylabel = L"$\Gamma_S / \Delta_0$ "
    yticks = ([0, 1, 2], [L"$1$", L"$10$", L"$10^2$"])

    row = 2

    for (col, file) in enumerate(files[row])
        ax = CairoMakie.Axis(fig[row, col]; xlabel, ylabel, yticks, xticks = xtickss[row][col], xlabelpadding = -10)
        dict_PD, dict_islands = load.(file)
        
        w = ws[col]
        R = w == 0 ? 68 : 70 - w
        plot_PD(ax, dict_PD; log = true, μlims = ylims[w])

        if col == 3
            dict_PD = load(file[1])
            ax2 = CairoMakie.Axis(fig[row, col])
            plot_MHC(ax2, files_MHC[row], dict_PD["xs"] .- 0.05, dict_PD["ys"]; log = true)
            xlims!(ax2, (17.77, 18.82))
            ylims!(ax2, (0, 3.6))
            hidedecorations!(ax2)
            hidespines!(ax2)
        end
  
        scatter!(ax, points[R].μ, points[R].τΓ |> log10; color = :white, markersize = 8)

        plot_islands(ax, dict_islands, log = true)
        
        xlims!(ax, (dict_PD["xs"][1], dict_PD["xs"][end]))
        ylims!(ax, (0, 2.1))
        col != 1 && hideydecorations!(ax, grid = false, ticks = false)
    end

    Colorbar(fig[1, 6]; colormap = :haline, ticks = ([0, 1], [L"$10^2$", L"$10^4$"]), ticklabelsvisible = true, ticksize = 2, ticklabelpad = 5, label = L"$\xi_M$ (nm)", labelpadding = -20, width = 10)
    Colorbar(fig[2, 6]; colormap = :haline, ticks = ([0, 1], [L"$10^2$", L"$10^4$"]), ticklabelsvisible = true, ticksize = 2, ticklabelpad = 5, label = L"$\xi_M$ (nm)", labelpadding = -20, width = 10)

    [colgap!(fig.layout, col, 15) for col in 1:4]
    colgap!(fig.layout, 5, 5)
    rowgap!(fig.layout, 1, 5)

    style = (fontsize = 20, font = "CMU Serif Bold")
    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -20, 0), style...)
    Label(fig[2, 1, TopLeft()], "f"; padding = (-40, 0, -20, 0), style...)

    Label(fig[1, 2, TopLeft()], "b"; padding = (-20, 0, -20, 0), style...)
    Label(fig[2, 2, TopLeft()], "g"; padding = (-20, 0, -20, 0), style...)

    Label(fig[1, 3, TopLeft()], "c"; padding = (-20, 0, -20, 0), style...)
    Label(fig[2, 3, TopLeft()], "h"; padding = (-20, 0, -20, 0), style...)

    Label(fig[1, 4, TopLeft()], "d"; padding = (-20, 0, -20, 0), style...)
    Label(fig[2, 4, TopLeft()], "i"; padding = (-20, 0, -20, 0), style...)

    Label(fig[1, 5, TopLeft()], "e"; padding = (-20, 0, -20, 0), style...)
    Label(fig[2, 5, TopLeft()], "j"; padding = (-20, 0, -20, 0), style...)


    return fig
end

## Generate Figure

p = "Data/Fig_TCM_PD/"
ws = [0, 10, 20, 30, 40]
files = [[[p*"TCM_w=$(w)_μα.jld2", p*"TCM_w=$(w)_μα_islands.jld2"]  for w in ws], [[p*"TCM_w=$(w)_μΓ.jld2", p*"TCM_w=$(w)_μΓ_islands.jld2"]  for w in ws]]

xtickss = [
    [[131, 134], [37, 41], [16, 20], [8, 12], [5, 10]],
    [[132.9, 133.4], [38.9, 39.7], [17.9, 18.7], [10.1, 10.9], [6.3, 7.3] ]
]

files_MHC = [
    "$(p)PhaseDiagram_MHC.csv",
    "$(p)PhaseDiagram_MHC_τΓ.csv"
]


points = load("$(p)TCM_params.jld2")["α=50"]

ylims = load("$(p)ylims.jld2")["μlims"]

fig = FigTCMPD(files, files_MHC, xtickss, points, ylims)
save("Figures/FigTCM_PD.pdf", fig)
fig


