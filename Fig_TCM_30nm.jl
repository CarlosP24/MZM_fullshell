using CairoMakie, Images, ImageMagick, JLD2, CSV, DataFrames

## Functions

function find_border(data, xrng, yrng, cond::Function)
    xlength, ylength = size(data)
    xs = range(xrng..., xlength)
    ys = range(yrng..., ylength)

    mesh = collect(Iterators.product(xs, ys))

    replace!(x -> cond(x) ? 0 : 1, data)

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end

function find_border_num(data, xs, ys; bounds = (0.1, 0.2), sigma = 1.4)
    mesh = collect(Iterators.product(xs, ys))

    edges = canny(data, bounds, sigma)
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end


function plot_2PD(ax, dict, data_a, nlim)
    xs, ys = dict["xs"], dict["ys"]
    data = dict["PD"]
    replace!(x -> x < nlim ? 0 : 1, data_a) 
    data_masked = getindex.(data, 1) .* data_a
    replace!(x -> x == 0 ? NaN : x, data_masked)
    data_plot = log10.(data_masked .* 5) 

    heatmap!(ax, xs, ys, data_plot; colormap = :haline,  rasterize = 5, colorrange = (2, 3.5))
    
    mesh = collect(Iterators.product(xs, ys))
    edges = canny(data_a, (0.1, 0.2))
    border = findall( x -> x == 1, edges)
    dborder = mesh[border]

    scatter!(ax, getindex.(dborder, 1), getindex.(dborder, 2), color = :orange, markersize = 4)

end

## Figure

function FigTCM_30nm(file, filePD, zrange; )
    fig = Figure(resolution = (0.9 * 600, 550), font = "CMU Serif Roman")
    
    ax = CairoMakie.Axis(fig[1,1]; xlabel = L"$\mu$ (meV)", ylabel = L"$\alpha$ (meVnm)", xticks = [36, 38, 40], yticks = [-200, -100, 0, 100, 200], xlabelpadding = -7
)



    dict = load(filePD[1])


    data_a = Matrix(CSV.read(filePD[2], DataFrame, header = false))


    plot_2PD(ax, dict, data_a, 0.66)

    dict_islands = load(filePD[3])
    μrng = dict_islands["xs"]
    αrng = dict_islands["ys"]
    PD_islands = dict_islands["PD"]
    border, μα_border = find_border_num(PD_islands, μrng, αrng, bounds = (0.1, 0.2), sigma = 1)
    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2), color = :red, markersize = 3)
    scatter!(ax, 39.2, 10; color = :white, markersize = 8)


    Colorbar(fig[1, 2], colormap = :haline, limits = (0, 1), label = L"$\xi_M$ (nm)", ticklabelsvisible = true, ticks = ([0,1], [L"10^2", L"10^4"]), labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)
    style = (font = "CMU Serif Bold", fontsize = 20)

    # Panel b - LDOS
    xlabel = L"\Phi/\Phi_0"
    ylabel = L"$\omega$ (meV)"
    label = L"$$LDOS (arb. units)"
    yticks = [-0.2, 0, 0.2]
    xlabelpadding = -7
    ax = CairoMakie.Axis(fig[2, 1]; xlabel, ylabel, yticks, xlabelpadding, )

    dict = load(file)
    mjdict = dict["mjdict"]
    heatmap!(ax, dict["xs"], real.(dict["ys"]), values(mjdict) |> sum; colormap = :thermal, colorrange = zrange, lowclip = :black, rasterize = 5, )
    
   

    Colorbar(fig[2, 2], colormap = :thermal, limits = (0, 1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)
    style = (font = "CMU Serif Bold", fontsize = 20)

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -15, 0); style...)
    Label(fig[2, 1, TopLeft()], "b", padding = (-40, 0, -15, 0); style...)

    rowgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 1, 5)
    return fig, μα_border
end

## Generate Figure 

p = "Data/Fig_TCM_30nm/"
try run(`mkdir Figures`) catch end

filePD = ["$(p)TCM_w=10nm_PD_2wedges.jld2",  "$(p)PhaseDiagram_30_μα.csv", "$(p)TCM_w=10nm_islands_2wedges.jld2"]

file = "$(p)w=10_nω.jld2"


zrange = (5e-4, 9e-3)

fig, μα_border = FigTCM_30nm(file, filePD, zrange)
save("Figures/FigTCM_30nm.pdf", fig)
fig
