using JLD2, CSV, DataFrames, CairoMakie, ColorSchemes, Images, ImageMagick

## Functions
rast = 5
function plot_Top(data, ax, xrange, yrange)
    xlength = size(data)[1]
    ys = range(yrange..., xlength)
    lines!(ax, data, ys; linewidth = 3, color = colorant"#697a5c", rasterize = rast) 
    band!(ax, data, yrange[2], ys;  color = :white, rasterize = rast)
end

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

## Figure
function plot_PD_sup(files)
    fig = Figure(resolution = (2/3 * 0.9 * 600, 2/3 *  2/3 * 650), font = "CMU Serif Roman")

    xlabel = L"$\tilde{\mu}$ (meV)"
    ylabel = L"$\alpha$ (meVnm)"

    μrng = (-5, 10)
    αrng = (-150, 150)

    xticks = [-2, 0, 2, 4]
    yticks = ([-100, -23.6, 0, 100], ["-100",L"\alpha_c", "0", "100"])

    ax = CairoMakie.Axis(fig[1, 1]; xlabel, ylabel, xticks, yticks)

    file = files["μ_top"]
    data = vec(Matrix(CSV.read(file, DataFrame, header = false)))
    plot_Top(data, ax, μrng, αrng)
    

    file = files["PD"]
    data_a = Matrix(CSV.read(file, DataFrame, header = false))
    border, μα_border = find_border(data_a, μrng, αrng, x -> x < 0.51  )
    xborder = getindex.(μα_border, 1)
    yborder = getindex.(μα_border, 2)
    scatter!(ax, xborder, yborder, color = :blue, markersize = 4)
    ax.backgroundcolor = (colorant"#d1f5b9", 0.5)

    sec_left = μα_border[findall(x -> getindex(x, 1) < 0, μα_border)]
    sec_right = μα_border[findall(x -> getindex(x, 1) > 0, μα_border)]

    #band!(ax, getindex.(sec_left, 1),αrng[1], getindex.(sec_left, 2),; color = (:blue, 0.5), rasterize = rast)
    #band!(ax, getindex.(sec_right, 1), getindex.(sec_right, 2), αrng[2],; color = (:blue, 0.5), rasterize = rast)

    hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
    hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
    hlines!(ax, -23.6; color = :black, linestyle = :dash)

    text!(ax, -2, 50; text = "Trivial", align = (:center, :center))
    text!(ax, 2, 50; text = "Topological", align = (:center, :center))
    xlims!(ax, (-3, 3))
    ylims!(ax, αrng)
    
    return fig 
end

## Generate Figure

p = "Data/Sup_Topology/"
try run(`mkdir Figures`) catch end

files = Dict(
    "μ_top" => "$(p)μ_top.csv",
    "PD" => "$(p)PhaseDiagram_μα_sup.csv"
)

fig = plot_PD_sup(files)
save("Figures/Sup_Topology.pdf", fig)
fig