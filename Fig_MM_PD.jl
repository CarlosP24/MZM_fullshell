using JLD2, CairoMakie, Images, ImageMagick

## Functions

rast = 5

function find_border_num(data, xs, ys)
    mesh = collect(Iterators.product(xs, ys))

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end

## Figure

function Fig_MM_PD(files; Vtop = -50, αtop = 100)
    fig = Figure(resolution = (0.9 * 600, 650), font = "CMU Serif Roman")

    # TCM - MM 
    # row 1 
    grid_up = fig[1, 1] = GridLayout()
    
    xlabel = L"$\tilde{\mu}$ (meV)"
    ylabel = L"$\alpha$ (meVnm)"

    # n = 1
    ax =  CairoMakie.Axis(grid_up[1, 1]; xlabel, ylabel)
    dict = load(files["TCM_n=1"])
    PD = dict["PD"]
    μrng = range(-2.5, 2.5, size(PD) |> first)
    αrng = range(-200, 200, size(PD) |> last)
    heatmap!(ax, μrng, αrng, PD; colormap = [(:white, 0), :red], rasterize = rast)

    PD0 = dict["PD0"]
    border, μα_border = find_border_num(PD0, μrng, αrng)
    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2); color = :navyblue, markersize = 2)

    text!(ax, -1.5, 2; text = L"$n=1$", align = (:center, :center))

    # n = 2
    ax = CairoMakie.Axis(grid_up[1, 2]; xlabel, ylabel)
    dict = load(files["TCM_n=2"])
    PD = dict["PD"]
    μrng = range(-2.5, 2.5, size(PD) |> first)
    αrng = range(-200, 200, size(PD) |> last)
    heatmap!(ax, μrng, αrng, -PD; colormap = [(:white, 0), :red], rasterize = rast)

    text!(ax, -1.5, 2; text = L"$n=2$", align = (:center, :center))


    hideydecorations!(ax, ticks = false, grid = false)

    Colorbar(fig[1, 2]; colormap = cgrad([:white, :red,], categorical = true), ticks = ([0.25, 0.75], ["0", "1"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$", labelpadding = -7)

    colgap!(grid_up, 1, 10)
    # SCM - Dome - MM 
    xlabel = L"$U_{min}$ (meV)"
    ylabel = L"$\langle \alpha \rangle$ (meVnm)"

    # n = 1
    ax = CairoMakie.Axis(fig[2, 1]; xlabel, ylabel, xreversed = true)
    dict = load(files["SCM_n=1"])
    PD = dict["PD"]
    Vrng = range(-1, -50, size(PD) |> first)
    αrng = range(-100, 100, size(PD) |> last)


    heatmap!(ax, Vrng, αrng, -PD; colormap = [(:white, 0), :red], rasterize = rast)

    PD0 = dict["PD0"]
    border, μα_border = find_border_num(PD0, Vrng, αrng)
    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2); color = :navyblue, markersize = 2)

    text!(ax, -3, 0; text = L"$n = 1$", rotation = π/2, align = (:center, :center))

    hidexdecorations!(ax, ticks = false, grid = false)

    Colorbar(fig[2, 2]; colormap = cgrad([:white, :red,], categorical = true), ticks = ([0.25, 0.75], ["0", "1"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$", labelpadding = -7)


    # n = 2 
    ax = CairoMakie.Axis(fig[3, 1]; xlabel, ylabel, xreversed = true)
    dict = load(files["SCM_n=2"])
    PD = dict["PD"]
    Vrng = range(-1, -50, size(PD) |> first)
    αrng = range(-100, 100, size(PD) |> last)


    heatmap!(ax, Vrng, αrng, -PD; colormap = [(:white, 0), :red], rasterize = rast)

    text!(ax, -3, 0; text = L"$n = 2$", rotation = π/2, align = (:center, :center))

    Colorbar(fig[3, 2]; colormap = cgrad([:white, :red,], categorical = true), ticks = ([0.25, 0.75], ["0", "1"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$", labelpadding = -7)

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 10)
    colgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)

    Label(fig[1, 1, Top()], "Tubular-core model"; padding = (0, 0, 10, 0))

    Label(fig[2, 1, Top()], "Solid-core model"; padding = (0, 0, 10, 0))


    Label(fig[1, 1, TopLeft()], "a";  padding = (-40, 0, -50, 0), style...)

    Label(fig[1, 1, Top()], "b";  padding = (0, 0, -50, 0), style...)

    Label(fig[2, 1, TopLeft()], "c";  padding = (-40, 0, -30, 0), style...)

    Label(fig[3, 1, TopLeft()], "d";  padding = (-40, 0, -30, 0), style...)
    return fig
end

## Generate Figure 

p = "Data/Fig_MM_PD/"
try run(`mkdir Figures`) catch end

files = Dict(
    "TCM_n=1" => "$(p)PD_TCM_n=1.jld2",
    "TCM_n=2" => "$(p)PD_TCM_n=2.jld2",
    "SCM_n=1" => "$(p)PD_SCM_dome_n=1.jld2",
    "SCM_n=2" => "$(p)PD_SCM_dome_n=2.jld2"
)

f = Fig_MM_PD(files)
save("Figures/FigMM_PD.pdf", f)
f