using CairoMakie, JLD2, Images, ImageMagick

## Functions
function plot_LDOS(ax, dict; colormap = :thermal, zlims = (0.02, 1), kw...)
    data = sum(values(dict["mjdict"]))
    #data = dict["mjdict"][3]
    xs, ys = dict["xs"], real.(dict["ys"])
    hmap = heatmap!(ax, xs, ys, data; colormap = colormap, colorrange = zlims, rasterize = 5, lowclip = :black, kw...)
    return hmap
end

function add_axis_inset(pos=fig[1, 1]; halign, valign, width=Relative(0.5),height=Relative(0.35), alignmode=Mixed(left=5, right=5), aspect = 1)
    inset_box = CairoMakie.Axis(pos; width, height, halign, valign, alignmode, aspect)
    # bring content upfront
    translate!(inset_box.scene, 0, 0, 10)
    return inset_box
end

function R(φ, harmonics;  R0 = 1)
    return (1 + sum([har * cos(n * φ) for (n, har) in harmonics]))*R0
end

function to_cart(r, φ)
    return (r * cos(φ), r * sin(φ))
end

function find_border_num(data, xs, ys)
    mesh = collect(Iterators.product(xs, ys))

    edges = canny(data, (0.1, 0.2))
    border = findall( x -> x == 1, edges)

    return border, mesh[border]
end

## Figure

rast = 5

function FigHexagon(files, zlims_dict,; )

    fig = Figure(resolution = (0.9 * 600,  0.9 * 4/3 * 650), font = "CMU Serif Roman")

    # Panel a - Phase phase_diagram
    xlabel = L"$\tilde{\mu}$ (meV)"
    ylabel = L"$\alpha$ (meVnm)"

    ax = CairoMakie.Axis(fig[1:2, 1:2]; xlabel, ylabel, xticks = [-2, -1, 0, 1, 2], yticks = [-200, -100, 0, 100, 200])
    dict = load(files["PD"])

    μrng = dict["μrng"]
    αrng = dict["αrng"]

    heatmap!(ax, μrng, αrng, dict["PD"]; colormap = [(:white, 0), :red, :darkred], rasterize = rast)

    PD0 = dict["PD0"]
    border, μα_border = find_border_num(PD0, μrng, αrng)
    scatter!(ax, getindex.(μα_border, 1), getindex.(μα_border, 2); color = :navyblue, markersize = 2)

    scatter!(ax, 0.35, 35; color = :yellow)
    scatter!(ax, 0.5, 60; color = :lightblue)
    scatter!(ax, 1.5, 35; color = :lightgreen)
    scatter!(ax, 1, 145; color = :orange)

    Colorbar(fig[1:2, 3]; colormap = cgrad([:white, :red, :darkred], categorical = true), ticks = ([0.15, 0.5, 0.85], ["0", "1", "2"]),   width = 10, ticksize = 2, ticklabelpad = 5, label = L"$N_M$")

    xlabel = L"\Phi/\Phi_0"
    ylabel = L"$\omega$ (meV)"
    label = L"$$LDOS (arb. units)"
    xticks = [0.6, 1, 1.4]
    yticks = [-0.2, 0, 0.2]

    # Panel b - Circle 
    ax = CairoMakie.Axis(fig[3, 1]; xlabel, ylabel, xticks = xticks, yticks = yticks)
    dict = load(files["circle"])
    plot_LDOS(ax, dict; colormap = :thermal, zlims = zlims_dict["circle"], rasterize = rast)
    hidexdecorations!(ax, ticks = false,)



    text!(ax, 1.2, 0; text = L"$\simeq$", align = (:center, :center), color = :yellow, fontsize = 20)

    inset_ax = add_axis_inset(fig[3,1]; halign = 0.6, valign = :center, height = Relative(0.2), width = Relative(0.2))
    lines!(inset_ax, [to_cart(R(φ, dict_harmonics["circle"]), φ) for φ in range(0, 2π, length = 1000)], color = "yellow", linewidth = 1.5)

    inset_ax2 = add_axis_inset(fig[3,1]; halign = 0.9, valign = :center, height = Relative(0.2), width = Relative(0.2))
    lines!(inset_ax2, [to_cart(R(φ, dict_harmonics["hexagon"]), φ) for φ in range(0, 2π, length = 1000)], color = "yellow", linewidth = 1.5)

    hidespines!(inset_ax)
    hidedecorations!(inset_ax)
    hidespines!(inset_ax2)
    hidedecorations!(inset_ax2)

    # Panel c - Two MZM 
    ax = CairoMakie.Axis(fig[3, 2]; xlabel, ylabel, xticks = xticks, yticks = yticks)
    dict = load(files["twoMZM"])
    plot_LDOS(ax, dict; colormap = :thermal, zlims = zlims_dict["twoMZM"], rasterize = rast)
    hidexdecorations!(ax, ticks = false,)
    hideydecorations!(ax, ticks = false,)



    inset_ax = add_axis_inset(fig[3,2]; halign = 1, valign = :center, height = Relative(0.3), width = Relative(0.3))
    lines!(inset_ax, [to_cart(R(φ, dict_harmonics["hexagon"]), φ) for φ in range(0, 2π, length = 1000)], color = :lightblue, linewidth = 1.5)
    hidespines!(inset_ax)
    hidedecorations!(inset_ax)

    Colorbar(fig[3, 3], limits = (0, 1), colormap = :thermal, label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

    # Panel d - New MZM 
    ax = CairoMakie.Axis(fig[4, 1]; xlabel, ylabel, xticks = xticks, yticks = yticks)
    dict = load(files["newMZM"])
    plot_LDOS(ax, dict; colormap = :thermal, zlims = zlims_dict["newMZM"], rasterize = rast)
   
    inset_ax = add_axis_inset(fig[4,1]; halign = 1, valign = :center, height = Relative(0.3), width = Relative(0.3))
    lines!(inset_ax, [to_cart(R(φ, dict_harmonics["hexagon"]), φ) for φ in range(0, 2π, length = 1000)], color = :lightgreen, linewidth = 1.5)
    hidespines!(inset_ax)
    hidedecorations!(inset_ax)

    # Panel e - Split MZM
    ax = CairoMakie.Axis(fig[4, 2]; xlabel, ylabel, xticks = xticks, yticks = yticks)
    dict = load(files["splitMZM"])
    plot_LDOS(ax, dict; colormap = :thermal, zlims = zlims_dict["splitMZM"], rasterize = rast)
    hideydecorations!(ax, ticks = false,)

    inset_ax = add_axis_inset(fig[4,2]; halign = 1, valign = :center, height = Relative(0.3), width = Relative(0.3))
    lines!(inset_ax, [to_cart(R(φ, dict_harmonics["hexagon"]), φ) for φ in range(0, 2π, length = 1000)], color = :orange, linewidth = 1.5)
    hidespines!(inset_ax)
    hidedecorations!(inset_ax)

    Colorbar(fig[4, 3], limits = (0, 1), colormap = :thermal, label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10,  width = 10, ticksize = 2, ticklabelpad = 5)

    colgap!(fig.layout, 1, 10)
    colgap!(fig.layout, 2, 5)

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 3, 5)

    rowsize!(fig.layout, 1, Relative(1/6))
    rowsize!(fig.layout, 2, Relative(1/6))

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a";  padding = (-40, 0, -22, 0), style...)

    Label(fig[3, 1, TopLeft()], "b";  padding = (-40, 0, -22, 0), style...)
    Label(fig[3, 2, TopLeft()], "c";  padding = (-10, 0, -22, 0), style...)

    Label(fig[4, 1, TopLeft()], "d";  padding = (-40, 0, -22, 0), style...)
    Label(fig[4, 2, TopLeft()], "e";  padding = (-10, 0, -22, 0), style...)


    return fig
end

## Generate Figure

dict_harmonics = Dict(
    "circle" => Dict(1 => 0),
    "hexagon" => Dict( 6 => 0.0528376/(1.81709/2), 12 => 0.0146504/(1.81709/2), 18 => 0.00666731/(1.81709/2), 24 => 0.0037836/(1.81709/2), 30 => 0.00243168/(1.81709/2), 36 => 0.00169257/(1.81709/2), 42 => 0.00124527/(1.81709/2), 48 => 0.000954284/(1.81709/2))
    )


p = "Data/Fig_Hexagon/"
try run(`mkdir Figures`) catch end

files = Dict(
    "PD" => "$(p)PD_Hexagon.jld2", 
    "circle" => "$(p)hexagon.jld2", 
    "splitMZM" => "$(p)splitMZM.jld2", 
    "splitMZM_noMM" => "$(p)splitMZM_noMM.jld2", 
    "newMZM" => "$(p)newMZM.jld2",
    "twoMZM" => "$(p)twoMZM.jld2"
    )
zlims_dict = Dict(
    "circle" => (2e-4, 2e-2), 
    "twoMZM" => (2e-4, 3e-2),
    "newMZM" => (2e-4, 3e-2), 
    "splitMZM" => (2e-3, 1e-1), 
    )

fig = FigHexagon(files, zlims_dict,)
save("Figures/FigHexagon.pdf", fig)
fig

