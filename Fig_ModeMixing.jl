using Random, Distributions, CairoMakie, JLD2

## Functions 

prefactors = Dict( 
    1 => sqrt(6) / π,
    2 => 3*sqrt(10)/π^2)

function harmonics_array(σ, i, ℓmax; prefactors = prefactors)
    Random.seed!(123321)
    σ1 = prefactors[i] * σ
    d(ℓ) = Distributions.Normal(0, σ1/ℓ^i) 
    return [rand(d(ℓ)) * exp(2π * rand() * im) for ℓ in 1:ℓmax]
end
 
function plot_LDOS(ax, dict, zscale; colormap = :thermal, zlims = (0, 1), kw...)
    data = sum(values(dict["mjdict"]))
    m = zscale * maximum(data)
    data ./= m
    xs, ys = dict["xs"], real.(dict["ys"])
    hmap = heatmap!(ax, xs, ys, data; colormap = colormap, colorrange = zlims, rasterize = 5, kw...)
    return hmap
end

function plot_LDOS(ax, dict; colormap = :thermal, zlims = (0, 1), kw...)
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    hmap = heatmap!(ax, xs, ys, data; colormap = colormap, colorrange = zlims, rasterize = 5, kw...)
    return hmap
end

function add_axis_inset(pos=fig[1, 1]; halign, valign, width=Relative(0.5),height=Relative(0.35), alignmode=Mixed(left=5, right=5), aspect = 1)
    inset_box = CairoMakie.Axis(pos; width, height, halign, valign, alignmode, aspect)
    # bring content upfront
    translate!(inset_box.scene, 0, 0, 10)
    return inset_box
end

function R(φ, harmonics;  R0 = 1)
    return (1 + sum([real(har * exp(n * φ * im)) for (n, har) in harmonics]))*R0
end

function to_cart(r, φ)
    return (r * cos(φ), r * sin(φ))
end

function plot_profile(pos, harmonics; color = :yellow, linestyle = :solid)
    inset_ax = add_axis_inset(pos; halign = :center, valign = :center, height = Relative(0.6), width = Relative(0.6))
    dict_har = Dict(ℓ => har for (ℓ, har) in enumerate(harmonics))
    lines!(inset_ax, [to_cart(R(φ, dict_har), φ) for φ in range(0, 2π, length = 1000)], linewidth = 1.5; color, linestyle)
    hidedecorations!(inset_ax)
    hidespines!(inset_ax)
end

## Figure

harmonics1 = harmonics_array(0.2, 1, 300)

function FigModeMixing(files, zmins, zmaxs, file_harmonics;)
    fig = Figure(resolution = (1100,  650), font = "CMU Serif Roman")
    xlabel = L"\Phi/\Phi_0"
    ylabel = L"$\omega$ (meV)"
    label = L"$$LDOS (arb. units)"
    xticks = Dict(
        1 => [0, 1, 2],
        2 => [0],
        3 => [1],
        4 => [2])
    xranges = Dict(
        1 => (-0.499, 2.499),
        2 => (-0.499, 0.499),
        3 => (0.501, 1.499),
        4 => (1.501, 2.499)
    )
    yranges = Dict(
        1 => (-0.26, 0.26),
        2 => (-0.05, 0.05),
        3 => (-0.05, 0.05),
        4 => (-0.05, 0.05)
    )

    xlines = [
        (-0.499, -0.2),
        (2.3, 2.499)
    ]

    yticks_t = [-0.2, 0, 0.2]
    yticks_z = [-0.02, 0, 0.02]

    harmonics = load(file_harmonics)

    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        yticks = col == 1 ? yticks_t : yticks_z
        ax = CairoMakie.Axis(fig[row, col]; xlabel, ylabel,  xticks = xticks[col], yticks, )
        xlims!(ax, xranges[col])
        ylims!(ax, yranges[col])
        ax.backgroundcolor = :black
        dict = load(file)
        zmin = zmins[row][col]
        zmax = zmaxs[row][col]
        zlims = (zmin, zmax)
        #hmap = plot_LDOS(ax, dict, zscale; zlims, lowclip = :black)
        plot_LDOS(ax, dict; zlims, lowclip = :black)

        #col == 1 && row == 2 && text!(ax, 0, 0; text = L"$\frac{\sigma}{\ell}$", color = :white, align = (:center, :center))
        #col == 1 && row == 3 && text!(ax, 0, 0; text = L"$\frac{\sigma}{\ell^2}$", color = :white, align = (:center, :center))

        col == 2 && row == 1 &&  plot_profile(fig[row, col], zeros(12))
      
        col == 2 && row == 2 &&  plot_profile(fig[row, col], harmonics["mod2"])
        col == 2 && row == 3 &&  plot_profile(fig[row, col], harmonics1)


        col == 2 && row == 1 && text!(ax, 0, 0; text = L"$\sigma_\ell = 0$", color = :white, align = (:center, :center), fontsize = 14)
        col == 2 && row == 2 && text!(ax, 0, 0; text = L"$\sigma_\ell = \frac{\sigma_1}{\ell^2}$", color = :white, align = (:center, :center), fontsize = 14)
        col == 2 && row == 3 && text!(ax, 0, 0; text = L"$\sigma_\ell = \frac{\sigma_1}{\ell}$", color = :white, align = (:center, :center), fontsize = 14)

        if (col == 3) & (row ==1)
            text!(ax, 1.1, 0.04; text = L"+4", color = :white, align = (:center, :center))
            text!(ax, 1.1, -0.04; text = L"-4", color = :white, align = (:center, :center))
            text!(ax, 0.75, 0.015; text = L"+3", color = :white, align = (:center, :center))
            text!(ax, 0.75, -0.015; text = L"-3", color = :white, align = (:center, :center))
            text!(ax, 0.95, 0.02; text = L"0", color = :white, align = (:center, :center))
            text!(ax, 0.95, -0.02; text = L"0", color = :white, align = (:center, :center))
        end

        if (col == 4) & (row ==1)
            text!(ax, 2.4, 0.01; text = L"+\frac{9}{2}", color = :white, align = (:center, :center))
            text!(ax, 2.4, -0.01; text = L"-\frac{9}{2}", color = :white, align = (:center, :center))
            text!(ax, 2.35, 0.035; text = L"+\frac{7}{2}", color = :white, align = (:center, :center))
            text!(ax, 2.35, -0.035; text = L"-\frac{7}{2}", color = :white, align = (:center, :center))
            text!(ax, 1.9, 0.015; text = L"+\frac{5}{2}", color = :white, align = (:center, :center))
            text!(ax, 1.9, -0.015; text = L"-\frac{5}{2}", color = :white, align = (:center, :center))
            text!(ax, 1.7, 0.03; text = L"+\frac{1}{2}", color = :white, align = (:center, :center))
            text!(ax, 1.7, -0.03; text = L"-\frac{1}{2}", color = :white, align = (:center, :center))
        
        end


        row != 3 && hidexdecorations!(ax, ticks = false)
        (col != 1) && (col != 2) && hideydecorations!(ax, ticks = false)
        col == 2 && hideydecorations!(ax, ticks = false, ticklabels = false)
   

        col == 4 && Colorbar(fig[row, 5], colormap = :thermal, limits = (0,1), label = label, ticklabelsvisible = true, ticks = [0,1], labelpadding = -10, width = 10, ticksize = 2, ticklabelpad = 5)
    end

    Label(fig[1, 2, Top()], L"$n=0$", padding = (0, 0, 5, 0))
    Label(fig[1, 3, Top()], L"$n=1$", padding = (0, 0, 5, 0))
    Label(fig[1, 4, Top()], L"$n=2$", padding = (0, 0, 5, 0))



    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (-30, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 3, TopLeft()], "c", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 4, TopLeft()], "d", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(fig[2, 1, TopLeft()], "e", padding = (-40, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[2, 2, TopLeft()], "f", padding = (-30, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[2, 3, TopLeft()], "g", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[2, 4, TopLeft()], "h", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)

    Label(fig[3, 1, TopLeft()], "i", padding = (-40, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[3, 2, TopLeft()], "j", padding = (-30, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[3, 3, TopLeft()], "k", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)
    Label(fig[3, 4, TopLeft()], "l", padding = (-10, 0, -25, 0), font = "CMU Serif Bold", fontsize = 20)

    colsize!(fig.layout, 1, Relative(0.5))
    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 10)
    colgap!(fig.layout, 3, 10)
    colgap!(fig.layout, 4, 5)

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

## Generate Figure

p = "Data/Fig_ModeMixing/"
try run(`mkdir Figures`) catch end


files = [
    ["none", "none_zoom0", "none_zoom1", "none_zoom2"],
    ["mod2_total", "mod2_zoom0", "mod2_zoom1", "mod2_zoom2"],
    ["mod1_total", "mod1_zoom0", "mod1_zoom1", "mod1_zoom2"],
]

files = map(file -> "$(p)".*file.*".jld2", files)
file_harmonics = "$(p)harmonics.jld2"

zmins = [
    [2e-3, 2e-4, 2e-4, 2e-4],
    [2e-4, 2e-4, 2e-4, 2e-4],
    [1e-4, 2e-4, 1e-3, 6e-4],
]

zmaxs = [
    [5e-2, 1e-2, 3e-2, 3e-2],
    [5e-2, 1e-2, 3e-2, 3e-2],
    [5e-2, 1e-2, 1e-2, 3e-2],
]

fig = FigModeMixing(files,zmins, zmaxs, file_harmonics;)
save("Figures/FigModeMixing.pdf", fig)
fig

