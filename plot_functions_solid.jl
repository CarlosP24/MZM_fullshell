using JLD2, CSV, DataFrames, CairoMakie, ColorSchemes, Colors, LinearAlgebra
my_dcolors = [:blue, :purple, :red]

## LDOS
function plot_LDOS(dict, ax, zscale; zlims = (0, 1), colorscale = identity, kw...)
    local hmap 
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    m = zscale * maximum(data)

    data ./= m
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = rast, colorscale, kw...)
    ax.backgroundcolor = :black
    return hmap
end

function plot_radial(dict, ax, zscale; kw...)
    local hmap
    data = sum(values(dict["mjdict"]))
    xs, ys = dict["xs"], real.(dict["ys"])
    m = zscale * maximum(data)
    data ./= m
    zlims = (0, 1)
    data ./= xs
    data .*= maximum(xs)
    data .= sqrt.(data)
    data = [reverse(data[:, 2:end], dims = 2) data]
    ys = [-reverse(ys[2:end]); ys]
    
    hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = rast, kw...)
    xlims!(ax, (0, 70))
    ax.backgroundcolor = :black
    return hmap
end

## Bands
cmap = Makie.ColorSchemes.colorschemes[:Set1_9]
ncolors = length(cmap)

function plot_b(dict, ax)
    Zs = sort(collect(keys(dict)))
    for Z in Zs
        color = cmap[mod1(Z+1, ncolors)]
        pts = dict[Z]
        foreach(p -> lines!(ax, p; color, label = p == first(pts) ? string(Z) : nothing, linewidth = 3), pts)
    end
end

## WF
function plot_wf_v2(dict, ax)
    WF_array =  vcat(collect(values(dict))...)
    ravs= map( x -> sum(x.psi .* x.rs), WF_array)
    bins = minimum(ravs):2e-1:maximum(ravs)
    c_Dict = Dict([bin => findall( x -> abs(x - bin) < 2e-1, ravs) for bin in bins])
    filter!(p -> isempty(p.second) ? false : true, c_Dict )
    color_Dict = Dict([ϵ => cmap[mod1.(unique(floor.(Int, map(x -> x.mJ, WF_array[c]))).+1, ncolors)] for (ϵ, c) in c_Dict])
    for (ϵ, c) in c_Dict
        rs = WF_array[c[1]].rs 
        psi = WF_array[c[1]].psi
        nrs = length(rs)
        colD = color_Dict[ϵ]
        colors = repeat(colD, 10)
        ncol = length(colors)
        l = ceil(Int, nrs/ncol) 
        n1 = 2-l
        n2 = 0
        for col in colors
            n1 = n1 + l  
            n2 = minimum([n1 + l, nrs])
            lines!(ax, rs[n1:n2], psi[n1:n2], color = col, linewidth = 3)
        end
    end
end

## Dome profile

dome_colors = Dict(:sml => colorant"#697a5c", :smb => colorant"#d1f5b9", :sc => colorant"#769dfb", :μ => colorant"#fd0e10")

function plot_dome(dict, ax)
    rs = dict["rs"]
    ϵs = dict["ϵs"]
    mes = minimum(ϵs)
    Mes = maximum(ϵs)
    lower = findall(x -> x < 0, ϵs)
    lrs = rs[lower]
    lϵs = ϵs[lower]

    band!(ax, lrs, lϵs, 0, color = (dome_colors[:smb], 0.5))
    band!(ax, [70, 80], mes - 10, 0, color = (dome_colors[:smb], 0.5))
    lines!(ax, rs, ϵs, color = dome_colors[:sml], linewidth = 3)
   
    lines!(ax, [0, 80], [0], color = :black, linestyle = :dash)
    lines!(ax, [70], [mes - 8, last(ϵs)], color = dome_colors[:sml], linewidth = 3)
    lines!(ax, [70], [mes, 0], color = dome_colors[:sml], linewidth = 3, linestyle = :dash)
    arrows!(ax, [70], [mes - 2], [0], [-6], linewidth = 3, color = dome_colors[:sml], arrowsize = 10 )

    ylims!(ax, mes - 10, Mes + 10)
    xlims!(ax, 0, 80)
    text!(ax, 20, mes + 5, text = "SM", color = :black, fontsize = 10)
    text!(ax, 72, mes + 5, text = "SC", color = :black, fontsize = 10 )
    
end

## Schematics

function schematics(ax)
    poly!(ax, Circle(Point2f(0,0), 10f0), color = colorant"#B9C8F7", strokewidth = 2)
    poly!(ax, Circle(Point2f(0,0), 8f0), color = colorant"#FBF5CD", strokewidth = 2)
    poly!(ax, Circle(Point2f(-8, 8), 0.5f0), color = :white, strokewidth = 2)
    poly!(ax, Circle(Point2f(-8, 8), 0.1f0), color = :black, )
    arrows!([0,  -8.3, -9.7], [0,  0, 0], [7.5,  -1.2, 1.2], [0,  0, 0], color = :black, linewidth = 3, arrowsize = 7)
    
    text!(ax, 3, 0.5; text = L"R")
    text!(ax, -9.7, 0.5; text = L"d") 
    text!(ax, -10, 7.2; text = L"z") 

    hidespines!(ax)
    hidedecorations!(ax)
    return
end

